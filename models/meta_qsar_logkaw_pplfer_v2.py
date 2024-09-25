"""Meta QSAR for logKaw"""
import numpy as np
value_names = ('logKaw',)
version = 2
endpoint = 'Log of air-water partition coefficient (Henry\'s Law Constant) - updated for PFAS'
citation = 'Brown, T. N.; Armitage, J. M.; Sangion, A.; Arnot, J. A.; '\
           'Improved prediction of PFAS partitioning with PPLFERs and QSPRs. '\
           'Environ. Sci.: Process. Impacts, 2024, Accepted.'
round_digits = 2
units = 'log L[w]/L[a]'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = [('S', 3), ('A', 3), ('B', 3), ('Vf', 1), ('L', 3)]
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    logKaw = -2.127 * solutedependencies[0]['S'][0] \
             -3.690 * solutedependencies[0]['A'][0] \
             -4.783 * solutedependencies[0]['B'][0] \
             +2.505 * solutedependencies[0]['Vf'][0] \
             -0.445 * solutedependencies[0]['L'][0] \
             +0.504
    logKawUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKawUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKawUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKawUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKawUL += 3**2
    if ecount+ucount < 4:
        logKawUL = int(np.ceil((logKawUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKawUL = 5
    logKawerr = (solutedependencies[0]['Vf'][0] * 0.047)**2 + \
                (solutedependencies[0]['L'][0] * -0.445)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.013 / -0.445)**2) + \
                0.027**2
    if solutedependencies[0]['S'][0] != 0:
        logKawerr += (solutedependencies[0]['S'][0] * -2.127)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.038 / -2.127)**2)
    if solutedependencies[0]['A'][0] != 0:
        logKawerr += (solutedependencies[0]['A'][0] * -3.690)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.038 / -3.690)**2)
    if solutedependencies[0]['B'][0] != 0:
        logKawerr += (solutedependencies[0]['B'][0] * -4.783)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.037 / -4.783)**2)
    logKawerr = logKawerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKawULconcat = []
            if ecount > 0:
                logKawULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKawULconcat.append('U')
                noteconcat.append('user')
            logKawULconcat.append(str(logKawUL))
            noteconcat.append('predicted values')
            if logKawUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKawUL = ''.join(logKawULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKawULconcat = []
            if ecount > 0:
                logKawULconcat.append('E')
            if ucount > 0:
                logKawULconcat.append('U')
            logKawUL = ''.join(logKawULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKawUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKawUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKawerr *= errorscale

    return round(logKaw, round_digits), logKawUL, round(logKawerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

