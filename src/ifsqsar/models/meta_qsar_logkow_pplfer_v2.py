"""Meta QSAR for logKow (wet solvent)"""
import numpy as np
value_names = ('logKow',)
version = 2
endpoint = 'Log of wet (practical) octanol-water partition coefficient (log P) - updated for PFAS'
citation = 'Brown, T. N.; Armitage, J. M.; Sangion, A.; Arnot, J. A.; '\
           'Improved prediction of PFAS partitioning with PPLFERs and QSPRs. '\
           'Environ. Sci.: Process. Impacts, 2024, Accepted.'
round_digits = 2
units = 'log L[w]/L[wet o]'
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
    logKow = -1.219 * solutedependencies[0]['S'][0] \
             -0.058 * solutedependencies[0]['A'][0] \
             -3.579 * solutedependencies[0]['B'][0] \
             +2.702 * solutedependencies[0]['Vf'][0] \
             +0.341 * solutedependencies[0]['L'][0] \
             +0.326
    logKowUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKowUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKowUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKowUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKowUL += 3**2
    if ecount+ucount < 4:
        logKowUL = int(np.ceil((logKowUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKowUL = 5
    logKowerr = solutedependencies[0]['Vf'][0]**2 * (0.047**2) + \
                (solutedependencies[0]['L'][0] * (0.341))**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.012**2) / (0.341)**2) + \
                (0.025**2)
    if solutedependencies[0]['S'][0] != 0:
        logKowerr += (solutedependencies[0]['S'][0] * (-1.219))**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.035**2) / (-1.219)**2)
    if solutedependencies[0]['A'][0] != 0:
        logKowerr += (solutedependencies[0]['A'][0] * (-0.058))**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.028**2) / (-0.058)**2)
    if solutedependencies[0]['B'][0] != 0:
        logKowerr += (solutedependencies[0]['B'][0] * (-3.579))**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.034**2) / (-3.579)**2)
    logKowerr = logKowerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKowULconcat = []
            if ecount > 0:
                logKowULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKowULconcat.append('U')
                noteconcat.append('user')
            logKowULconcat.append(str(logKowUL))
            noteconcat.append('predicted values')
            if logKowUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKowUL = ''.join(logKowULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKowULconcat = []
            if ecount > 0:
                logKowULconcat.append('E')
            if ucount > 0:
                logKowULconcat.append('U')
            logKowUL = ''.join(logKowULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKowUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKowUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKowerr *= errorscale

    return round(logKow, round_digits), logKowUL, round(logKowerr, round_digits), ', '.join(domainnotes), citation, units, endpoint

