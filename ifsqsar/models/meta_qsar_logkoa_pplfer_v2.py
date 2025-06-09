"""Meta QSAR for logKoa (dry solvent)"""
import numpy as np
value_names = ('logKoa',)
version = 2
endpoint = 'Log of octanol-air partition coefficient - updated for PFAS'
citation = 'Brown, T. N.; Armitage, J. M.; Sangion, A.; Arnot, J. A.; '\
           'Improved prediction of PFAS partitioning with PPLFERs and QSPRs. '\
           'Environ. Sci.: Process. Impacts, 2024, Accepted.'
round_digits = 2
units = 'log L[a]/L[o]'
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
    # recal. with pfas
    logKoa = +0.475 * solutedependencies[0]['S'][0] \
             +3.566 * solutedependencies[0]['A'][0] \
             +0.885 * solutedependencies[0]['B'][0] \
             +0.109 * solutedependencies[0]['Vf'][0] \
             +0.892 * solutedependencies[0]['L'][0] \
             -0.166
    logKoaUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKoaUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKoaUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKoaUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKoaUL += 3**2
    if ecount+ucount < 4:
        logKoaUL = int(np.ceil((logKoaUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKoaUL = 5
    logKoaerr = (solutedependencies[0]['Vf'][0] * 0.059)**2 + \
                (solutedependencies[0]['L'][0] * 0.892)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.016 / 0.892)**2) + \
                0.028**2
    if solutedependencies[0]['S'][0] != 0:
        logKoaerr += (solutedependencies[0]['S'][0] * 0.475)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.050 / 0.475)**2)
    if solutedependencies[0]['A'][0] != 0:
        logKoaerr += (solutedependencies[0]['A'][0] * 3.566)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.052 / 3.566)**2)
    if solutedependencies[0]['B'][0] != 0:
        logKoaerr += (solutedependencies[0]['B'][0] * 0.885)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.049 / 0.885)**2)
    logKoaerr = logKoaerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKoaULconcat = []
            if ecount > 0:
                logKoaULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKoaULconcat.append('U')
                noteconcat.append('user')
            logKoaULconcat.append(str(logKoaUL))
            noteconcat.append('predicted values')
            if logKoaUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKoaUL = ''.join(logKoaULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKoaULconcat = []
            if ecount > 0:
                logKoaULconcat.append('E')
            if ucount > 0:
                logKoaULconcat.append('U')
            logKoaUL = ''.join(logKoaULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKoaUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKoaUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKoaerr *= errorscale

    return round(logKoa, round_digits), logKoaUL, round(logKoaerr, round_digits), ', '.join(domainnotes), citation, units, endpoint

