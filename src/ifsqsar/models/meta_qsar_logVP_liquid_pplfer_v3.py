"""Meta QSAR for logVPliquid"""
import numpy as np
value_names = ('logVPliquid',)
version = 3
endpoint = 'Log of vapor pressure of liquid or super-cooled liquid predicted by PPLFER at 298K - updated for PFAS'
citation = 'Brown, T. N.; Armitage, J. M.; Sangion, A.; Arnot, J. A.; '\
           'Improved prediction of PFAS partitioning with PPLFERs and QSPRs. '\
           'Environ. Sci.: Process. Impacts, 2024, Accepted.'
round_digits = 2
units = 'log Pa'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = [('S', 3), ('A', 3), ('B', 3), ('Vf', 1), ('L', 3), ('state', 1), ('MVliquid', 2)]
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {'O': (round(np.log10(3169.0), round_digits), 0, np.nan, 'experimental value used', 'well known value', units, endpoint)}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    domainnotes = [propagated_domain_notes]
    # calculate AB and ABerr
    AB = (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])**0.5
    if solutedependencies[0]['A'][0] > 0 and solutedependencies[0]['B'][0] > 0:
        ABerr = solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0] * ((solutedependencies[0]['A'][2]/solutedependencies[0]['A'][0])**2 + (solutedependencies[0]['B'][2]/solutedependencies[0]['B'][0])**2)**0.5
        ABerr = AB * 0.5 * ABerr / (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])
    else:
        ABerr = 0
    # calculate logSa with pplfer equation from Brown et al. 2023
    logSa = -1.309 * solutedependencies[0]['S'][0] \
            -0.983 * solutedependencies[0]['A'][0] \
            -0.558 * solutedependencies[0]['B'][0] \
            -1.629 * AB \
            -0.779 * solutedependencies[0]['Vf'][0] \
            -0.655 * solutedependencies[0]['L'][0] \
            +0.690
    # calculate UL and err
    logSaUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logSaUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logSaUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logSaUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logSaUL += 3**2
    if ecount+ucount < 4:
        logSaUL = int(np.ceil((logSaUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logSaUL = 5
    logSaerr = (solutedependencies[0]['Vf'][0] * 0.132)**2 + \
               (solutedependencies[0]['L'][0] * -0.655)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.035 / -0.655)**2) + \
                0.068**2
    if solutedependencies[0]['S'][0] != 0:
        logSaerr += (solutedependencies[0]['S'][0] * -1.309)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.098 / -1.309)**2)
    if solutedependencies[0]['A'][0] != 0:
        logSaerr += (solutedependencies[0]['A'][0] * -0.983)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.217 / -0.983)**2)
    if solutedependencies[0]['B'][0] != 0:
        logSaerr += (solutedependencies[0]['B'][0] * -0.558)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.118 / -0.558)**2)
    if AB != 0:
        logSaerr += (AB * -1.629)**2 * ((ABerr / AB)**2 + (0.257 / -1.629)**2)
    logSaerr = logSaerr ** 0.5
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logSaULconcat = []
            if ecount > 0:
                logSaULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logSaULconcat.append('U')
                noteconcat.append('user')
            logSaULconcat.append(str(logSaUL))
            noteconcat.append('predicted values')
            if logSaUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logSaUL = ''.join(logSaULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logSaULconcat = []
            if ecount > 0:
                logSaULconcat.append('E')
            if ucount > 0:
                logSaULconcat.append('U')
            logSaUL = ''.join(logSaULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logSaUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    # apply error scaling
    if logSaUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
        if solutedependencies[0]['state'][3] in ['likely solid', 'maybe solid']:
            errorscale *= 5/3
        if logSaUL in [2, 3, 5]:
            errorscale *= 1.25
    logSaerr *= errorscale
    # convert to VP and cap at atmospheric pressure
    logVP = logSa + np.log10(8314.46261815324) + np.log10(298.15)
    if logVP > np.log10(101325):
        domainnotes = ['Predicted vapor pressure ({}) capped at atmospheric pressure, UL set to 6; original aggregate UL: {}'.format(round(logVP, round_digits), logSaUL)] + domainnotes
        logVP = np.log10(101325)
        logSaUL = 6
    return round(logVP, round_digits), logSaUL, round(logSaerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

