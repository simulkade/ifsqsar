"""
ifsqsar/ifsqsar.py
developed by Trevor N. Brown
Implements a python API for applying IFS QSARs and a simple GUI that built on top of the API
The primary elements of the API are the functions apply_qsars_to_molecule and apply_qsars_to_molecule_list
"""

from openbabel import openbabel as ob
import numpy as np
from . import smiles_norm
import re

chargedatom = re.compile('[\[].+?[-+][\]]')
mixturespec = re.compile('(\{.*?\})')

def apply_qsars_to_molecule(qsarlist,
                            smiles,  # SMILES as string
                            converter=None,  # OBConversion
                            values=('insmi',
                                    'normsmi',
                                    'sminote',
                                    'OBMol',
                                    'endpoint',
                                    'units',
                                    'qsarpred',
                                    'UL',
                                    'error',
                                    'ULnote',
                                    'citation'
                                    ),  # iterable of values to be output
                            outformat='rows',  # 'dict', 'columns', 'rows'
                            header=True,  # True or False
                            separator='\t',  # any string
                            endline='\n',  # any string
                            ):
    """Apply a list of QSARs to a molecule as a SMILES and return formatted output.

    Required Arguments:
        qsarlist -- list of QSAR objects obtained from get_qsar_list function of the models subpackage
        smiles -- structure of the molecule to be predicted as a SMILES string

    Optional Arguments:
        converter -- openbabel OBConversion instance, saves a little load time if passed
        values -- tuple of outputs to be returned, all are included by default:
            "insmi" -- input SMILES
            "normsmi" -- normalized SMILES
            "sminote" -- warnings or errors from SMILES normalization
            "OBMol" -- openbabel OBMol instance of normsmi, only for dict format
            "endpoint" -- a description of the endpoint the model predicts
            "units" -- units of the predicted value
            "qsarpred" -- predicted value
            "UL" -- Uncertainty Level (UL) assigned by applicability domain checks
            "error" -- estimated prediction uncertainty
            "ULnote" -- applicability domain warnings
            "citation" -- literature to cite for the predicted value
        outformat -- "rows" (default) or "columns" for formatted text output, or "dict" for a dict
        header -- include header line in formatted text output, default=True
        separator -- column separator for formatted text output, default="\\t" (tab)
        endline -- row separator for formatted text output, default="\\n" (newline)

    """
    # initialize results dict from values iterable
    result = {'SMILES success': True, 'QSAR list': []}
    for val in ('insmi', 'normsmi', 'sminote'):
        if val in values:
            result[val] = ''
    if 'OBMol' in values:
        result['OBMol'] = None
    # check inputs
    assert type(smiles) == str
    if converter is not None:
        assert type(converter) == ob.OBConversion
    # smiles is a single molecule so pass to smiles_norm
    countchems = [0, 0, 0]
    if re.search(mixturespec, smiles) is None:
        chemismixture = False
        # generate normalized OBMol
        molecule, normsmiles, conversionnote = smiles_norm.convertsmiles(smiles, converter)
        # check conversion results and output
        if normsmiles == '':
            result['SMILES success'] = False
            if 'insmi' in values:
                result['insmi'] = smiles
            if 'sminote' in values:
                result['sminote'] = conversionnote
        else:
            if 'insmi' in values:
                result['insmi'] = smiles
            if 'normsmi' in values:
                result['normsmi'] = normsmiles
            if 'sminote' in values:
                result['sminote'] = conversionnote
        # make sure that the smiles note does not contain any separators or endlines
        if 'sminote' in values:
            seplist = [',', ';', '|', '~']
            if separator in seplist:
                seplist.remove(separator)
                result['sminote'] = result['sminote'].replace(separator, seplist[0])
            if endline in seplist:
                seplist.remove(endline)
                result['sminote'] = result['sminote'].replace(endline, seplist[1])
        # save the OBMol if required
        if result['SMILES success'] and 'OBMol' in values:
            result['OBMol'] = molecule
        # create list of solutes and solvents for passing to QSARs
        solutelist = (molecule,)
        solventlist = tuple()
        componentlist = tuple()
        solutef = (('u', '1'),)
        solventf = tuple()
        componentf = tuple()
    # smiles is a mixture so split into solutes and solvents
    else:
        chemismixture = True
        # split the input into component type flags and smiles
        mixturesplit = re.split(mixturespec, smiles)
        while '' in mixturesplit:
            mixturesplit.remove('')
        # parse through split smiles input
        solutelist = []
        solventlist = []
        componentlist = []
        solutef = []
        solventf = []
        componentf = []
        nextissolutesmiles = False
        nextissolventsmiles = False
        nextiscomponentsmiles = False
        normsmiles = ''
        sminote = ''
        for ms in mixturesplit:
            if nextissolutesmiles:
                countchems[0] += 1
                # generate normalized OBMol for solutes
                nextissolutesmiles = False
                molecule, normsmi, conversionnote = smiles_norm.convertsmiles(ms, converter)
                if normsmi == '':
                    result['SMILES success'] = False
                solutelist.append(molecule)
                normsmiles = ''.join([normsmiles, '{{solute {}}}'.format(countchems[0]), normsmi])
                if conversionnote != '':
                    notelist = [sminote, ''.join(['Notes for Solute (', ms, '): ', conversionnote])]
                    if '' in notelist:
                        notelist.remove('')
                    sminote = ', '.join(notelist)
            elif nextissolventsmiles:
                countchems[1] += 1
                # generate normalized OBMol for solvents
                nextissolventsmiles = False
                molecule, normsmi, conversionnote = smiles_norm.convertsmiles(ms, converter)
                if normsmi == '':
                    result['SMILES success'] = False
                solventlist.append(molecule)
                normsmiles = ''.join([normsmiles, '{{solvent {}}}'.format(countchems[1]), normsmi])
                if conversionnote != '':
                    notelist = [sminote, ''.join(['Notes for Solvent (', ms, '): ', conversionnote])]
                    if '' in notelist:
                        notelist.remove('')
                    sminote = ', '.join(notelist)
            elif nextiscomponentsmiles:
                countchems[2] += 1
                # generate normalized OBMol for solvents
                nextiscomponentsmiles = False
                molecule, normsmi, conversionnote = smiles_norm.convertsmiles(ms, converter)
                if normsmi == '':
                    result['SMILES success'] = False
                componentlist.append(molecule)
                normsmiles = ''.join([normsmiles, '{{component {}}}'.format(countchems[2]), normsmi])
                if conversionnote != '':
                    notelist = [sminote, ''.join(['Notes for Component (', ms, '): ', conversionnote])]
                    if '' in notelist:
                        notelist.remove('')
                    sminote = ', '.join(notelist)
            else:
                splitms = ms.lstrip('{').rstrip('}').split(',')
                if splitms[0] == 'solute':
                    nextissolutesmiles = True
                    if len(splitms) > 1:
                        ftype, fvalue = splitms[1].split(':')
                        solutef.append((ftype, fvalue))
                    else:
                        solutef.append(('u', '0'))
                elif splitms[0] == 'solvent':
                    nextissolventsmiles = True
                    if len(splitms) > 1:
                        ftype, fvalue = splitms[1].split(':')
                        solventf.append((ftype, fvalue))
                    else:
                        solventf.append(('u', '1'))
                elif splitms[0] == 'component':
                    nextiscomponentsmiles = True
                    if len(splitms) > 1:
                        ftype, fvalue = splitms[1].split(':')
                        componentf.append((ftype, fvalue))
                    else:
                        componentf.append(('u', '1'))
                else:
                    notelist = [sminote, 'SMILES error: invalid component type specification']
                    if '' in notelist:
                        notelist.remove('')
                    sminote = ', '.join(notelist)
        # convert lists to tuples
        solutelist = tuple(solutelist)
        solventlist = tuple(solventlist)
        componentlist = tuple(componentlist)
        solutef = tuple(solutef)
        solventf = tuple(solventf)
        componentf = tuple(componentf)
        # check conversion results and output
        if not result['SMILES success']:
            if 'insmi' in values:
                result['insmi'] = smiles
            if 'sminote' in values:
                result['sminote'] = sminote
        else:
            if 'insmi' in values:
                result['insmi'] = smiles
            if 'normsmi' in values:
                result['normsmi'] = normsmiles
            if 'sminote' in values:
                result['sminote'] = sminote
        # make sure that the smiles note does not contain any separators or endlines
        if 'sminote' in values:
            seplist = [',', ';', '|', '~']
            if separator in seplist:
                seplist.remove(separator)
                result['sminote'] = result['sminote'].replace(separator, seplist[0])
            if endline in seplist:
                seplist.remove(endline)
                result['sminote'] = result['sminote'].replace(endline, seplist[1])
    # if chemical input is a mixture generate column headers
    qsarpredmixcolumns = []
    ULmixcolumns = []
    errormixcolumns = []
    if chemismixture:
        for mc in range(countchems[0]):
            qsarpredmixcolumns.append('qsarpred solute {}'.format(mc+1))
            ULmixcolumns.append('UL solute {}'.format(mc+1))
            errormixcolumns.append('error solute {}'.format(mc+1))
        for mc in range(countchems[2]):
            qsarpredmixcolumns.append('qsarpred component {}'.format(mc+1))
            ULmixcolumns.append('UL component {}'.format(mc+1))
            errormixcolumns.append('error component {}'.format(mc+1))
    # parse through the list of QSARs applying each to the molecule
    for qsar in qsarlist:
        # load the model
        qsar.load()
        # initialize dict of calculated results
        result['QSAR list'].append(qsar.model_name)
        result[qsar.model_name] = {}
        if 'endpoint' in values:
            result[qsar.model_name]['endpoint'] = ''
        if 'units' in values:
            result[qsar.model_name]['units'] = ''
        if 'qsarpred' in values:
            if qsar.ismixture and chemismixture:
                for mc in qsarpredmixcolumns:
                    result[qsar.model_name][mc] = np.nan
            else:
                result[qsar.model_name]['qsarpred'] = np.nan
        if 'UL' in values:
            if qsar.ismixture and chemismixture:
                for mc in ULmixcolumns:
                    result[qsar.model_name][mc] = np.nan
            else:
                result[qsar.model_name]['UL'] = np.nan
        if 'error' in values:
            if qsar.ismixture and chemismixture:
                for mc in errormixcolumns:
                    result[qsar.model_name][mc] = np.nan
            else:
                result[qsar.model_name]['error'] = np.nan
        if 'ULnote' in values:
            result[qsar.model_name]['ULnote'] = ''
        if 'citation' in values:
            result[qsar.model_name]['citation'] = ''
        # continue if SMILES was not successfully converted
        if not result['SMILES success']:
            continue
        # apply model and store output
        qsar_prediction, uncertainty_level, error, note, citation, units, endpoint = qsar.apply_model(solutes=solutelist, solvents=solventlist, components=componentlist, solutef=solutef, solventf=solventf, componentf=componentf)
        if 'endpoint' in values:
            result[qsar.model_name]['endpoint'] = endpoint
        if 'units' in values:
            result[qsar.model_name]['units'] = units
        if 'qsarpred' in values:
            if qsar.ismixture and chemismixture:
                if type(qsar_prediction) == list:
                    assert len(qsarpredmixcolumns) == len(qsar_prediction)
                    for mc in range(len(qsarpredmixcolumns)):
                        result[qsar.model_name][qsarpredmixcolumns[mc]] = qsar_prediction[mc]
                else:
                    for mc in range(len(qsarpredmixcolumns)):
                        result[qsar.model_name][qsarpredmixcolumns[mc]] = qsar_prediction
            else:
                result[qsar.model_name]['qsarpred'] = qsar_prediction
        if 'UL' in values:
            if qsar.ismixture and chemismixture:
                if type(uncertainty_level) == list:
                    assert len(ULmixcolumns) == len(uncertainty_level)
                    for mc in range(len(ULmixcolumns)):
                        result[qsar.model_name][ULmixcolumns[mc]] = uncertainty_level[mc]
                else:
                    for mc in range(len(ULmixcolumns)):
                        result[qsar.model_name][ULmixcolumns[mc]] = uncertainty_level
            else:
                result[qsar.model_name]['UL'] = uncertainty_level
        if 'error' in values:
            if qsar.ismixture and chemismixture:
                if type(error) == list:
                    assert len(errormixcolumns) == len(error)
                    for mc in range(len(errormixcolumns)):
                        result[qsar.model_name][errormixcolumns[mc]] = error[mc]
                else:
                    for mc in range(len(ULmixcolumns)):
                        result[qsar.model_name][ULmixcolumns[mc]] = error
            else:
                result[qsar.model_name]['error'] = error
        if 'ULnote' in values:
            result[qsar.model_name]['ULnote'] = note
            # make sure that the note does not contain any separators or endlines
            seplist = [',', ';', '|', '~']
            if separator in seplist:
                seplist.remove(separator)
                result[qsar.model_name]['ULnote'] = result[qsar.model_name]['ULnote'].replace(separator, seplist[0])
            if endline in seplist:
                seplist.remove(endline)
                result[qsar.model_name]['ULnote'] = result[qsar.model_name]['ULnote'].replace(endline, seplist[1])
        if 'citation' in values:
            result[qsar.model_name]['citation'] = citation
            # make sure that the note does not contain any separators or endlines
            seplist = [',', ';', '|', '~']
            if separator in seplist:
                seplist.remove(separator)
                result[qsar.model_name]['citation'] = result[qsar.model_name]['citation'].replace(separator, seplist[0])
            if endline in seplist:
                seplist.remove(endline)
                result[qsar.model_name]['citation'] = result[qsar.model_name]['citation'].replace(endline, seplist[1])
    # return output as dict of values
    if outformat == 'dict':
        return result
    # return output as a string where the output values are in rows and the chemical is in the column
    elif outformat == 'columns':
        outstring = ''
        for val in values:
            if val in ('insmi', 'normsmi', 'sminote'):
                if header:
                    outstring = ''.join([outstring, val, separator, result[val], endline])
                else:
                    outstring = ''.join([outstring, result[val], endline])
        for qsar in result['QSAR list']:
            for val in values:
                if val in ('endpoint', 'units', 'ULnote', 'citation'):
                    if header:
                        outstring = ''.join([outstring, qsar, ' ', val, separator, str(result[qsar][val]), endline])
                    else:
                        outstring = ''.join([outstring, str(result[qsar][val]), endline])
                if val in ('qsarpred', 'UL', 'error'):
                    if val in result[qsar]:
                        if header:
                            if type(result[qsar][val]) != str and np.isnan(result[qsar][val]):
                                outstring = ''.join([outstring, qsar, ' ', val, separator, '', endline])
                            else:
                                outstring = ''.join([outstring, qsar, ' ', val, separator, str(result[qsar][val]), endline])
                        else:
                            if type(result[qsar][val]) != str and np.isnan(result[qsar][val]):
                                outstring = ''.join([outstring, '', endline])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][val]), endline])
                    for s in range(1, 201):
                        localval = ' '.join([val, 'solute', str(s)])
                        if localval not in result[qsar]:
                            break
                        if header:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, qsar, ' ', localval, separator, '', endline])
                            else:
                                outstring = ''.join([outstring, qsar, ' ', localval, separator, str(result[qsar][localval]), endline])
                        else:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, '', endline])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][localval]), endline])
                    for c in range(1, 201):
                        localval = ' '.join([val, 'component', str(c)])
                        if localval not in result[qsar]:
                            break
                        if header:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, qsar, ' ', localval, separator, '', endline])
                            else:
                                outstring = ''.join([outstring, qsar, ' ', localval, separator, str(result[qsar][localval]), endline])
                        else:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, '', endline])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][localval]), endline])
        return outstring
    # return output as a string where the output values are in columns and the chemical is in the row
    elif outformat == 'rows':
        outstring = ''
        # header
        if header:
            first = True
            for val in values:
                if val in ('insmi', 'normsmi', 'sminote'):
                    if first:
                        outstring = ''.join([outstring, val])
                        first = False
                    else:
                        outstring = ''.join([outstring, separator, val])
            for qsar in result['QSAR list']:
                for val in values:
                    if val in ('endpoint', 'units', 'ULnote', 'citation'):
                        if first:
                            outstring = ''.join([outstring, qsar, ' ', val])
                            first = False
                        else:
                            outstring = ''.join([outstring, separator, qsar, ' ', val])
                    if val in ('qsarpred', 'UL', 'error'):
                        if val in result[qsar]:
                            if first:
                                outstring = ''.join([outstring, qsar, ' ', val])
                                first = False
                            else:
                                outstring = ''.join([outstring, separator, qsar, ' ', val])
                        for s in range(1, 201):
                            localval = ' '.join([val, 'solute', str(s)])
                            if localval not in result[qsar]:
                                break
                            if first:
                                outstring = ''.join([outstring, qsar, ' ', localval])
                                first = False
                            else:
                                outstring = ''.join([outstring, separator, qsar, ' ', localval])
                        for c in range(1, 201):
                            localval = ' '.join([val, 'component', str(c)])
                            if localval not in result[qsar]:
                                break
                            if first:
                                outstring = ''.join([outstring, qsar, ' ', localval])
                                first = False
                            else:
                                outstring = ''.join([outstring, separator, qsar, ' ', localval])
            outstring = ''.join([outstring, endline])
        # output values
        first = True
        for val in values:
            if val in ('insmi', 'normsmi', 'sminote'):
                if first:
                    outstring = ''.join([outstring, result[val]])
                    first = False
                else:
                    outstring = ''.join([outstring, separator, result[val]])
        for qsar in result['QSAR list']:
            for val in values:
                if val in ('endpoint', 'units', 'ULnote', 'citation'):
                    if first:
                        outstring = ''.join([outstring, str(result[qsar][val])])
                        first = False
                    else:
                        outstring = ''.join([outstring, separator, str(result[qsar][val])])
                if val in ('qsarpred', 'UL', 'error'):
                    if val in result[qsar]:
                        if first:
                            if type(result[qsar][val]) != str and np.isnan(result[qsar][val]):
                                outstring = ''.join([outstring, ''])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][val])])
                            first = False
                        else:
                            if type(result[qsar][val]) != str and np.isnan(result[qsar][val]):
                                outstring = ''.join([outstring, separator, ''])
                            else:
                                outstring = ''.join([outstring, separator, str(result[qsar][val])])
                    for s in range(1, 201):
                        localval = ' '.join([val, 'solute', str(s)])
                        if localval not in result[qsar]:
                            break
                        if first:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, ''])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][localval])])
                            first = False
                        else:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, separator, ''])
                            else:
                                outstring = ''.join([outstring, separator, str(result[qsar][localval])])
                    for c in range(1, 201):
                        localval = ' '.join([val, 'component', str(c)])
                        if localval not in result[qsar]:
                            break
                        if first:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, ''])
                            else:
                                outstring = ''.join([outstring, str(result[qsar][localval])])
                            first = False
                        else:
                            if type(result[qsar][localval]) != str and np.isnan(result[qsar][localval]):
                                outstring = ''.join([outstring, separator, ''])
                            else:
                                outstring = ''.join([outstring, separator, str(result[qsar][localval])])
        outstring = ''.join([outstring, endline])
        return outstring


def apply_qsars_to_molecule_list(qsarlist,
                                 smileslist=None,  # list of SMILES as strings
                                 infilename=None,  # input file name
                                 inheaderrows=1,  # number of header lines
                                 inheadtrgtrow=1,  # header row to select from
                                 inheadersmiles='smiles',  # header value indicating SMILES
                                 inseparator='\t',  # any string
                                 inendline='\n',  # any string
                                 converter=None,  # OBConversion
                                 values=('insmi',
                                         'normsmi',
                                         'sminote',
                                         'OBMol',
                                         'endpoint',
                                         'units',
                                         'qsarpred',
                                         'UL',
                                         'error',
                                         'ULnote',
                                         'citation'
                                         ),  # iterable of values to be output
                                 outfilename=None,  # output file name
                                 outkeepdata=True,  # also output all of the input file contents
                                 outformat='rows',  # 'dict', 'columns', 'rows'
                                 outheader=True,  # True or False
                                 outseparator='\t',  # any string
                                 outendline='\n',  # any string
                                 ):
    """Apply a list of QSARs to a list of SMILES, pass them to apply_qsars_to_molecule, then concatenate the results.

    Required Arguments:
        qsarlist -- list of QSARs obtained from get_qsar_list function of the models subpackage
                    pass None to get the default list

    Optional Arguments:
        smileslist -- a list of SMILES, either smileslist or infilename must be specified
        infilename -- name of delimited input file containing SMILES
        inheaderrows -- number of lines of headers at the start of the input file, default=1
        inheadtrgtrow -- 1-indexed header row to use, default=1
        inheadersmiles -- name of column in inheadtrgtrow that contains SMILES, default='smiles'
                          (not case sensitive)
        inseparator -- column separator for input file, default="\\t" (tab)
        inendline -- row separator for input file, default="\\n" (newline)
        converter -- openbabel OBConversion instance, saves a little load time if passed
        values -- tuple of outputs to be returned, all are included by default:
            "insmi" -- input SMILES
            "normsmi" -- normalized SMILES
            "sminote" -- warnings or errors from SMILES normalization
            "OBMol" -- openbabel OBMol instance of normsmi, only for dict format
            "endpoint" -- a description of the endpoint the model predicts
            "units" -- units of the predicted value
            "qsarpred" -- predicted value
            "UL" -- Uncertainty Level (UL) assigned by applicability domain checks
            "error" -- estimated prediction uncertainty
            "ULnote" -- applicability domain warnings
            "citation" -- literature to cite for the predicted value
        outfilename -- output file name, default=None which returns concatenated output
        outkeepdata -- include all of the input file contents in formatted text output, default=True
        outformat -- "rows" (default) or "columns" for formatted text output, or "dict" for a dict
        outheader -- include header line in formatted text output, default=True
        outseparator -- column separator for formatted text output, default="\\t" (tab)
        outendline -- row separator for formatted text output, default="\\n" (newline)
        """
    # load data from file
    filelines = None
    if infilename is not None:
        filetext = ''
        with open(infilename, 'r') as infile:
            filetext = infile.read()
        filelines = filetext.split(inendline)
        # remove empty lines
        while '' in filelines:
            filelines.remove('')
        # find the index of the column with SMILES
        smiles_index = 0
        splitline = filelines[inheadtrgtrow-1].split(inseparator)
        for s in range(len(splitline)):
            if splitline[s].lower().strip().rstrip() == inheadersmiles.lower():
                smiles_index = s
                break
        # extract SMILES into smileslist
        smileslist = []
        for i in range(inheaderrows, len(filelines)):
            splitline = filelines[i].split(inseparator)
            smileslist.append(splitline[smiles_index])
    # instantiate a converter if needed
    if converter is None:
        converter = ob.OBConversion()
        converter.SetInAndOutFormats('smi', 'can')
    # initialize dict to store output
    result = {'QSAR list':[]}
    for val in ('insmi', 'normsmi', 'sminote', 'OBMol'):
        if val in values:
            result[val] = []
    for qsar in qsarlist:
        # initialize dict of calculated results
        result['QSAR list'].append(qsar.model_name)
        result[qsar.model_name] = {}
        for val in values:
            if val in ('endpoint', 'units', 'ULnote', 'citation'):
                result[qsar.model_name][val] = []
    # parse through structures
    orderedcolumnlist = []
    for smiles in range(len(smileslist)):
        # apply qsars to this structure
        singleresult = apply_qsars_to_molecule(qsarlist,
                                               smileslist[smiles],
                                               converter=converter,
                                               values=values,
                                               outformat='dict',
                                               separator=outseparator,
                                               endline=outendline,
                                               )
        # concatenate to output dict
        for val in values:
            if val in ('insmi', 'normsmi', 'sminote', 'OBMol'):
                if val not in orderedcolumnlist:
                    orderedcolumnlist.append(val)
                result[val].append(singleresult[val])
        for qsar in result['QSAR list']:
            for val in values:
                if val in ('endpoint', 'units', 'ULnote', 'citation'):
                    if (qsar, val) not in orderedcolumnlist:
                        orderedcolumnlist.append((qsar, val))
                    elif val == 'citation' and (qsar, val) in orderedcolumnlist:
                        orderedcolumnlist.remove((qsar, val))
                        orderedcolumnlist.append((qsar, val))
                    result[qsar][val].append(singleresult[qsar][val])
                if val in ('qsarpred', 'UL', 'error'):
                    if val in singleresult[qsar]:
                        if val not in result[qsar]:
                            orderedcolumnlist.append((qsar, val))
                            result[qsar][val] = []
                            for i in range(smiles):
                                result[qsar][val].append(np.nan)
                        result[qsar][val].append(singleresult[qsar][val])
                    for s in range(1, 201):
                        localval = ' '.join([val, 'solute', str(s)])
                        if localval not in singleresult[qsar]:
                            break
                        if localval not in result[qsar]:
                            orderedcolumnlist.append((qsar, localval))
                            result[qsar][localval] = []
                            for i in range(smiles):
                                result[qsar][localval].append(np.nan)
                        result[qsar][localval].append(singleresult[qsar][localval])
                    for c in range(1, 201):
                        localval = ' '.join([val, 'component', str(c)])
                        if localval not in singleresult[qsar]:
                            break
                        if localval not in result[qsar]:
                            orderedcolumnlist.append((qsar, localval))
                            result[qsar][localval] = []
                            for i in range(smiles):
                                result[qsar][localval].append(np.nan)
                        result[qsar][localval].append(singleresult[qsar][localval])
        # top up any result lists that aren't full
        for qsar in result['QSAR list']:
            for val in ('qsarpred', 'UL', 'error'):
                if val in values:
                    if val in result[qsar]:
                        for i in range(1+smiles-len(result[qsar][val])):
                            result[qsar][val].append(np.nan)
                    for s in range(1, 201):
                        localval = ' '.join([val, 'solute', str(s)])
                        if localval not in result[qsar]:
                            break
                        for i in range(1+smiles-len(result[qsar][localval])):
                            result[qsar][localval].append(np.nan)
                    for c in range(1, 201):
                        localval = ' '.join([val, 'component', str(c)])
                        if localval not in result[qsar]:
                            break
                        for i in range(1+smiles-len(result[qsar][localval])):
                            result[qsar][localval].append(np.nan)
    # return results dict if outformat is dict
    outline = ''
    if outformat == 'dict':
        return result
    # create column output
    elif outformat == 'columns':
        for column in orderedcolumnlist:
            if column == 'OBMol':
                continue
            if type(column) == str:
                if outheader:
                    outline = ''.join([outline, column])
                for chem in result[column]:
                    if type(chem) != str and np.isnan(chem):
                        outline = outseparator.join([outline, ''])
                    else:
                        outline = outseparator.join([outline, str(chem)])
                outline = ''.join([outline, outendline])
            elif type(column) == tuple:
                if outheader:
                    outline = ''.join([outline, ' '.join([column[0], column[1]])])
                for chem in result[column[0]][column[1]]:
                    if type(chem) != str and np.isnan(chem):
                        outline = outseparator.join([outline, ''])
                    else:
                        outline = outseparator.join([outline, str(chem)])
                outline = ''.join([outline, outendline])
    # create rows output
    elif outformat == 'rows':
        if outheader:
            first = True
            for column in orderedcolumnlist:
                if column == 'OBMol':
                    continue
                if type(column) == str:
                    if first:
                        outline = column
                        first = False
                    else:
                        outline = outseparator.join([outline, column])
                elif type(column) == tuple:
                    if first:
                        outline = ' '.join([column[0], column[1]])
                        first = False
                    else:
                        outline = outseparator.join([outline, ' '.join([column[0], column[1]])])
            outline = ''.join([outline, outendline])
        for chem in range(len(smileslist)):
            first = True
            for column in orderedcolumnlist:
                if column == 'OBMol':
                    continue
                if type(column) == str:
                    if first:
                        outline = ''.join([outline, str(result[column][chem])])
                        first = False
                    else:
                        outline = outseparator.join([outline, str(result[column][chem])])
                elif type(column) == tuple:
                    if first:
                        if type(result[column[0]][column[1]][chem]) != str and np.isnan(result[column[0]][column[1]][chem]):
                            outline = ''.join([outline, ''])
                        else:
                            outline = ''.join([outline, str(result[column[0]][column[1]][chem])])
                        first = False
                    else:
                        if type(result[column[0]][column[1]][chem]) != str and np.isnan(result[column[0]][column[1]][chem]):
                            outline = outseparator.join([outline, ''])
                        else:
                            outline = outseparator.join([outline, str(result[column[0]][column[1]][chem])])
            outline = ''.join([outline, outendline])
    # if not outputting to file return result
    if outfilename is None:
        if outformat == 'dict:':
            return result
        else:
            return outline
    else:
        assert outformat == 'rows'
        # replace tokens for empty delimiters
        if outseparator == '<nosep>':
            outline = outline.replace('<nosep>', '')
        if outendline == '<noend>':
            outline = outline.replace('<noend>', '')
        # output with data from input file
        if outkeepdata and filelines is not None:
            resultlines = outline.split(outendline)
            # first append header line if kept
            infileoffset = inheaderrows
            resultoffset = 0
            if outheader:
                resultoffset = 1
                # append output header row to target row from input file
                filelines[inheadtrgtrow-1] = outseparator.join([filelines[inheadtrgtrow-1], resultlines[0]])
                # if there are more input header rows add empty output lines
                for i in range(infileoffset):
                    if i == inheadtrgtrow-1:
                        continue
                    filelines[i] = outseparator.join([filelines[i], outseparator * resultlines[0].count(outseparator)])
            for i in range(infileoffset, len(filelines)):
                filelines[i] = outseparator.join([filelines[i], resultlines[i-infileoffset+resultoffset]])
            with open(outfilename, 'w') as outfile:
                outfile.write(outendline.join(filelines))
        # output result without input data
        else:
            with open(outfilename, 'w') as outfile:
                outfile.write(outline)


class IFSGUIClass:
    """A GUI interface for reading in structures as SMILES and applying QSARs to the structures."""

    tk = __import__('tkinter')

    class _ReadOnlyText(tk.Text):
        """Subclass of tk.Text that is read-only."""

        def __init__(self, *args, **kwargs):
            """Replace insert and delete bindings."""
            # subclass to tk.Text
            import tkinter as tk
            tk.Text.__init__(self, *args, **kwargs)
            self.SEL = tk.SEL
            self.END = tk.END
            self.INSERT = tk.INSERT
            from idlelib.redirector import WidgetRedirector
            self.redirector = WidgetRedirector(self)
            # freeze user changes
            self.insert = self.redirector.register('insert', lambda *args, **kw: 'break')
            self.delete = self.redirector.register('delete', lambda *args, **kw: 'break')
            # bind ctrl-a as select all
            self.bind("<Control-Key-a>", self.select_all)

        def select_all(self, event):
            """Select all event bound to ctrl-a."""
            # select all text
            self.tag_add(self.SEL, '1.0', self.END)
            self.mark_set(self.INSERT, '1.0')
            self.see(self.INSERT)
            return 'break'

    def __init__(self):
        """GUI enters single mode by default. Available QSARs are loaded."""
        # initiate root
        self.root = self.tk.Tk()
        self.root.wm_title('IFSQSAR')
        self.filedialog = __import__('tkinter.filedialog', fromlist=[''])
        # create frame
        self.frame = self.tk.Frame(self.root)
        self.frame.pack_propagate(0)
        self.frame.pack()
        # initiate mode variable and import models
        self.mixturemode = self.tk.StringVar()
        self.mixturemode.set('purechemical')
        from . import models
        self.pure_qsarmodels = models.get_qsar_list(qsarlist=['fhlb', 'hhlb', 'hhlt', 'HLbiodeg',
                                                              'dsm', 'tmconsensus', 'tbpplfer',
                                                              'logKow', 'logKowdry', 'logKoa', 'logKaw', 'logKoo',
                                                              'logVPliquid', 'logSwliquid', 'logSoliquid',
                                                              'MVliquid', 'densityliquid', 'MVsolid', 'densitysolid', 'MW',
                                                              'state',
                                                              'E', 'S', 'A', 'B', 'V', 'Vf', 'L',
                                                              's', 'a', 'b', 'v', 'l', 'c'])
        self.mixture_qsarmodels = models.get_qsar_list(qsarlist=['logKsa'])
        # setup openbabel converter
        self.obcon = ob.OBConversion()
        self.obcon.SetInAndOutFormats('smi', 'can')
        # set to single mode
        self.setup_single_mode()
        # start up gui
        self.root.mainloop()

    def setup_single_mode(self):
        """Delete any batch mode widgets present and load single mode widgets."""
        # delete batch mode widgets
        if hasattr(self, 'buttoninput'):
            self.buttoninput.destroy()
            delattr(self, 'buttoninput')
        if hasattr(self, 'textinput'):
            self.textinput.destroy()
            delattr(self, 'textinput')
        if hasattr(self, 'framebatchradio'):
            self.framebatchradio.destroy()
            delattr(self, 'framebatchradio')
        if hasattr(self, 'modebatchradiopure'):
            self.modebatchradiopure.destroy()
            delattr(self, 'modebatchradiopure')
        if hasattr(self, 'modebatchradiomixt'):
            self.modebatchradiomixt.destroy()
            delattr(self, 'modebatchradiomixt')
        if hasattr(self, 'buttonoutput'):
            self.buttonoutput.destroy()
            delattr(self, 'buttonoutput')
        if hasattr(self, 'textoutput'):
            self.textoutput.destroy()
            delattr(self, 'textoutput')
        if hasattr(self, 'buttongotosingle'):
            self.buttongotosingle.destroy()
            delattr(self, 'buttongotosingle')
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.destroy()
            delattr(self, 'buttoninfo')
        if hasattr(self, 'framebatch'):
            self.framebatch.destroy()
            delattr(self, 'framebatch')
        if hasattr(self, 'buttoncalcbatch'):
            self.buttoncalcbatch.destroy()
            delattr(self, 'buttoncalcbatch')
        # single mode - text box to enter smiles
        self.labelsmiles = self.tk.Label(self.frame, text='Enter a SMILES', font='TkDefaultFont 10')
        self.labelsmiles.grid(row=0, column=0)
        self.entrysmiles = self.tk.Entry(self.frame, width=50, font='TkTextFont 10')
        self.entrysmiles.grid(row=0, column=1)
        # single mode - radio buttons to select input type
        self.framesingleradio = self.tk.Frame(self.frame)
        self.framesingleradio.grid(row=1, column=1)
        self.modesingleradiopure = self.tk.Radiobutton(self.framesingleradio, text='pure chemical', font='TkDefaultFont 10', variable=self.mixturemode, value='purechemical')
        self.modesingleradiopure.grid(row=0, column=0)
        self.modesingleradiomixt = self.tk.Radiobutton(self.framesingleradio, text='solute-solvent/mixture', font='TkDefaultFont 10', variable=self.mixturemode, value='mixture')
        self.modesingleradiomixt.grid(row=0, column=1)
        # single mode - text box to show results
        self.labelresult = self.tk.Label(self.frame, text='Model Results', font='TkDefaultFont 10')
        self.labelresult.grid(row=2, column=0)
        self.frametext = self.tk.Frame(self.frame)
        self.frametext.grid(row=2, column=1)
        self.scrbar = self.tk.Scrollbar(self.frametext)
        self.scrbar.pack(side=self.tk.RIGHT, fill=self.tk.Y)
        self.textresult = IFSGUIClass._ReadOnlyText(self.frametext, height=5, width=48, font='TkTextFont 10')
        self.textresult.pack(side=self.tk.LEFT)
        self.scrbar.config(command=self.textresult.yview)
        self.textresult.config(yscrollcommand=self.scrbar.set)
        # single mode - buttons to switch to batch mode and display info
        self.framesingle = self.tk.Frame(self.frame)
        self.framesingle.grid(row=3, column=0)
        self.buttongotobatch = self.tk.Button(self.framesingle, text='Batch Mode', font='TkDefaultFont 10',
                                         height=1, width=15, command=self.setup_batch_mode)
        self.buttongotobatch.grid(row=0, column=0)
        self.buttoninfo = self.tk.Button(self.framesingle, text='Info', font='TkDefaultFont 10',
                                    height=1, width=15, command=self.info)
        self.buttoninfo.grid(row=1, column=0)
        # single mode - button to calculate results
        self.buttoncalcsingle = self.tk.Button(self.frame, text='Apply IFS QSARs', font='TkDefaultFont 10',
                                          height=2, width=25, command=self.calculate_single)
        self.buttoncalcsingle.grid(row=3, column=1)

    def setup_batch_mode(self):
        """Delete any single mode widgets present and load batch mode widgets."""
        # delete single mode widgets
        if hasattr(self, 'labelsmiles'):
            self.labelsmiles.destroy()
            delattr(self, 'labelsmiles')
        if hasattr(self, 'entrysmiles'):
            self.entrysmiles.destroy()
            delattr(self, 'entrysmiles')
        if hasattr(self, 'framesingleradio'):
            self.framesingleradio.destroy()
            delattr(self, 'framesingleradio')
        if hasattr(self, 'modesingleradiopure'):
            self.modesingleradiopure.destroy()
            delattr(self, 'modesingleradiopure')
        if hasattr(self, 'modesingleradiomixt'):
            self.modesingleradiomixt.destroy()
            delattr(self, 'modesingleradiomixt')
        if hasattr(self, 'labelresult'):
            self.labelresult.destroy()
            delattr(self, 'labelresult')
        if hasattr(self, 'textresult'):
            self.textresult.destroy()
            delattr(self, 'textresult')
        if hasattr(self, 'scrbar'):
            self.scrbar.destroy()
            delattr(self, 'scrbar')
        if hasattr(self, 'frametext'):
            self.frametext.destroy()
            delattr(self, 'frametext')
        if hasattr(self, 'buttongotobatch'):
            self.buttongotobatch.destroy()
            delattr(self, 'buttongotobatch')
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.destroy()
            delattr(self, 'buttoninfo')
        if hasattr(self, 'framesingle'):
            self.framesingle.destroy()
            delattr(self, 'framesingle')
        if hasattr(self, 'buttoncalcsingle'):
            self.buttoncalcsingle.destroy()
            delattr(self, 'buttoncalcsingle')
        # clear filenames if they exist
        self.inputfilename = ''
        self.outputfilename = ''
        # batch mode - button to select input file
        self.buttoninput = self.tk.Button(self.frame, text='Select Input File', font='TkDefaultFont 10',
                                     height=1, width=15, command=self.select_input_file)
        self.buttoninput.grid(row=0, column=0)
        self.textinput = IFSGUIClass._ReadOnlyText(self.frame, height=1, width=50, font='TkTextFont 10')
        self.textinput.grid(row=0, column=1)
        # batch mode - radio buttons to select input type
        self.framebatchradio = self.tk.Frame(self.frame)
        self.framebatchradio.grid(row=1, column=1)
        self.modebatchradiopure = self.tk.Radiobutton(self.framebatchradio, text='pure chemical', font='TkDefaultFont 10', variable=self.mixturemode, value='purechemical')
        self.modebatchradiopure.grid(row=0, column=0)
        self.modebatchradiomixt = self.tk.Radiobutton(self.framebatchradio, text='solute-solvent/mixture', font='TkDefaultFont 10', variable=self.mixturemode, value='mixture')
        self.modebatchradiomixt.grid(row=0, column=1)
        # batch mode - button to select output file
        self.buttonoutput = self.tk.Button(self.frame, text='Select Ouput File', font='TkDefaultFont 10',
                                      height=1, width=15, command=self.select_output_file)
        self.buttonoutput.grid(row=2, column=0)
        self.textoutput = IFSGUIClass._ReadOnlyText(self.frame, height=1, width=50, font='TkTextFont 10')
        self.textoutput.grid(row=2, column=1)
        # batch mode - buttons to switch to single mode and display info
        self.framebatch = self.tk.Frame(self.frame)
        self.framebatch.grid(row=3, column=0)
        self.buttongotosingle = self.tk.Button(self.framebatch, text='Single Mode', font='TkDefaultFont 10',
                                          height=1, width=15, command=self.setup_single_mode)
        self.buttongotosingle.grid(row=0, column=0)
        self.buttoninfo = self.tk.Button(self.framebatch, text='Info', font='TkDefaultFont 10',
                                    height=1, width=15, command=self.info)
        self.buttoninfo.grid(row=3, column=0)
        # batch mode - button to calculate results
        self.buttoncalcbatch = self.tk.Button(self.frame, text='Apply IFS QSARs', font='TkDefaultFont 10',
                                         height=2, width=25, command=self.calculate_batch)
        self.buttoncalcbatch.grid(row=3, column=1)

    def toggle_disabled(self, setstate):
        """Disable or enable all widgets."""
        if hasattr(self, 'buttoninput'):
            self.buttoninput.config(state=setstate)
            self.buttoninput.update()
        if hasattr(self, 'textinput'):
            self.textinput.config(state=setstate)
            self.textinput.update()
        if hasattr(self, 'buttonoutput'):
            self.buttonoutput.config(state=setstate)
            self.buttonoutput.update()
        if hasattr(self, 'textoutput'):
            self.textoutput.config(state=setstate)
            self.textoutput.update()
        if hasattr(self, 'buttongotosingle'):
            self.buttongotosingle.config(state=setstate)
            self.buttongotosingle.update()
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.config(state=setstate)
            self.buttoninfo.update()
        if hasattr(self, 'buttoncalcbatch'):
            self.buttoncalcbatch.config(state=setstate)
            self.buttoncalcbatch.update()
        if hasattr(self, 'labelsmiles'):
            self.labelsmiles.config(state=setstate)
            self.labelsmiles.update()
        if hasattr(self, 'entrysmiles'):
            self.entrysmiles.config(state=setstate)
            self.entrysmiles.update()
        if hasattr(self, 'labelresult'):
            self.labelresult.config(state=setstate)
            self.labelresult.update()
        if hasattr(self, 'buttongotobatch'):
            self.buttongotobatch.config(state=setstate)
            self.buttongotobatch.update()
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.config(state=setstate)
            self.buttoninfo.update()
        if hasattr(self, 'buttoncalcsingle'):
            self.buttoncalcsingle.config(state=setstate)
            self.buttoncalcsingle.update()

    def info(self):
        """Display the module README.md in a popup window."""
        # create a popup to show information about the program
        popup = self.tk.Tk()
        popup.wm_title('Info')
        scrbar = self.tk.Scrollbar(popup)
        scrbar.pack(side=self.tk.RIGHT, fill=self.tk.Y)
        text = IFSGUIClass._ReadOnlyText(popup, font='Consolas 10')
        import os
        thispath = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(thispath, 'README.md'), 'r') as readmefile:
            infostring = readmefile.read()
        infostring = infostring.replace('<pre>', '').replace('</pre>', '').replace('\[', '[')
        text.insert('1.0', infostring)
        text.pack(side=self.tk.LEFT)
        scrbar.config(command=text.yview)
        text.config(yscrollcommand=scrbar.set)
        popup.mainloop()

    def calculate_single(self):
        """Get the SMILES from the GUI and apply QSARs."""
        # get smiles from the gui, apply models and write results to gui
        smiles = self.entrysmiles.get()
        self.toggle_disabled(self.tk.DISABLED)
        if self.mixturemode.get() == 'purechemical':
            localqsarlist = self.pure_qsarmodels
        elif self.mixturemode.get() == 'mixture':
            localqsarlist = self.mixture_qsarmodels
        results = apply_qsars_to_molecule(localqsarlist, smiles=smiles, converter=self.obcon, outformat='columns')
        # display results
        self.textresult.delete('1.0', self.tk.END)
        self.textresult.insert('1.0', results)
        self.toggle_disabled(self.tk.NORMAL)

    def calculate_batch(self):
        """Load SMILES from input file, apply QSARs, and write to output file."""
        # check for input and output files
        if self.inputfilename == '' or self.outputfilename == '':
            print('file not selected')
            return
        self.toggle_disabled(self.tk.DISABLED)
        # set input file type
        inextension = self.inputfilename.split('.')[-1]
        if inextension == 'txt':
            inseparator = '\t'
            inendline = '\n'
        elif inextension == 'csv':
            inseparator = ','
            inendline = '\n'
        else:
            print('file type not recognized: ', inextension)
            return
        # set output file type
        outextension = self.outputfilename.split('.')[-1]
        if outextension == 'txt':
            outseparator = '\t'
            outendline = '\n'
        elif outextension == 'csv':
            outseparator = ','
            outendline = '\n'
        else:
            print('file type not recognized: ', outextension)
            return
        if self.mixturemode.get() == 'purechemical':
            localqsarlist = self.pure_qsarmodels
        elif self.mixturemode.get() == 'mixture':
            localqsarlist = self.mixture_qsarmodels
        apply_qsars_to_molecule_list(localqsarlist,
                                     infilename=self.inputfilename,
                                     inheaderrows=1,
                                     inheadtrgtrow=1,
                                     inheadersmiles='smiles',
                                     inseparator=inseparator,
                                     inendline=inendline,
                                     converter=self.obcon,
                                     values=('insmi',
                                             'normsmi',
                                             'sminote',
                                             'OBMol',
                                             'endpoint',
                                             'units',
                                             'qsarpred',
                                             'UL',
                                             'error',
                                             'ULnote',
                                             'citation'),
                                     outfilename=self.outputfilename,
                                     outkeepdata=True,
                                     outformat='rows',
                                     outheader=True,
                                     outseparator=outseparator,
                                     outendline=outendline)
        self.toggle_disabled(self.tk.NORMAL)
        self.textinput.delete('1.0', self.tk.END)
        self.textoutput.delete('1.0', self.tk.END)

    def select_input_file(self):
        """Spawn a load file popup."""
        # create popup to select output file
        self.inputfilename = self.filedialog.askopenfilename(title='Select Input File',
                                                             filetypes=[('tab-delimited txt files', '*.txt'),
                                                                        ('csv files', '*.csv')])
        self.textinput.delete('1.0', self.tk.END)
        self.textinput.insert('1.0', self.inputfilename)

    def select_output_file(self):
        """Spawn a save file popup."""
        # create popup to select input file
        self.outputfilename = self.filedialog.asksaveasfilename(title='Select Output File',
                                                                defaultextension='.txt',
                                                                filetypes=[('tab-delimited txt files', '*.txt'),
                                                                           ('csv files', '*.csv')])
        self.textoutput.delete('1.0', self.tk.END)
        self.textoutput.insert('1.0', self.outputfilename)


def main():
    """Main loop for tk GUI"""
    app_manager = IFSGUIClass()

