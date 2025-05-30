"""Clone of QSPR for PPLFER solute descriptor V - F adjusted"""
import numpy
value_names = ('Vf',)
version = 1
endpoint = 'Abraham PPLFER solute descriptor V - McGowan Volume - adjusted by Goss et al.'
citation = 'Mcgowan, J. C., '\
           'The Estimation of Solubility Parameters and Related Properties of Liquids. '\
           'J Chem Tech Biot A 1984, 34 (1), 38-42. '\
           'Adjustment for F: '\
           'Goss K-U, Bronner G, Harner T, Hertel M, Schmidt TC, ' \
           'The Partition Behavior of Fluorotelomer Alcohols and Olefins. ' \
           'Environ Sci Technol 40 (11):3572-3577.'
round_digits = 3
units = '0.01 cm^3/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
molecule_format = 'v1.0.0'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('[#6]', 'no description', 0.0),
                            ('[#7]', 'no description', 0.0),
                            ('[#8]', 'no description', 0.0),
                            ('[#9]', 'no description', 0.0),
                            ('[#14]', 'no description', 0.0),
                            ('[#15]', 'no description', 0.0),
                            ('[#16]', 'no description', 0.0),
                            ('[#17]', 'no description', 0.0),
                            ('[#5]', 'no description', 0.0),
                            ('[#32]', 'no description', 0.0),
                            ('[#33]', 'no description', 0.0),
                            ('[#34]', 'no description', 0.0),
                            ('[#35]', 'no description', 0.0),
                            ('[#50]', 'no description', 0.0),
                            ('[#51]', 'no description', 0.0),
                            ('[#52]', 'no description', 0.0),
                            ('[#53]', 'no description', 0.0),
                            ('[*H1]', 'no description', 0.0),
                            ('[*H2]', 'no description', 0.0),
                            ('[*H3]', 'no description', 0.0),
                            ('[*H4]', 'no description', 0.0),
                            ('[*H5]', 'no description', 0.0),
                            ('[*H6]', 'no description', 0.0),
                            ('[!#1]~[!#1]', 'no description', 0.0),
                            ], dtype=[('smarts', 'S11'), ('description', 'S14'), ('fragstdev', float)])
coefficientarrays = numpy.array([(0.1635,),
                                 (0.1439,),
                                 (0.1243,),
                                 (0.1248,),
                                 (0.2683,),
                                 (0.2487,),
                                 (0.2291,),
                                 (0.2095,),
                                 (0.1832,),
                                 (0.3102,),
                                 (0.2942,),
                                 (0.2781,),
                                 (0.2621,),
                                 (0.3935,),
                                 (0.3744,),
                                 (0.3614,),
                                 (0.3453,),
                                 (0.0215,),
                                 (0.043,),
                                 (0.0645,),
                                 (0.086,),
                                 (0.1075,),
                                 (0.129,),
                                 (-0.0656,),
                                 ], dtype=float)

stored = {'O': (0.1673, numpy.nan, numpy.nan, 'reliable value used', citation, units, endpoint)}

def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    return local_prediction, error

