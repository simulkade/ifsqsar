"""Partial clone of QSRP for MV, intended for internal use"""
import numpy
value_names = ('MVmlrRings',)
version = 1
endpoint = 'Molar volume - fraction of model for internal use'
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'cm^3/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
molecule_format = 'v1.0.0'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([
                            ('sssr', 'table ', 0.0),
                            ('C1C3CC2CC(CC1C2)C3', 'table 5', 0.0),
                            ('C12C3C4C1C5C2C3C45', 'table 5', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([
                                 (-3.89,),
                                 (3.50+(3.5-6.98),),
                                 (38.24+(38.24-33.94),)
                                 ], dtype=float)

stored = {}

def post_processing(prediction, error):
    if prediction != 0:
        prediction += 15.15
    prediction = max(-15.97, prediction)
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

