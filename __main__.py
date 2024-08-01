"""
ifsqsar/__main__.py
developed by Trevor N. Brown
Implements a command line interface (CLI) that accesses the functionality of the ifsqsar python API
"""

if __name__ == "__main__":
    from . import ifsqsar
    from . import models
    import argparse
    import os

    thispath = os.path.dirname(os.path.abspath(__file__))
    # define options for commandline interface
    argparser = argparse.ArgumentParser(prog='ifsqsar',
                                        usage='%(prog)s [options]',
                                        description='IFSQSAR Command Line Interface, see full docs for usage examples',
                                        epilog='Invoke ifsqsar with no options start the GUI')
    # output full docs
    argparser.add_argument('-d',
                           '--docs',
                           action='store_true',
                           help='Prints the full documentation to standard output')
    # input must be either from file OR from a list of SMILES
    inputgroup = argparser.add_mutually_exclusive_group()
    inputgroup.add_argument('-s',
                            '--smiles',
                            metavar='smiles',
                            action='store',
                            type=str,
                            help='Comma-separated list of SMILES')
    inputgroup.add_argument('-i',
                            '--infile',
                            metavar='infile',
                            action='store',
                            type=str,
                            help='Path and name of input file')
    # input format specification; only used if input is from file
    argparser.add_argument('-r',
                           '--inheaderrows',
                           metavar='inheaderrows',
                           action='store',
                           type=int,
                           default=1,
                           help='Input file format: number of header rows, default = 1')
    argparser.add_argument('-t',
                           '--inheadtrgtrow',
                           metavar='inheadtrgtrow',
                           action='store',
                           type=int,
                           default=1,
                           help='Input file format: 1-indexed row number in which to look for SMILES label, '
                                'default = 1')
    argparser.add_argument('-c',
                           '--inheadersmiles',
                           metavar='inheadersmiles',
                           action='store',
                           type=str,
                           default='SMILES',
                           help='Input file format: label in header that indicates column containing SMILES, '
                                'default = SMILES')
    argparser.add_argument('-p',
                           '--inseparator',
                           metavar='inseparator',
                           action='store',
                           type=str,
                           default='\t',
                           help='Input file format: character used to separate columns, default = \\t (tab)')
    argparser.add_argument('-e',
                           '--inendline',
                           metavar='inendline',
                           action='store',
                           type=str,
                           default='\n',
                           help='Input file format: character used to separate rows, default = \\n (newline)')
    # output
    argparser.add_argument('-o',
                           '--outfile',
                           metavar='outfile',
                           action='store',
                           type=str,
                           help='Path and name of output file, if not specified print to standard output')
    # output format specification
    argparser.add_argument('-u',
                           '--outtossdata',
                           action='store_false',
                           default=True,
                           help='Output format: toss data from input file, '
                                'default = include data from input file in the output')
    argparser.add_argument('-a',
                           '--outspprsshead',
                           action='store_false',
                           default=True,
                           help='Output format: suppress header in output, '
                                'default = include header in the output')
    argparser.add_argument('-f',
                           '--outformat',
                           metavar='outformat',
                           action='store',
                           type=str,
                           default='rows',
                           help='Output format: chemical structures in rows, or in columns, '
                                'default = rows')
    argparser.add_argument('-m',
                           '--outseparator',
                           metavar='outseparator',
                           action='store',
                           type=str,
                           nargs='?',
                           default='\t',
                           const='',
                           help='Output format: character used to separate columns, default = \\t (tab)')
    argparser.add_argument('-n',
                           '--outendline',
                           metavar='outendline',
                           action='store',
                           type=str,
                           nargs='?',
                           default='\n',
                           const='',
                           help='Output format: character used to separate rows, default = \\n (newline)')
    # model selection
    argparser.add_argument('-q',
                           '--qsars',
                           metavar='qsars',
                           action='store',
                           type=str,
                           nargs='?',
                           default='default',
                           const='',
                           help='Comma-separated list of qsars to apply.'
                                'Single chemical QSARs:'
                                'fhlb, hhlb, hhlt, HLbiodeg, dsm, tm, tmpplfer, tmconsensus, '
                                'E, S, A, B, L, V, Vf, s, a, b, v, l, c, logKow, logKowdry, logKoa, '
                                'logKaw, logKoo, logVPliquid, logSwliquid, logSoliquid, '
                                'MVliquid, MVsolid, densityliquid, densitysolid, '
                                'MW, state. '
                                '"pure" for all QSARs that handle one SMILES input. '
                                'Mixture QSARs:'
                                'logKsa, '
                                '"mixture" for all QSARs that handle mixture SMILES input. '
                                'See full docs for explanation'
                           )
    # output value selection
    argparser.add_argument('-v',
                           '--values',
                           metavar='values',
                           action='store',
                           type=str,
                           nargs='?',
                           default='insmi,normsmi,sminote,endpoint,units,qsarpred,UL,error,ULnote,citation',
                           const='',
                           help='Comma-separated list of values to return. Full list: '
                                 'insmi, normsmi, sminote, eendpoint, units, qsarpred, UL, ULnote, error, citation. See full docs for explanation'
                           )
    # parse the options passed then decide actions
    args = argparser.parse_args()
    # no input specified, start the GUI
    if args.docs:
        with open(os.path.join(thispath, 'README.md'), 'r') as readmefile:
            infostring = readmefile.read()
        infostring = infostring.replace('<pre>', '').replace('</pre>', '').replace('\[', '[')
        print(infostring)
    elif args.infile is None and args.smiles is None:
        ifsqsar.main()
    else:
        # load QSARs
        if args.qsars == '':
            qsarmodels = []
        elif args.qsars == 'default':
            qsarmodels = models.get_qsar_list(qsarlist=['fhlb', 'hhlb', 'hhlt', 'HLbiodeg',
                                                        'dsm', 'tmconsensus', 'tbpplfer',
                                                        'logKow', 'logKowdry', 'logKoa', 'logKaw', 'logKoo',
                                                        'logVPliquid', 'logSwliquid', 'logSoliquid',
                                                        'MVliquid', 'densityliquid', 'MVsolid', 'densitysolid', 'MW',
                                                        'state',
                                                        'E', 'S', 'A', 'B', 'V', 'L',
                                                        's', 'a', 'b', 'v', 'l', 'c',
                                                        'logKsa'])
        elif args.qsars == 'pure':
            qsarmodels = models.get_qsar_list(qsarlist=['fhlb', 'hhlb', 'hhlt', 'HLbiodeg',
                                                        'dsm', 'tmconsensus', 'tbpplfer',
                                                        'logKow', 'logKowdry', 'logKoa', 'logKaw', 'logKoo',
                                                        'logVPliquid', 'logSwliquid', 'logSoliquid',
                                                        'MVliquid', 'densityliquid', 'MVsolid', 'densitysolid', 'MW',
                                                        'state',
                                                        'E', 'S', 'A', 'B', 'V', 'L',
                                                        's', 'a', 'b', 'v', 'l', 'c'])
        elif args.qsars == 'mixture':
            qsarmodels = models.get_qsar_list(qsarlist=['logKsa'])
        else:
            qsarmodels = models.get_qsar_list(qsarlist=args.qsars.split(','))
        # choose values
        values = []
        if args.values != '':
            values = args.values.split(',')
            while values[-1] == '':
                values.pop(-1)
        # SMILES list has been passed, apply models to list
        if args.smiles is not None:
            smiles = args.smiles.split(',')
            filename = None
        # infile name has been passed, apply models to file
        elif args.infile is not None:
            smiles = None
            filename = args.infile
        else:
            smiles = ''
            filename = ''
        # apply models
        if args.outseparator == '':
            outsep = '<nosep>'
        else:
            outsep = args.outseparator
        if args.outendline == '':
            outend = '<noend>'
        else:
            outend = args.outendline
        result = ifsqsar.apply_qsars_to_molecule_list(qsarmodels,
                                                      smiles,
                                                      infilename=filename,
                                                      inheadtrgtrow=args.inheadtrgtrow,
                                                      inheaderrows=args.inheaderrows,
                                                      inheadersmiles=args.inheadersmiles,
                                                      inseparator=args.inseparator,
                                                      inendline=args.inendline,
                                                      values=values,
                                                      outfilename=args.outfile,
                                                      outkeepdata=args.outtossdata,
                                                      outheader=args.outspprsshead,
                                                      outformat=args.outformat,
                                                      outseparator=outsep,
                                                      outendline=outend,
                                                      )
        # print to screen if no outfile
        if args.outfile is None:
            if args.outseparator == '':
                result = result.replace('<nosep>', '')
            if args.outendline == '':
                result = result.replace('<noend>', '')
            print(result, end='')


