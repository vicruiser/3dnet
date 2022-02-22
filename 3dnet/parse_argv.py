# -*- coding: utf-8 -*-
import argparse
import os


def parse_commandline():
    '''
    Parse inputs from command line.

    Returns
    -------
    args
        arguments to give to the functions
    '''
    description = '''
    -----------------------------
     ____  _____             _   
    |___ \|  __ \           | |  
      __) | |  | |_ __   ___| |_ 
    |__ <| |  | | '_ \ / _ \ __|
    ___) | |__| | | | |  __/ |_ 
    |____/|_____/|_| |_|\___|\__|    
    -----------------------------
    '''


    # innit parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description)

    # pdb, cif or model
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-pdb', nargs='+', metavar="<String>", dest="pdb",
                               help='One or multiple PDB files provided \
                                    via command line or from a file')
    input_group.add_argument('-cif', nargs='+', metavar="<String>", dest="cif",
                               help='One or multiple CIF files provided \
                                    via command line or from a file')                                
    # create default output directory
    parser.add_argument("-o", "--outdir", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out=".")
    
    # force overwrite
    parser.add_argument('-f', "--force", dest="force", action='store_true',
                        help="force to owerwrite? Inactive by default", default=False)
    
    # filter by distance (applicable to interfaces)
    parser.add_argument('-plddt', dest="plddt", metavar="<float>", default = None,
                        help="pLDDT threshold if input is a AF2 model")

    # filter by distance (applicable to interfaces)
    parser.add_argument('-cd', "--centroid_dist", dest="centroid_dist", metavar="<float>", default = 7, 
                        help="Threshold of maximum distance allowed in angstroms between two centroids ")

    # filter by distance (applicable to interfaces)
    parser.add_argument('-seqsep', "--sequence_separation", dest="sequence_separation", metavar="<float>", default = 1,
                        help="Minimum linear sequence separation between two aminoacids")         
                       

    # filter by distance (applicable to interfaces)
    parser.add_argument('-cs', "--community_size", dest="comm_size", metavar="<float>", default = 1,
                        help="Minimum number of elements to be considered part of the community")         
    
    # create default output directory
    parser.add_argument('-pdbcomm', dest="pdbcomm", action='store_true',
                        default=False,
                        help="Save PDB file with B-factor is the number of the community? Default False")  
    # create default output directory
    parser.add_argument('-p', "--parallel", dest="parallel", action='store_true',
                        default=False,
                        help="Parallelize process")

    # interfaces database file
    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel")
    parser.set_defaults(njobs=1)


    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
