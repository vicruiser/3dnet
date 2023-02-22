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
                               help='one or multiple PDB files provided \
                                    via command line or from a file')
    input_group.add_argument('-pdbcode', nargs='+', metavar="<String>", dest="pdbcode",
                               help='one or multiple PDB codes to fetch from PDB provided \
                                    via command line or from a file')                                    
    input_group.add_argument('-cif', nargs='+', metavar="<String>", dest="cif",
                               help='one or multiple CIF files provided \
                                    via command line or from a file')                                        
                                    
    parser.add_argument("--file", dest="file", action='store_true',
                        help="take input from FILE. Default is False", default=False)                           
    # create default output directory
    parser.add_argument("-o", "--outdir", metavar="<String>", dest="out",
                        help="output directory")
    parser.set_defaults(out=".")
  
    # force overwrite
    parser.add_argument('-f', "--force", dest="force", action='store_true',
                        help="force to owerwrite? Inactive by default", default=False)
    
    parser.add_argument('-plddt', dest="plddt", metavar="<float>", default = None,
                        help="pLDDT threshold if input is a AF2 model")

    parser.add_argument('-cd', "--centroid_dist", dest="centroid_dist", metavar="<float>", default = 7, 
                        help="threshold of maximum distance allowed in angstroms between two residue centroids (< 7A by default) ")

    parser.add_argument('-seqsep', "--sequence_separation", dest="sequence_separation", metavar="<int>", default = 1,
                        help="minimum linear sequence separation between two aminoacids (> 1 by default)")         
                       
    parser.add_argument('-cs', "--community_size", dest="comm_size", metavar="<int>", default = 1,
                        help="minimum number of elements to be considered part of the community ( >1 by default)")    
   
    parser.add_argument('-g', "--graph", dest="graph", action='store_true',
                        help="Save graph? Default is false", default=False)

    parser.add_argument('-gf', "--graph_format", dest="graph_format", metavar="<String>", default = 'edgelist',
                        help="If save graph opiton is active, select graph output format 'edgelist', 'pajek', 'ncol', 'gl',\
                         'graphml', 'dimacs', 'gml', 'dot', 'leda'. Default is 'edgelist' ")      
    
    parser.add_argument('-pdbcomm', dest="pdbcomm", action='store_true',
                        default=False,
                        help="save PDB file with B-factor is the number of the community? Default False")  

    parser.add_argument('-p', "--parallel", dest="parallel", action='store_true',
                        default=False,
                        help="parallelize process")

    parser.add_argument("-j", "--jobs", dest="njobs", metavar="<int>",
                        help="number of jobs to run in parallel")
    parser.set_defaults(njobs=1)


    # store arguments into variable
    args = parser.parse_args()

    # clean up (recommended)
    del(parser)

    return args
