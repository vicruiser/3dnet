# Load modules 
import numpy as np
import pandas as pd
import os
import sys
import math
from igraph import *
from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from joblib import Parallel, delayed, parallel_backend
from .parse_argv import parse_commandline
from .input_isfile import isfile


# parse PDB or mmCIF file
def parsePDB(pdb_file):
    parser = PDBParser(PERMISSIVE=1)
    id = os.path.splitext(os.path.basename(pdb_file))[0]
    structure = parser.get_structure(id,pdb_file)
    return(structure)

def parseMMCIF(cif_file): 
    parser = MMCIFParser()
    id = os.path.splitext(os.path.basename(cif_file))[0]
    structure  = parser.get_structure(id, cif_file)
    return(structure)

def downloadPDB(pdb_id, outdir): 
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir = os.path.join(outdir, 'pdbs'))
    structure = parseMMCIF(os.path.join(outdir,'pdbs', pdb_id + '.cif' ))
    return(structure)
    
# calculate centroids
def getCentroids(chain, plddt):
    if plddt is None: 
        plddt = -1
    centroids = dict()
    for residue in chain:
        tags = residue.get_full_id()
        if tags[3][0] == " ":
            arr = np.empty((0,3), int)
            for atom in residue:
                if atom.bfactor > float(plddt):
                    x,y,z = atom.get_coord()
                    arr = np.append(arr, np.array([[x, y, z]]), axis=0)
                    resnum = list(residue.id)[1]
                    centroids[resnum] = np.array(np.mean(arr, axis=0))
    return(centroids)

# build graph
def buildGraph(centroids,  sequence_separation, centroid_dist, outdir, pdb, chain):
    res_i_edges = []
    res_j_edges = []
    for res_i in centroids.keys():
        for res_j in centroids.keys():
            if (res_i < res_j) & ( (res_j - res_i ) > int(sequence_separation)): 
                dist = np.linalg.norm(centroids[res_i] - centroids[res_j])
                if float(dist) < float(centroid_dist): 
                    res_i_edges.append('X'+ str(res_i))
                    res_j_edges.append('X' + str(res_j))
    edges= pd.DataFrame({'ri': res_i_edges,'rj':res_j_edges})
    g = Graph.DataFrame(edges, directed=False)
    return(g)

def getLouvainPartition(g):
    louvain_partition = g.community_multilevel()
    return(louvain_partition)

def getSubgraph(g, louvain_partition, comm_size):
    clusters = [i for i in louvain_partition if len(i) > int(comm_size)]
    v = list(set(sum(clusters , [])))
    g2 = g.induced_subgraph(v)
    return(g2)

def changeBfactor(communities, chain):
    names = communities.name.unique()
    communities = communities.set_index('name') 
    for residue in chain: 
        if ('X' + str(residue.id[1])) in names:
            for atom in residue: 
                atom.set_bfactor(communities.at['X' + str(residue.id[1]),'membership'] )
        else: 
            for atom in residue: 
                atom.set_bfactor(0)
def overwrite(force, fname):
    if force is True:
        try:
            os.remove(fname)
        except:
            pass
    else: 
        if os.path.exists(fname):
            print('File ' + fname + '. To overwrite, select option -f/--force.')
            exit()

def saveCommPDB(pdb, fname):
    io = PDBIO()
    io.set_structure(pdb)
    io.save(fname)


def getPDBCommunities(pdb, outdir, centroid_dist, sequence_separation, pdbcomm, comm_size, plddt, force): 
    # iterate over all models in a pdb
    for model in pdb:
        model_id = str(model.id + 1 )
        # iterate over all chains
        for chain in model: 
            chain_length = len([_ for _ in chain.get_residues() if PDB.is_aa(_)])
            if chain_length > 0 :
                # calculate centroids 
                cs = getCentroids(chain, plddt)
                # build grapth based on set distance
                g = buildGraph(cs,  sequence_separation, centroid_dist, outdir, pdb, chain)
                # get partitions
                louvain_partition = getLouvainPartition(g)
                # remove nodes in community size < than selected threshold
                g2 = getSubgraph(g, louvain_partition, comm_size)
                # get new communities
                louvain_partition2 = getLouvainPartition(g2)    
                membership=louvain_partition2.membership
                # store as a data frame                
                communities = pd.DataFrame({'name': g2.vs['name'], 'membership':[x+1 for x in membership]})
                # save results 
                if plddt is None:    
                    fname = os.path.join(outdir, str(pdb.id) + '_model' + model_id +
                                                           '_chain' + chain.id +
                                                           '_centrdist' + str(centroid_dist) +
                                                           '_seqsep' + str(sequence_separation) +
                                                           '_mincommsize' + str(comm_size) +
                                                            '_communities.txt')
                else: 
                    fname = os.path.join(outdir, str(pdb.id) + '_model' + model_id +
                                                           '_chain' + chain.id +
                                                           '_centrdist' + str(centroid_dist) +
                                                           '_seqsep' + str(sequence_separation) +
                                                           '_mincommsize' + str(comm_size) +
                                                           '_plddt' + str(plddt) +
                                                            '_communities.txt')                                                        
                overwrite(force, fname)
                communities.to_csv(fname,
                                    sep = "\t",
                                    index= False, 
                                    header = False)
                # store pdb file with communities in B-factor?
                if pdbcomm is True: 
                    changeBfactor(communities, chain)
                    # names = communities.name.unique()
                    # communities = communities.set_index('name') 
                    # for residue in chain: 
                    #     if ('X' + str(residue.id[1])) in names:
                    #         for atom in residue: 
                    #             atom.set_bfactor(communities.at['X' + str(residue.id[1]),'membership'] )
                    #     else: 
                    #         for atom in residue: 
                    #             atom.set_bfactor(0)
    if pdbcomm is True:  
        if plddt is None:     
            fname = os.path.join(outdir, str(pdb.id) + '_model' + model_id +
                                                           '_chain' + chain.id +
                                                           '_centrdist' + str(centroid_dist) +
                                                           '_seqsep' + str(sequence_separation) +
                                                           '_mincommsize' + str(comm_size) +
                                                            '_communities.pdb')   
        else: 
            fname = os.path.join(outdir, str(pdb.id) + '_model' + model_id +
                                                           '_chain' + chain.id +
                                                           '_centrdist' + str(centroid_dist) +
                                                           '_seqsep' + str(sequence_separation) +
                                                           '_mincommsize' + str(comm_size) +
                                                           '_plddt' + str(plddt) +
                                                            '_communities.pdb')                                                        
        overwrite(force, fname)                                                            
        saveCommPDB(pdb, fname)                            
    # io = PDBIO()
    # io.set_structure(pdb)
    # io.save(os.path.join(outdir, str(pdb.id) +  '_communities.pdb'))

def execute_parallel(pdb_file, pdb_type, outdir, centroid_dist, sequence_separation, pdbcomm,comm_size, plddt, force):
    if pdb_type =="pdb":
        pdb = parsePDB(pdb_file)
    if pdb_type == "cif":
        pdb = parseMMCIF(file)
    getPDBCommunities(pdb, outdir, centroid_dist, sequence_separation, pdbcomm,comm_size, plddt, force)

def main():
    # parse command line options
    args = parse_commandline()

    if args.pdb is not None:
        input = args.pdb
    if args.cif is not None:
        input = args.cif
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    # handle input 
    print(isfile(input))
    if isfile(input) == 'yes':
        with open(input) as list_pdbs:
            print(list_pdbs)
            if args.parallel is False: 
                for file in list_pdbs:
                    if args.pdb is not None: 
                        pdb = parsePDB(file)
                    if args.cif is not None: 
                        pdb = parseMMCIF(file)
                    getPDBCommunities(pdb, args.out, args.centroid_dist, args.sequence_separation, args.pdbcomm, args.comm_size, args.plddt, args.force)
            else:
                if args.pdb is not None: 
                    pdb_type = "pdb"
                if args.cif is not None: 
                    pdb_type = "cif"
                Parallel(n_jobs=int(args.njobs))(delayed(execute_parallel)(file, 
                                                                 pdb_type,
                                                                   args.out,
                                                                   args.centroid_dist,
                                                                   args.sequence_separation,
                                                                   args.pdbcomm,
                                                                   args.comm_size, 
                                                                   args.plddt,
                                                                   args.force)
                                                       for file in list_pdbs)
    else: 
        if args.parallel is False: 
            if args.pdb is not None: 
                for file in input:
                    pdb = parsePDB(file)
                    getPDBCommunities(pdb, args.out, args.centroid_dist, args.sequence_separation, args.pdbcomm, args.comm_size, args.plddt, args.force)
            if args.cif is not None: 
                for file in input:
                    pdb = parseMMCIF(file)
                    getPDBCommunities(pdb, args.out, args.centroid_dist, args.sequence_separation, args.pdbcomm,args.comm_size,  args.plddt, args.force)
        else:
            if args.pdb is not None: 
                pdb_type = "pdb"
            if args.cif is not None: 
                pdb_type = "cif"
            Parallel(n_jobs=int(args.njobs))(delayed(execute_parallel)(file, 
                                                                 pdb_type,
                                                                   args.out,
                                                                   args.centroid_dist,
                                                                   args.sequence_separation,
                                                                   args.pdbcomm,
                                                                   args.comm_size, 
                                                                   args.plddt,
                                                                   args.force)
                                                       for file in input)