from Bio.PDB import Superimposer
import numpy as np
from pathlib import  Path
import json
import pandas as pd
import math
import itertools
from task1_defs import *

#
#
# import numpy as np
# from Bio.PDB import PDBList, Superimposer
# from Bio.PDB.PDBParser import PDBParser
#
# pdb_id = '2k0e'
#
# # Fetch a PDB file to the current dir
# pdbl = PDBList()
# pdbl.retrieve_pdb_file(pdb_id, pdir='data/', file_format='pdb') # Will save to pdbXXXX.ent
#
# # Load the structure
# structure = PDBParser(QUIET=True).get_structure(pdb_id, "data/pdb{}.ent".format(pdb_id))

lines = []
with open(Path("data/output/output_c.json")) as f:
    for row in f.readlines():
        if row.strip().startswith("//"):
            continue
        lines.append(row)
content = json.loads("\n".join(lines))
ped_id = content["ped_id"]
# print(ped_id)
# structure_RMSD(structure)
def read_feature(feature_conf):
    lines = []
    with open(Path("data/output/{}/feature_conf_{}.json".format(ped_id,feature_conf))) as f:
        for row in f.readlines():
            if row.strip().startswith("//"):
                continue
            lines.append(row)
    return json.loads("\n".join(lines))
feature_conf_rg = read_feature("rg")
feature_conf_asa = read_feature("asa")
feature_conf_ss = read_feature("ss")
feature_conf_dm = read_feature("dm")

lines = []
with open(Path("data/RMSD_each_position/{}.json".format(ped_id))) as f:
    for row in f.readlines():
        if row.strip().startswith("//"):
            continue
        lines.append(row)
feature_ens_rmsd = json.loads("\n".join(lines))
# feature_ens_rmsd = json.\
# print(feature_ens_rmsd)
# a.1 radius of gyration for each conformation in one ensemble
def ens_radius_gyration(ensemble_id):
    """
    Gives an M (number of models for this ensemble) vector with radious of gyration values
    :param ensemble_id: id of one ensemble
    :return: M vector with radius of gyration values for each conformation in one ensemble
    """
    ens_rg = []  # list storing ensemble rg
    M = len(feature_conf_rg[ensemble_id])  # number of conformations in one ensemble
    for model_index in range(M):
        ens_rg.append(feature_conf_rg[ensemble_id][str(model_index)])
    return ens_rg

# a.2 SS entropy across ensemble conformations
def ss_entropy(ensemble_id):
    """
    Gives an N (number of conformations for this ensemble) vector with ss entropy
    :param ensemble_id: id of one ensemble
    :return: N vector with ss entropy for each position ensemble conformations
    """
    M = len(feature_conf_ss[ensemble_id])  # number of conformations in one ensemble
    N = len(feature_conf_ss[ensemble_id]['0'])  # number of residues in one conformation
    ss_ent = []  # list storing ss entropy

    for res_index in range(N):  # loop all residues
        ss_list = []  # list of ss for each residue position
        for model_index in range(M):  # loop all models
            ss = feature_conf_ss[ensemble_id][str(model_index)][res_index]
            if ss is not None:
                ss_list.append(ss)
        if len(ss_list) == 0:  # the first and last residue
            ss_ent.append(None)
        else:  # other residues
            ss_list_df = pd.DataFrame(ss_list)
            prob = ss_list_df[0].value_counts() / len(ss_list_df)  # compute frequency for each ss type
            # print(prob)
            shannon_entropy = 0
            for ss_type_index in range(len(prob)):
                shannon_entropy += (-prob[ss_type_index] * math.log(prob[ss_type_index], 2))
            ss_ent.append(shannon_entropy)

    return ss_ent

# a.3 Median solvent accessibility across ensemble conformations
def median_asa(ensemble_id):
    M = len(feature_conf_asa[ensemble_id])  # number of conformations in one ensemble
    N = len(feature_conf_asa[ensemble_id]['0'])  # number of residues in one conformation
    median_asa_list = []  # list storing median asa

    for residue_index in range(N):  # loop all residues
        asa_list = []  # list of asa for each residue position
        for model_index in range(M):  # loop all models
            asa_list.append(feature_conf_asa[ensemble_id][str(model_index)][residue_index])
        median_asa_list.append(np.median(asa_list))  # find the median of these distances across all models
    return np.array(median_asa_list, dtype=float)

# a.4 Median RMSD for each position across conformations
def RMSD_median(ensemble_id, window_size=5):
    """
    Calculate the median rmsd vector across the ensembles
    :param structure:
    :param window_size: def=5, is the window in which 5 residues are considered in order to form a structure fragment
    :return: median rmsd vector of the ensemble
    """
    M = len(feature_ens_rmsd[ensemble_id])  # number of conformations in one ensemble
    N = len(feature_ens_rmsd[ensemble_id]['1'])  # number of residues in one conformation

    median_rmsd_list = []  # list storing median asa
    for residue_index in range(N):  # loop all residues

        rmsd_list = []  # list of rmsd for each residue position
        for model_index in range(1, M ):  # loop all models
            rmsd_list.append(feature_ens_rmsd[ensemble_id][str(model_index)][residue_index])
        median_rmsd_list.append(np.median(rmsd_list))  # find the median of rmsd across all models

    return np.array(median_rmsd_list, dtype=float)

# a.5 Median distance of each pair of equivalent positions across conformations
def median_distance_matrix(ensemble_id):
    M = len(feature_conf_dm[ensemble_id])  # number of conformations in one ensemble
    N = len(feature_conf_dm[ensemble_id]['0'])  # number of residues in one conformation
    median_distances = [[[] for residue in range(N)] for residue in range(N)]  # list storing median distance

    # loop all comparisons between position pairs
    for comparison_index, (res1, res2) in enumerate(itertools.product(range(N), range(N))):
        distances_res1_res2 = []  # list of distances between position res1 and res2 across all models in one ensemble
        for model_index in range(M):  # loop all models
            distances_res1_res2.append(feature_conf_dm[ensemble_id][str(model_index)][res1][res2])
        median_distances[res1][res2] = np.median(
            distances_res1_res2)  # find the median of these distances across all models
    return np.array(median_distances, dtype=float)


 # a.6  distance of each pair of equivalent positions across conformations
def sd_distance_matrix(ensemble_id):
    M = len(feature_conf_dm[ensemble_id])  # number of conformations in one ensemble
    N = len(feature_conf_dm[ensemble_id]['0'])  # number of residues in one conformation
    sd_distances = [[[] for residue in range(N)] for residue in range(N)]  # list storing sd of distance

    # loop all comparisons between position pairs
    for comparison_index, (res1, res2) in enumerate(itertools.product(range(N), range(N))):
        distances_res1_res2 = []  # list of distances between position res1 and res2 across all models in one ensemble
        for model_index in range(M):  # loop all models
            distances_res1_res2.append(feature_conf_dm[ensemble_id][str(model_index)][res1][res2])
        sd_distances[res1][res2] = np.std(distances_res1_res2)  # find the sd of these distances across all models
    return np.array(sd_distances, dtype=float)