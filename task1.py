import requests
import tarfile
import os
from pathlib import Path
import json
import random
import sys
sys.path.append("/lib/python3/dist-packages")
import pymol
from pymol import cmd
from Bio.PDB.PDBParser import PDBParser
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from task1_defs import *


# from definations import *
# from task2 import *

def downloadfiles(ped_id, num_ensembles=None):
    """
    Download the pdb files of the first num_ensembles of ped_id protein from the "proteinensemble.org" database

    :param ped_id: ProteinEnsembleDatabase id used to identify protein in the proteinensemble database
    :param num_ensembles: number of ensembles to download
    :return: ensembles_id downloaded, path to their pdb files
    """
    # global path_pdb, path_tar
    url = "https://proteinensemble.org/api/" + ped_id
    resp_json = requests.get(url).json()
    # print(resp_json["title"])
    ensembles_ids = []
    need2download = False  # check if all ensembles are already downloaded
    for curr_ensemble in resp_json["ensembles"][:num_ensembles]:
        ensembles_ids.append(curr_ensemble["ensemble_id"])
        if not need2download:
            n2d = True
            for files in os.listdir(path_pdb):
                if curr_ensemble["ensemble_id"] in files:
                    n2d = False
                    break
            if n2d:
                need2download = True

    ensembles_ids.sort()

    if need2download:
        # get direct link to the downloadable file (url as string)
        url = "https://proteinensemble.org/api/download"
        parameters = {
            "ensemble_id": ensembles_ids
        }
        download_link = requests.get(url, params=parameters).text
        resp_file = requests.get(download_link.replace('"', ''))

        # download ensembles
        with open(path_tar / '{}.tar'.format(ped_id), 'wb') as f:
            f.write(resp_file.content)

        fname = ped_id + '.tar'
        tar = tarfile.open(path_tar / fname)
        tar.extractall(path_tar)
        tar.close()

        dir = os.listdir(path_tar)
        for filename in dir:
            if filename == fname or ped_id not in filename:
                continue
            f = path_tar / filename
            tar2 = tarfile.open(f)
            tar2.extractall(path_pdb)
            tar2.close()

    return ensembles_ids, [path_pdb / f for f in sorted(os.listdir(path_pdb)) if ped_id in f]
    # # the main issue here is that  creates a subfolder with the .pdb file of our interest,
    # # exists an elegant way to extract directly ?


# # input of TASK_1
ped_id = input("Please input the ensemble id:") # ped_id input
num_ensembles = 5


###################################################################################################
###################################################################################################
# paths we need in two tasks
path_tar = Path("data/.tar_folder") / ped_id
path_pdb = Path("data/pdb_files") / ped_id
path_output = Path("data/output") / ped_id
path_RMSD = Path("data/RMSD_each_position")
path_pdb.mkdir(parents=True, exist_ok=True)
path_tar.mkdir(parents=True, exist_ok=True)
path_output.mkdir(parents=True, exist_ok=True)
path_RMSD.mkdir(parents=True, exist_ok=True)

# paths which we use to write features of conformations in json files
path_conf_rg = Path("data/output/{}/feature_conf_{}.json".format(ped_id,"rg"))
path_conf_asa = Path("data/output/{}/feature_conf_{}.json".format(ped_id,"asa"))
path_conf_ss = Path("data/output/{}/feature_conf_{}.json".format(ped_id,"ss"))
path_conf_dm = Path("data/output/{}/feature_conf_{}.json".format(ped_id,"dm"))
path_ens_RMSD = Path("data/RMSD_each_position/{}.json".format(ped_id))

###################################################################################################
###################################################################################################


###################################################################################################
###################################### output of TASK_1 ###########################################
###################################################################################################

ensemble_ids, path_ensembles = downloadfiles(ped_id, num_ensembles) # get list of ensemble ids and paths
# directory to contain features
dict_out_1_rg = {}
dict_out_1_asa = {}
dict_out_1_ss = {}
dict_out_1_dm = {}
dict_out_1_RMSD = {}



######################################### output a ################################################
# calculate features of conformations
# loop all ensembles
for ensemble_counter, [ens_id, path_ens] in enumerate(zip(ensemble_ids, path_ensembles)):
    with open(path_ens, 'r') as file:
        structure = PDBParser(QUIET=True).get_structure(ens_id, file) # load structure of each ensemble

    ########## make every conformations in one ensemble as one file ###########
    with open(path_ens, 'r') as f:
        my_pdb = f.read()

    # find all MODELs
    my_MODELs = []
    for each_letter in range(len(my_pdb) - len("MODEL")):
        letters = my_pdb[each_letter: each_letter + len("MODEL")]
        if letters == "MODEL":
            my_MODELs.append(each_letter)

    # find all ENDMDLs
    my_ENDMDLs = []
    for each_letter in range(len(my_pdb) - len("ENDMDL")):
        # print(each_letter)
        letters = my_pdb[each_letter: each_letter + len("ENDMDL")]
        if letters == "ENDMDL":
            my_ENDMDLs.append(each_letter + len("ENDMDL"))

    # make every conformations in one ensemble as one file
    # models = []
    for model_i in range(len(my_ENDMDLs)):
        each_model = my_pdb[my_MODELs[model_i]:my_ENDMDLs[model_i]]
        # models.append(each_model)
        with open(Path("data/pdb_files/{}/{}_model{}.pdb".format(ped_id, ens_id,  model_i + 1)), 'w') as f:
            f.write(each_model)


    # containing feature for all ensembles
    dict_conf_rg = {}
    dict_conf_asa = {}
    dict_conf_ss = {}
    dict_conf_dm = {}



    rmsd_per_ens = structure_RMSD(structure)
    # print(rmsd_per_ens)
    # loop all conformations within one ensemble
    for model_counter, model in enumerate(structure):
        # out_model = {} # containing feature for each ensemble
        rg = radius_gyration(model) # a.1
        asa = calculateASA(model, Path("data/pdb_files/{}/{}_model{}.pdb".format(ped_id, ens_id, model_counter + 1)))
        ss_list = secondary_structure(model,rama_ss_ranges)[0]
        distance_matrix = get_distance_matrix(model)

        dict_conf_rg[model_counter] = rg
        dict_conf_asa[model_counter] = asa.tolist()
        dict_conf_ss[model_counter] = ss_list
        dict_conf_dm[model_counter] = distance_matrix.tolist()

    dict_out_1_rg[ens_id] = dict_conf_rg
    dict_out_1_asa[ens_id] = dict_conf_asa
    dict_out_1_ss[ens_id] = dict_conf_ss
    dict_out_1_dm[ens_id] = dict_conf_dm
    dict_out_1_RMSD[ens_id] = rmsd_per_ens

# write features of conformations in json files
with open(path_conf_rg, 'w') as f:
    json.dump(dict_out_1_rg, f, indent=2)
with open(path_conf_asa, 'w') as f:
    json.dump(dict_out_1_asa, f, indent=2)
with open(path_conf_ss, 'w') as f:
    json.dump(dict_out_1_ss, f, indent=2)
with open(path_conf_dm, 'w') as f:
    json.dump(dict_out_1_dm, f, indent=2)
with open(path_ens_RMSD, 'w') as f:
    json.dump(dict_out_1_RMSD, f, indent=2)

######################################### output b ################################################
# read features from files generated in TASK_1
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
all_models = []

cluster_asa = {} # mean asa value for clustering models
for ens_id in ensemble_ids:
    M = len(feature_conf_asa[ens_id])
    for model_index in range(M):
        len_0 = 4 - len(str(model_index))
        model_number = '0' * len_0 + str(model_index)
        my_model = "{}_{}".format(ens_id, model_number)
        all_models.append(my_model)
        cluster_asa[my_model] = np.average(feature_conf_asa[ens_id][str(model_index)])
data = cluster_asa.values()
my_cluster = k_means_label(data,8)[1] # 8 is the largest k to cluster

# get the clustering outcome with k which is calculated automatically before
cluster = {}
for i in range(len(all_models)):
    if str(my_cluster[i]) not in cluster.keys():
        cluster[str(my_cluster[i])] = []
    cluster[str(my_cluster[i])].append(all_models[i])
sorted_cluster = dict(sorted(cluster.items(), key = lambda x:int(x[0]))) # the clustering outcome

# get the representative conformations
my_confs = []
for cluster_index in range(len(sorted_cluster)):
    my_confs.append(random.choice(sorted_cluster[str(cluster_index)]))

rel_dm = np.average([feature_conf_dm[model[:-5]][str(int(model[-4:]))] for model in my_confs], axis = 0)
# calculate distance between 2 models
def RMSD(model_1, model_2):
    dm_former = feature_conf_dm[model_1[:-5]][str(int(model_1[-4:]))]
    if model_2 == "rel":
        dm_latter = rel_dm
    else:
        dm_latter = feature_conf_dm[model_2[:-5]][str(int(model_2[-4:]))]
    N = len(dm_former)
    # print(N)
    Diff_dis = np.zeros((N, N))

    n = 0  # no. of pairs
    for i in range(N):
        for j in range(i, N):
            if not np.isnan(dm_former[i][j]):
                n += 1
                Diff_dis[i][j] = dm_former[i][j] - dm_latter[i][j]
    return math.sqrt(np.sum(Diff_dis ** 2) / n)
my_confs_rmsd = [RMSD(model, "rel") for model in my_confs]
# my_confs_rmsd_standard = (np.array(my_confs_rmsd) - np.mean(my_confs_rmsd))/ np.std(my_confs_rmsd)
my_confs_rg = [feature_conf_rg[model[:-5]][str(int(model[-4:]))] for model in my_confs]
my_confs_asa = [np.average(feature_conf_asa[model[:-5]][str(int(model[-4:]))]) for model in my_confs]
my_confs_coord = list(zip(my_confs_rmsd, my_confs_rg, my_confs_asa))

fig1=plt.figure()
ax1 = Axes3D(fig1)
ax1.scatter3D(my_confs_rmsd,my_confs_rg,my_confs_asa, cmap='Blues')
for i in range(len(my_confs_rmsd)):
    for j in range(i, len(my_confs_rmsd)):
        ax1.plot3D((my_confs_rmsd[i], my_confs_rmsd[j]),(my_confs_rg[i], my_confs_rg[j]), (my_confs_asa[i], my_confs_asa[j]),'gray')
plt.savefig('data/output/{}/task1_output_b_{}_1.png'.format(ped_id,ped_id), bbox_inches='tight')
ax1.view_init(30, 30)
plt.savefig('data/output/{}/task1_output_b_{}_2.png'.format(ped_id,ped_id), bbox_inches='tight')
ax1.view_init(-90, -90)
plt.savefig('data/output/{}/task1_output_b_{}_3.png'.format(ped_id,ped_id), bbox_inches='tight')


######################################### output c ################################################

# confs_rg = [feature_conf_rg[model[:-5]][str(int(model[-4:]))] for model in my_confs]
confs_asa = [feature_conf_asa[model[:-5]][str(int(model[-4:]))] for model in my_confs]
# confs_rg_standard = []
confs_asa_standard = []
for i in range(len(confs_asa)):
    # confs_rg_standard.append(np.array(confs_rg - np.mean(confs_rg)) / np.std(confs_rg))
    confs_asa_standard.append(np.array(confs_asa[i] - np.mean(confs_asa[i])) / np.std(confs_asa))
# print(confs_asa_standard)
res_low_var = np.std(confs_asa_standard, axis = 0)
var_dict = dict(zip(range(len(confs_asa[0])), res_low_var))
del var_dict[0]
del var_dict[len(var_dict.keys())]
low_vari_resi = sorted(var_dict.items(), key = lambda x:x[1])[0][0]
order_models = []
for asa in confs_asa:
    order_dict = dict(zip(range(len(confs_asa[0])), np.array(asa) - np.std(asa)))
    del order_dict[0]
    del order_dict[len(order_dict.keys())]
    order = dict(sorted(order_dict.items(), key=lambda x: x[1]), reverse = True)
    del order["reverse"]
    order_models.append(list(order.keys()))


path_output_c = Path("data/output/output_c.json")
content = {}
content["my_confs"] = my_confs
content["low_vari_resi"] = low_vari_resi
content["order_models"] = order_models
content["ensemble_ids"] = ensemble_ids
content["ped_id"] = ped_id
content["path_ensembles"] = [str(a) for a in path_ensembles]
# print(content)
with open(path_output_c, 'w') as f:
    json.dump(content, f, indent=2)
