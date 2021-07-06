import json
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.cluster import hierarchy
from task2_defs import *

lines = []
with open(Path("data/output/output_c.json")) as f:
    for row in f.readlines():
        if row.strip().startswith("//"):
            continue
        lines.append(row)
content = json.loads("\n".join(lines))
ped_id = content["ped_id"]
ensemble_ids = content["ensemble_ids"]
path_ensembles = content["path_ensembles"]

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
with open(Path("data/output/output_c.json")) as f:
    for row in f.readlines():
        if row.strip().startswith("//"):
            continue
        lines.append(row)
content = json.loads("\n".join(lines))
######################################### TASK 2_output a ################################################
# paths which we use to write features of ensembles in json files
path_ens_rg = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"rg"))
path_ens_ss_entropy = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"ss_entropy"))
path_ens_med_asa = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"med_asa"))
path_ens_med_rmsd = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"med_rmsd"))
path_ens_med_dm = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"med_dm"))
path_ens_sd_dm = Path("data/output/{}/feature_ens_{}.json".format(ped_id,"sd_dm"))

# directory to contain features
dict_out_2_rg = {}
dict_out_2_ss_entropy = {}
dict_out_2_med_asa = {}
dict_out_2_med_rmsd = {}
dict_out_2_med_dm = {}
dict_out_2_sd_dm = {}

#calculate features of ensembles
#loop all ensembles
for ensemble_counter, [ens_id, path_ens] in enumerate(zip(ensemble_ids, path_ensembles)):


    rg = ens_radius_gyration(ens_id) # a.1
    ss_ent = ss_entropy(ens_id)  # a.2
    med_asa = median_asa(ens_id) # a.3
    med_rmsd = RMSD_median(ens_id)  # a.4
    med_dm = median_distance_matrix(ens_id)  # a.5
    sd_dm = sd_distance_matrix(ens_id)  # a.6

        # out_model = {"rg": rg, "asa": asa.tolist(), "ss": ss_list, "dm": distance_matrix.tolist()}
    dict_out_2_rg[ens_id] = rg
    dict_out_2_ss_entropy[ens_id] = ss_ent
    dict_out_2_med_asa[ens_id] = med_asa.tolist()
    dict_out_2_med_rmsd[ens_id] = med_rmsd.tolist()
    dict_out_2_med_dm[ens_id] = med_dm.tolist()
    dict_out_2_sd_dm[ens_id] = sd_dm.tolist()

# write features of ensembles in json files
with open(path_ens_rg, 'w') as f:
    json.dump(dict_out_2_rg, f, indent=2)
with open(path_ens_ss_entropy, 'w') as f:
    json.dump(dict_out_2_ss_entropy, f, indent=2)
with open(path_ens_med_asa, 'w') as f:
    json.dump(dict_out_2_med_asa, f, indent=2)
with open(path_ens_med_rmsd, 'w') as f:
    json.dump(dict_out_2_med_rmsd, f, indent=2)
with open(path_ens_med_dm, 'w') as f:
    json.dump(dict_out_2_med_dm, f, indent=2)
with open(path_ens_sd_dm, 'w') as f:
    json.dump(dict_out_2_sd_dm, f, indent=2)

######################################### TASK 2_output b ################################################

def ens_dRMS(ensemble_id_1, ensemble_id_2):
    median_dis_former = median_distance_matrix(ensemble_id_1)
    median_dis_latter = median_distance_matrix(ensemble_id_2)
    N = len(median_dis_former)
    # print(N)
    Diff_dis = np.zeros((N, N))

    n = 0  # no. of pairs
    for i in range(N):
        for j in range(i, N):
            if not np.isnan(median_dis_former[i][j]):
                n += 1
                Diff_dis[i][j] = median_dis_former[i][j] - median_dis_latter[i][j]
    return math.sqrt(np.sum(Diff_dis ** 2) / n)


heatmap_matrix = []
for ens_id_i in ensemble_ids:
    row = []
    for ens_id_j in ensemble_ids:
        row.append(ens_dRMS(ens_id_i , ens_id_j))
    heatmap_matrix.append(row)
# print(heatmap_matrix)
fig2, ax2 = plt.subplots(figsize=(12, 12))
im = ax2.imshow(heatmap_matrix)
fig2.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('data/output/{}/heatmap_distance_{}.png'.format(ped_id,ped_id), bbox_inches='tight')


# Set ticks
ax2.xaxis.set_major_locator(MultipleLocator(10))
ax2.xaxis.set_minor_locator(AutoMinorLocator(10))
ax2.yaxis.set_major_locator(MultipleLocator(10))
ax2.yaxis.set_minor_locator(AutoMinorLocator(10))


upper_tri_dis = []
for i in range(0,5):
    for j in range(i+1,5):
        upper_tri_dis.append(heatmap_matrix[i][j])
ytdist = np.array(upper_tri_dis)
Z = hierarchy.linkage(ytdist, 'single')
plt.figure()
dn = hierarchy.dendrogram(Z)
hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
plt.savefig('data/output/{}/dendrogram_distance_{}.png'.format(ped_id,ped_id))


######################################### TASK 2_output c ################################################
def read_feature_ens(feature_ens):
    lines = []
    with open(Path("data/output/{}/feature_ens_{}.json".format(ped_id,feature_ens))) as f:
        for row in f.readlines():
            if row.strip().startswith("//"):
                continue
            lines.append(row)
    return json.loads("\n".join(lines))

feature_ens_rg = read_feature_ens("rg")
feature_ens_ss_ent = read_feature_ens("ss_entropy")
feature_ens_med_asa = read_feature_ens("med_asa")
feature_ens_med_rmsd = read_feature_ens("med_rmsd")
feature_ens_med_dm = read_feature_ens("med_dm")
feature_ens_sd_dm = read_feature_ens("sd_dm")

# ensemble_ids = ["PED00142e025","PED00142e026","PED00142e027","PED00142e028","PED00142e029"]
N = len(feature_ens_ss_ent[ensemble_ids[0]])
def calcu_var(ens_feature):

    ens_var = []
    for res_index in range(1, N - 1):
        value_list = []
        for ensemble_counter, ens_id in enumerate(ensemble_ids):
            value_list.append(ens_feature[ens_id][res_index])
        var = np.var(np.array(value_list))
        # print(var)
        ens_var.append(var)
        # print(ens_var)
    # var =
    return np.array(ens_var)

ss_ent_var = calcu_var(feature_ens_ss_ent)
ss_ent_var = (ss_ent_var-np.mean(ss_ent_var))/np.std(ss_ent_var)
med_asa_var = calcu_var(feature_ens_med_asa)
med_asa_var = (med_asa_var-np.mean(med_asa_var))/np.std(med_asa_var)
med_rmsd_var = calcu_var(feature_ens_med_rmsd)
med_asa_var = (med_rmsd_var-np.mean(med_rmsd_var))/np.std(med_rmsd_var)

final_var = ss_ent_var + med_asa_var + med_rmsd_var
# print(final_var)
dict_var = dict(zip(range(1, N-1), [round(x,2) for x in final_var]))
sorted_var = sorted(dict_var.items(), key = lambda x:x[1])
my_var = sorted_var[:3]+sorted_var[-3:]
# print(my_var)

fig3, ax3 = plt.subplots(figsize=(12, 12))
ax3.scatter(range(1,N-1), final_var)
for x in my_var:
    plt.text(x[0], x[1] , x, fontsize=8)
# im = ax.imshow(heatmap_matrix)
# fig.colorbar(im, fraction=0.03, pad=0.05)
plt.savefig('data/output/{}/local_score_{}.png'.format(ped_id,ped_id), bbox_inches='tight')