import sys
sys.path.append("/lib/python3/dist-packages")
import pymol
from pymol import cmd
from pathlib import Path
import numpy as np
import json
# from colour import Color
import random

######################################### TASK 1_output c ################################################
lines = []
with open(Path("data/output/output_c.json")) as f:
    for row in f.readlines():
        if row.strip().startswith("//"):
            continue
        lines.append(row)
content = json.loads("\n".join(lines))

my_confs = content["my_confs"]
low_vari_resi = content["low_vari_resi"]
order_models = content["order_models"]
ensemble_ids = content["ensemble_ids"]
ped_id = content["ped_id"]

color_ribbons = ["gray", "olive", "purple", "green", "marine", "orange", "gray", "red", "cyan", "pink"]
# colors = list(Color("blue").range_to(Color("red"), 120))
# for color_index in range(120):
#     colors[color_index] = colors[color_index].hex_l
sizes = list(np.linspace(0.1, 1.2, len(order_models[0])))

# order_model = []
# for i in range(len(model_input)):
#     orders = list(range(120))
#     random.shuffle(orders)
#     order_model.append(orders)
# from new_output import ensemble_ids


pymol.finish_launching()  # Open Pymol
# path_ped = "data/pdb_files/PED00142/PED00142e026.pdb"
for ens_id in ensemble_ids:
# ens_id = "PED00142e025"
    path_ped = "data/pdb_files/PED00142/{}.pdb".format(ens_id)
    cmd.load(path_ped, multiplex=1)
cmd.hide(representation="everything")

def pymol_show(model_input, low_vari_resi, order):
    sele_name_model = "sele_models"
    cmd.select(sele_name_model, ",".join(model_input))
    cmd.show(representation="ribbon", selection=sele_name_model)

    # cmd.set("ribbon_color", "red", sele_name_model)
    for i in range(len(model_input)):
        cmd.set("ribbon_color", color_ribbons[i % 10], model_input[i])
        if i > 0:
            cmd.align("{} & resi {}".format(model_input[i], low_vari_resi), "{} & resi {}".format(model_input[0], low_vari_resi))
        # cmd.align("PED00142e026_0049 & resi 40", "PED00142e026_0029 & resi 40")
    #
    #      # show sphere, resi 40
        for resi_index in range(1, len(order[i])+1):

    #
            sele_name = "{}_resi_{}".format(model_input[i], resi_index)
            res = '{} & resi {} & name CA'.format(model_input[i], resi_index)
    # #         # res2 = 'PED00142e026_0045 & resi {} & name CA'.format(resi_index)
    # #         # res3 = 'PED00142e026_0049 & resi {} & name CA'.format(resi_index)
            cmd.select(sele_name, res)
            cmd.show("spheres", sele_name)
            cmd.color(color_ribbons[i % 10], sele_name)
            profile = order[i][resi_index-1]-1

    #
    #         # cmd.color("0x" + str(colors[profile])[1:], sele_name)
    #
            cmd.set("sphere_scale", sizes[profile], sele_name)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.zoom()
    cmd.orient(sele_name_model)
    cmd.png("data/output/{}/pymol_image.png".format(ped_id), width=5000, dpi=1000, ray=1)
    cmd.quit()


pymol_show(my_confs,low_vari_resi, order_models)

