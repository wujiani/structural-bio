############# TASK_1. Analysis of one single ensemble #############

import numpy as np
import math
import os
from Bio.PDB import DSSP, PPBuilder, Polypeptide, Superimposer
from sklearn import metrics
from sklearn.cluster import KMeans
# from task1 import feature_conf_dm

# a. Features of single conformation
# a.1 Radius of gyration
def radius_gyration(model):
    """
    Calculates the Radius of Gyration (Rg) of a protein in Angstroms.
    Does not use mass and assume heavy atoms have the same mass.

    https://en.wikipedia.org/wiki/Radius_of_gyration  (formula considering mass)                   # root mean square is existed
    https://link.springer.com/article/10.1134/S0026893308040195  (formula without mass)

    :param model: one conformation of the structure
    :return: single value of radius of gyration of one conformation
    """
    chain = model['A']
    # Heavy atoms coordinates
    coord = list()
    for atom in chain.get_atoms():
        if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
            coord.append(atom.get_coord())
    coord = np.array(coord)  # N X 3

    barycenter = np.sum(coord, axis=0) / coord.shape[0]  # center of mass is more correct

    # Calculate distance of each atom from the barycenter
    dist = coord - barycenter
    dist = dist * dist
    dist = np.sqrt(np.sum(dist, axis=1))  # rms(root mean square)
    # print(dist)

    return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)


# a.2 ASA
path_dssp = os.system("which mkdssp") # directory of dssp used to calculate asa in a.2
def calculateASA(model, model_path):
    """
    calculate ASA using dssp module

    :param model: one conformation of the structure
    :param model_path: pdb file of each model
    :return: a list containing ASA value for each residue of the conformation
    """
    dssp = DSSP(model, model_path, dssp=str(path_dssp))
    dssp_dict = dict(dssp)
    # ('A', (' ', index, ' ')) :
    # (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
    # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
    # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)

    return np.array([term[3] for term in dssp_dict.values()], dtype=float)


# a.3 SS
# parameters for ramachandran plot used to distinguish ss in a.3
rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                  (-180, 50, 80, 130, 'E', 'blue'),
                  (-100, -180, 100, 60, 'P', 'green'),
                  (-100, 50, 100, 130, 'P', 'green'),
                  (-180, -120, 180, 170, 'H', 'red'),
                  (0, -180, 180, 360, 'L', 'yellow')]

def secondary_structure(model, rama_ss_ranges):
    np.seterr(all='raise')
    """
    calculate ss according to phi, psi

    :param model: one conformation of the structure
    :param rama_ss_ranges: regions to distinguish ss
    :return: a list containing ss for each residue of the conformation
    """

    ppb = PPBuilder()  # PolyPeptideBuilder
    rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }
    ss = []
    # Calculate PSI and PHI
    for chain in model:
        pp = Polypeptide.Polypeptide(chain)
        phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
        for i, residue in enumerate(pp):
            rama.setdefault(chain.id, [[], [], []])
            rama[chain.id][0].append(residue)
            if phi_psi[i][0] is not None :
                rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
            else:
                rama[chain.id][1].append(phi_psi[i][0])
            if phi_psi[i][1] is not None:
                rama[chain.id][2].append(math.degrees(phi_psi[i][1]))
            else:
                rama[chain.id][2].append(phi_psi[i][1])

    for chain_id in rama:
        for residue, phi, psi in zip(*rama[chain_id]):
            ss_class = None
            if phi is not None and psi is not None:
                for x, y, width, height, ss_c, color in rama_ss_ranges:
                    if x <= phi < x + width and y <= psi < y + height:
                        ss_class = ss_c
                        break
            ss.append(ss_class)
    return ss, rama

## regions according to the paper. Convention from matplot lib for rectangles
## n.b not using for all the dssp convenction: P polyproline and L for left handed


# a.4 distance matrix
def get_distance_matrix(model, seq_sep=6):
    """
    Calculate the distance matrix for each residue
    :param model: one conformation of the structure
    :param seq_sep: def=6, calculate distance if sequence separation > seq_sep else dist=None
    :return: a matrix containing distances between pairs of residues
    """
    chain = model['A']
    distances = []
    for residue1 in chain:
        if residue1.id[0] == " ":  # Exclude hetero/water residues
            row = []
            for residue2 in chain:
                if residue2.id[0] == " ":  # Exclude hetero/water residues
                    if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                        row.append(residue1["CA"] - residue2["CA"])
                    else:
                        row.append(None)
            distances.append(row)
    return np.array(distances, dtype=float)

# RMSD of one ensemble
def structure_RMSD(structure, window_size = 5):
    # structure_rmsd= {}  # RMSD, no_models-1 X no_fragments X fragment_size
    structure_rmsd = []
    super_imposer = Superimposer()
    ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]  # coordinates of CA of this conformation
    # print(structure.get_models())
    for i, model in enumerate(structure):
        if i > 0:
            model_rmsd = []  # RMSD, no_fragment X fragment_size
            alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model

            # Iterate fragments
            for start in range(len(ref_model) - window_size + 1):
                end = start + window_size
                ref_fragment = ref_model[start:end]
                alt_fragment = alt_model[start:end]

                # Calculate rotation/translation matrices
                super_imposer.set_atoms(ref_fragment, alt_fragment)

                # Rotate-translate coordinates
                alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
                alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
                alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]
                # Calculate RMSD
                # https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
                ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
                dist = ref_fragment_coord - alt_fragment_coord
                # rmsd_fragment = np.sqrt(np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
                rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment
                model_rmsd.append(rmsd_res)
                # print(model_rmsd)
            model_rmsd = np.array(model_rmsd)  # no_fragments X fragment_size
            # Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
            model_rmsd = np.pad(model_rmsd, ((0, 0), (0, len(ref_model) - window_size)))
            # Roll the fragments one by one (add heading zeros)
            for i, row in enumerate(model_rmsd):
                model_rmsd[i] = np.roll(row, i)
            model_rmsd = np.array(model_rmsd)
            model_rmsd_average = np.average(model_rmsd, axis=0)
            structure_rmsd.append(model_rmsd_average)
            # print(structure_rmsd)
    rmsd_dict = {}
    for model_index in range(len(structure_rmsd)):
        rmsd_dict[str(model_index+1)] = structure_rmsd[model_index].tolist()

    return rmsd_dict

# clustering
def k_means_label(a,k_max):
    def km_index(k):
        pv = list(a)
        gf = np.array([pv]).T
        # from sklearn.cluster import KMeans
        y_pred = KMeans(n_clusters=k, random_state=9).fit_predict(gf)
        # print(y_pred)
        index = metrics.silhouette_score(gf, y_pred, metric='euclidean')
        return (index, (k,y_pred))
    cs = list(range(3, k_max+1))
    df = dict(map(km_index, cs))
    sorted_df = sorted(df.items(), key=lambda x: x[0], reverse=True)
    return sorted_df[0][1]

