import sys
import os
from Bio.PDB import PDBList
import networkx as nx
import biographs as bg
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pickle
import networkx.algorithms.clique as clique
import tools.sidechain_length as scl
import tools.amino_acids_conversion as aaconv
import glob
from biopandas.pdb import PandasPdb
from shutil import copyfile


num_atoms = {
    "ALA": 5,
    "ARG": 11,
    "ASN": 8,
    "ASP": 8,
    "CYS": 6,
    "GLU": 9,
    "GLN": 9,
    "GLY": 4,
    "HIS": 10,
    "ILE": 8,
    "LEU": 8,
    "LYS": 9,
    "MET": 8,
    "PHE": 11,
    "PRO": 7,
    "SER": 6,
    "THR": 7,
    "TRP": 14,
    "TYR": 12,
    "VAL": 7,
}
def number_atoms(amino_acid):
    return num_atoms[amino_acid]

def download_pdb(pdb_id, pdbs_path):
    """Downloads a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    """
    pdbl = PDBList(obsolete_pdb=True)
    pdbl.download_pdb_files(pdb_id, file_format="pdb", pdir=pdbs_path)


def get_pdb_path(pdb_id, pdbs_path="pdbs"):
    """Gets the location of a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    
    Returns
    -------
    str
        pdb file path
    str
        True if the pdb has been downloaded
    """

    if not pdbs_path:
        pdbs_path = "pdbs"

    downloaded = False

    package_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
    abs_path = os.path.join(package_path, "data", pdbs_path)
    abs_file_path = os.path.join(abs_path, pdb_id + ".*")

    if len(glob.glob(abs_file_path)) == 0:
        abs_file_path = os.path.join(abs_path, "pdb" + pdb_id + ".*")
        downloaded = True

        if len(glob.glob(abs_file_path)) == 0:
            os.makedirs(abs_path, exist_ok=True)
            download_pdb([pdb_id], abs_path)

        else:
            pdb_id = "pdb" + pdb_id
            abs_file_path = os.path.join(abs_path, pdb_id + ".*")
    pdb_path = glob.glob(abs_file_path)[0]

    return pdb_path, downloaded


def remove_hydrogens(pdb_file, save_copy=False):
    """Remove hydrogens from a pdb file
    
    Parameters
    ----------
    pdb_file : str
        pdb file path
    save_copy : boolean, optional
        if True, save a copy of the original file (default is False)
    """

    pdb_id = pdb_file.rsplit(".", 1)[0]
    copyfile(pdb_file, pdb_file.rsplit(".", 1)[0] + "_ORIG.pdb")

    with open(pdb_file, "r") as f:
        pdb = [line.rsplit("\n")[0] for line in f]
    pdb_new = []
    for line in pdb:
        if line[0:4] == "ATOM" and (line[13] == "H" or line[12] == "H"):
            pass
        else:
            pdb_new.append(line)
    with open(pdb_file, "w") as f:
        f.writelines([line + "\n" for line in pdb_new])


def list_relations(dim="all"):
    """Returns the list of relations to analyze
    
    Parameters
    ----------
    dim : str, optional
        "1D", "2D", "1-2D", "3-4D", "all" or "" (default is "all")
        "" is interpreted as "all"
    
    Returns
    -------
    list
        list of relations to analyze
    """

    if dim == "2D":
        rel_list = ["2D2", "2D3", "2D4"]
    elif dim == "1D":
        rel_list = ["1D"]
    elif dim == "not1D":
        rel_list = ["2D2", "2D3", "2D4", "3D", "4D"]
    elif dim == "1-2D":
        rel_list = ["1D", "2D2", "2D3", "2D4"]
    elif dim == "3-4D":
        rel_list = ["3D", "4D"]
    elif dim == "intramolecular":
        rel_list = ["1D", "2D2", "2D3", "2D4", "3D"]
    elif dim == "all" or dim == "":
        rel_list = ["1D", "2D2", "2D3", "2D4", "3D", "4D"]
    else:
        rel_list = [dim]

    return rel_list


def assign_secondary_structure(pdb):
    """Returns the secondary structure elements of a pdb
    
    Parameters
    ----------
    pdb : str
        pdb file path
    
    Returns
    -------
    dict
        dictionary of secondary structure elements
    """

    ppdb = PandasPdb().read_pdb(pdb)

    secondary_structure = {}

    helices_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "HELIX"][
        "entry"
    ]
    for helix in helices_from_pdb:
        identifier_h = helix[5:8].strip()
        initial_chain_h = helix[13].strip()
        initial_pos_h = helix[16:19].strip()
        final_pos_h = helix[28:31].strip()
        for i in range(int(initial_pos_h), int(final_pos_h) + 1):
            secondary_structure[initial_chain_h + str(i)] = (
                "helix" + identifier_h + "-" + initial_chain_h
            )

    sheets_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "SHEET"][
        "entry"
    ]
    for sheet in sheets_from_pdb:
        identifier_s = sheet[6:8].strip()
        initial_chain_s = sheet[15].strip()
        initial_pos_s = sheet[17:20].strip()
        final_pos_s = sheet[28:31].strip()
        for i in range(int(initial_pos_s), int(final_pos_s) + 1):
            secondary_structure[initial_chain_s + str(i)] = (
                "sheet" + identifier_s + "-" + initial_chain_s
            )

    mol = bg.Pmolecule(pdb)
    net = mol.network()

    residues_type = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname
        res_pos = residue.parent.id + str(residue.id[1])
        residues_type[res_pos] = res_type

    residues = list(net.nodes)  # assume they are ordered
    last_structure = None
    last_chain = None
    i = 0
    for residue in residues:
        chain = residue[0]
        try:
            structure = secondary_structure[residue]
            if structure != last_structure:
                i += 1
        except KeyError:
            if chain != last_chain:
                i += 1
            structure = "loop" + str(i)
            secondary_structure[residue] = structure
        last_structure = structure
        last_chain = chain

    return secondary_structure


def get_neighbor_structure_relation(secondary_structure, u, v):
    """Returns the relation (1D, 2D, 3D, 4D) between to neighboring nodes
    
    Parameters
    ----------
    secondary_structure : dict
        dictionary of secondary structure elements
    u: str
        node label
    v: str
        node label
    
    Returns
    -------
    str
        relation (1D, 2D, 3D, 4D)
    """

    chain_u = u[0]
    chain_v = v[0]
    pos_u = u[1::]
    pos_v = v[1::]
    struct_u = secondary_structure[u]
    struct_v = secondary_structure[v]

    if chain_u == chain_v:
        dist = int(pos_v) - int(pos_u)
        if np.abs(dist) == 1:
            relation = "1D"
        elif struct_u == struct_v:
            if np.abs(dist) < 5:
                relation = "2D" + str(np.abs(dist))
            else:
                relation = "3D"
        else:
            relation = "3D"
    else:
        relation = "4D"
        dist = 999

    return relation, dist


def list_aa(sort_by_size=False, sort_by_sidechain=False):
    """Imports the list of amino acids.
        Requires "aminoacids_size.txt" to be present in the directory.
    
    Parameters
    ----------
    sort_by_size : boolean, optional
        if True, sorts the amino acids by number of atoms (default is False)
    sort_by_sidechain : boolean, optional
        if True, sorts the amino acids by side chain length (default is False)
    
    Returns
    -------
    list
        list of amino acids
    """

    package_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]

    if sort_by_size:
        aa_path = os.path.join(package_path, "data", "aminoacids_size.txt")
    elif sort_by_sidechain:
        aa_path = os.path.join(package_path, "data", "aminoacids_sidechain.txt")
    else:
        aa_path = os.path.join(package_path, "data", "aminoacids.txt")

    with open(aa_path, "r") as f:
        amino_acids = [aa.rsplit("\n")[0] for aa in f]

    return amino_acids


def dict_aa_schl(amino_acids):
    """Creates a dictionary of amino acids side chain lengths.
    
    Parameters
    ----------
    amino_acids : list
        list of amino acids
    
    Returns
    -------
    dict
        dictionary of amino acids side chain lengths
    """

    amino_acids1l = [aaconv.three2one(aa) for aa in amino_acids]
    schl_dict = scl.dict_classif
    amino_acids_schl = [schl_dict[aa] for aa in amino_acids1l]

    return amino_acids_schl


def get_chains(mol):
    """Gets the chain of a Pmolecule object.
        
    Parameters
    ----------
    mol : Pmolecule
        Biographs Pmolecule object
    
    Returns
    -------
    list
        list of chains id's
    """

    chains = mol.model.child_list
    chains_id = [c.get_id() for c in chains]

    return chains_id


def create_aa_network(
    pdb_id,
    rel_list,
    folder_path,
    selected_positions=None,
    selected_chains=None,
    cutoff=5,
    pdbs_path="pdbs",
    remove_hydrogen_atoms=True,
):
    """Creates the amino acid network from a pdb id.
        
    Parameters
    ----------
    pdb_id : str
        pdb id of the protein
    rel_list: list
        list of relation (1D, 2D, 3D, 4D) to consider.
    folder_path: str
        path of the output folder
    selected_positions: None or list, optional
        list of sequence positions to consider. If None, all positions are considered (default is None)
    cutoff: int, optional
        cutoff threshold for the connection of nodes in the amino acid network (dafault is 5).
    remove_hydrogen_atoms: boolean, optional
        if True, saves removes the hydrogen atoms from the pdb file (default is True).
    Returns
    -------
    Graph
        NetworkX Graph object
    dict
        Node labels
    """

    if not pdbs_path:
        pdbs_path = "pdbs"

    amino_acids = list_aa()

    pdb, downloaded = get_pdb_path(pdb_id, pdbs_path)

    if remove_hydrogen_atoms:
        remove_hydrogens(pdb)

    mol = bg.Pmolecule(pdb)
    net = mol.network(cutoff=cutoff)

    edges_relation_dict = {}

    # take only selected positions:
    if selected_positions:
        for node in list(net.nodes):
            pos = int(node[1::])
            if pos not in selected_positions:
                net.remove_node(node)
    else:
        positions = [int(node[1::]) for node in list(net.nodes)]
        pos_min = min(positions)
        pos_max = max(positions)
        selected_positions = range(pos_min, pos_max + 1)

    if selected_chains:
        for node in list(net.nodes):
            chain = node[0]
            if chain not in selected_chains:
                net.remove_node(node)
    else:
        selected_chains = sorted(set([node[0] for node in list(net.nodes)]))

    secondary_structure = assign_secondary_structure(pdb)

    residues_dict = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname.strip()
        if len(res_type) < 3:
            res_type = aaconv.one2three(res_type)
        res_pos = residue.parent.id + str(residue.id[1])
        residues_dict[res_pos] = res_type

    for residue in mol.model.get_residues():
        adj_vector = [0] * 20
        weight_vector = [0] * 20
        node_name = residue.parent.id + str(residue.id[1])
        deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
        if deg == 0:
            net.remove_node(residue.parent.id + str(residue.id[1]))
        else:
            seqpos = residue.id[1]
            if seqpos not in selected_positions or residue.parent.id not in selected_chains:
                continue

            for neighbor in list(nx.neighbors(net, node_name)):
                relation, distance = get_neighbor_structure_relation(
                    secondary_structure, node_name, neighbor
                )
                # select only edges of desired relations
                if relation in rel_list:
                    edges_relation_dict[(node_name, neighbor)] = relation
                else:
                    net.remove_edge(neighbor, node_name)
    nx.set_edge_attributes(net, edges_relation_dict, 'relation')

    node_labels = {node: "%s%s:%s"%(aaconv.three2one(residues_dict[node]), node[1::], node[0]) for node in net.nodes}

    return net, node_labels

def plot_cutoffscreening(
    N,
    n_rows,
    n_cols,
    nodes_list,
    node_labels,
    cutoff_start,
    cutoff_stop,
    x,
    label,
    folder_path,
    file_name=None,
    ylabel=None,
    bars = False,
    barwidth=.2,
    markers=[""],
    linestyles=["-"],
    mask = True
):

    n_subfig = n_cols * n_rows

    cutoff_range = range(cutoff_start, cutoff_stop)

    if type(x) == list:
        x = np.array(x)

    if type(markers) == str:
        markers = [markers]
    if type(linestyles) == str:
        linestyles = [linestyles]
        
    

    if len(x.shape) > 2:
        single = False
        if not file_name or not ylabel:
            raise (
                Exception,
                "please provide a file_name and a ylabel when passing multiple datasets to plot.",
            )
        if len(markers) < x.shape[0]:
            m = markers[0]
            markers = [m] * x.shape[0]
        if len(linestyles) < x.shape[0]:
            l = linestyles[0]
            linestyles = [l] * x.shape[0]
    else:
        single = True
        if not file_name:
            file_name = label
        if not ylabel:
            ylabel = label

    schl_dict = scl.dict_classif
    schl_dict_num = scl.dict_size

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2, n_rows * 2))

    maxx = np.nanmax(x.flatten())
    minx = np.nanmin(x.flatten())

    i = 0
    j = 0
    for n in range(N):
        if j == n_cols:
            j = 0
            i += 1
        ax = axs[i][j]
        j += 1

        if single:
            if bars:
                ax.bar(cutoff_range, x[:, n], width=barwidth)
            else:
                ax.plot(cutoff_range, x[:, n], marker=markers[0], linestyle=linestyles[0])
        else:
            if bars:
                size_data = x.shape[0]
                if barwidth * size_data >= 1:
                    barwidth = 0.7 /  size_data
                shift = - (barwidth * size_data) / 2
                for xk, labelk in zip(x, label):
                    ax.bar(np.array(cutoff_range) + shift, xk[:, n], width=barwidth, label=labelk)
                    shift += barwidth
                minor_ticks = np.concatenate([np.array(cutoff_range) - 0.5, [cutoff_range[-1] + 0.5]])
                ax.set_xticks(minor_ticks, minor=True)
                ax.grid(True, axis='x', which='minor')
            else:
                for xk, labelk, markerk, linestylesk in zip(x, label, markers, linestyles):
                    ax.plot(cutoff_range, xk[:, n], label=labelk, marker=markerk, linestyle=linestylesk)

        node = nodes_list[n]
        node_label = node_labels[node]
        node_type = node_label[0]
        node_schl = schl_dict[node_type]
        if node_schl == "s":
            node_label_color = "#1f77b4"
        elif node_schl == "m":
            node_label_color = "#ff7f0e"
        else:
            node_label_color = "#2ca02c"
        ax.text(5, 0.8 * maxx, node_labels[node], color=node_label_color)
        ax.set_xticks(cutoff_range)
        ax.axes.set_xlim(cutoff_start - 0.5, cutoff_stop - 0.5)
        ax.axes.set_ylim(min(0, minx), 1.2 * maxx)
        ax.plot(cutoff_range, [0] * len(cutoff_range), "k", linewidth=.5)

        if mask:
            node_schl_num = schl_dict_num[node_type]
            diffmax = round(schl_dict_num["W"] - node_schl_num)
            ax.axvspan(5 + diffmax, cutoff_stop - 0.5, alpha=0.3, color="gray")

        if n == 0:
            ax.set_xlabel("cutoff")
            ax.set_ylabel(ylabel)
            if not single:
                ax.legend(loc=1, labelspacing=0, fontsize="small")

    for n in range(N + 1, n_subfig + 1):
        ax = axs[i][j]
        j += 1
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(folder_path, file_name + "_vs_cutoff_nodes.pdf"))
    plt.close()

def save_csv_cutoff(x, nodes_list, node_labels, cutoff_start, cutoff_stop, folder_path, file_name):
    with open(os.path.join(folder_path, file_name + ".csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['amino acid'] + ["cutoff %s" %c for c in range(cutoff_start, cutoff_stop)])
        for n, node in enumerate(nodes_list):
            label = node_labels[node]
            writer.writerow([label] + list(x[:, n]))

