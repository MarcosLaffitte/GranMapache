################################################################################
#                                                                              #
#  - Analysis of Hydrogen Extensions                                           #
#                                                                              #
#  - Isomorphism of ITS graphs                                                 #
#                                                                              #
#  - for AMB - BMC - extension of WABI 2025                                    #
#                                                                              #
#  - Made by Marcos Laffitte - Github @MarcosLaffitte                          #
#                                                                              #
#  - 3 - Relabeling strategy with NetworkX                                     #
#                                                                              #
################################################################################




# dependencies #################################################################




# already in python ------------------------------------------------------------
import sys
import math
import time
import random
import pickle
from copy import deepcopy
from itertools import permutations
from operator import eq, itemgetter




# additional dependencies ------------------------------------------------------
import cython
import scipy
import numpy as np
import pandas as pd
import networkx as nx
import pysmiles as ps
import cmasher as cmr
import seaborn as sns
from rdkit import Chem
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D




# networkx isomorphism with node and edge match --------------------------------
from operator import eq
from networkx.algorithms import isomorphism
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match




# long's synkit ----------------------------------------------------------------
from synkit.IO import rsmi_to_its, rsmi_to_graph, load_from_pickle
from synkit.Graph.Hyrogen._misc import check_hcount_change
from synkit.Graph.ITS.its_decompose import get_rc, its_decompose




# custom dependencies ----------------------------------------------------------
import gmapache as gm




# warnings ---------------------------------------------------------------------
import logging
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)




# parameters ###################################################################




# control ----------------------------------------------------------------------
print_all = False




# input ------------------------------------------------------------------------




# functions ####################################################################




# function: determine if a given atom map fully maps all atoms -----------------
def check_eveness(mapped_reaction):

    # local variables
    atoms_R = 0
    atoms_P = 0
    products = []
    reactants = []
    mapping_R = []
    mapping_P = []
    aam_attribute = dict()

    # get mapping on reactants side
    reactants = list((mapped_reaction.split(">>")[0]).split("."))
    for each_reactant in reactants:
        mol = ps.read_smiles(each_reactant, strict = False)
        aam_attribute = nx.get_node_attributes(mol, "class")
        for atom in list(mol.nodes()):
            atoms_R = atoms_R + 1
            if(atom in list(aam_attribute.keys())):
                mapping_R.append(aam_attribute[atom])
    mapping_R.sort()

    # get mapping on products side
    products = list((mapped_reaction.split(">>")[1]).split("."))
    for each_product in products:
        mol = ps.read_smiles(each_product, strict = False)
        aam_attribute = nx.get_node_attributes(mol, "class")
        for atom in list(mol.nodes()):
            atoms_P = atoms_P + 1
            if(atom in list(aam_attribute.keys())):
                mapping_P.append(aam_attribute[atom])
    mapping_P.sort()

    # determine consistency
    if(mapping_R == mapping_P):
        if((len(mapping_R) == atoms_R) and (len(mapping_R) == atoms_P)):
            return(True)

    # end of function
    return(False)




# function: get base-graph and hydro-graph representations of SMILES -----------
def get_graph_representations(mapped_reaction):

    # local variables
    hydro_balanced = True
    hydro_index = 0
    hydro_count_R = 0
    hydro_count_P = 0
    temp_mol_str = ""
    products_side = ""
    reactants_side = ""
    mapped_atoms = []
    products_list = []
    reactants_list = []
    base_hcount = dict()
    atom_mapping = dict()
    hydro_mapping = dict()
    G = None
    H = None
    temp_mol = None
    full_hydro_G = None
    full_hydro_H = None
    each_product = None
    each_reactant = None
    molecule_hydro = None

    # split reaction
    reactants_side = (mapped_reaction.split(">>"))[0]
    products_side = (mapped_reaction.split(">>"))[1]

    # get reactants graph G
    reactants_list = reactants_side.split(".")
    G = nx.Graph()
    full_hydro_G = nx.Graph()

    # turn each reactant into a graph
    for each_reactant in reactants_list:

        # read smiles into networkx graph
        try:
            temp_mol = Chem.MolFromSmiles(each_reactant, sanitize = True)
            temp_mol_str = Chem.MolToSmiles(temp_mol,
                                            canonical = True,
                                            kekuleSmiles = True,
                                            allHsExplicit = True)
        except:
            temp_mol = Chem.MolFromSmiles(each_reactant, sanitize = False)
            temp_mol_str = Chem.MolToSmiles(temp_mol,
                                            canonical = True,
                                            kekuleSmiles = True,
                                            allHsExplicit = True)

        # get molecule with implicit hydrogens
        molecule = ps.read_smiles(temp_mol_str,
                                  explicit_hydrogen = False,
                                  reinterpret_aromatic = True,
                                  strict = False)

        # get molecule with explicit hydrogens
        molecule_hydro = ps.read_smiles(temp_mol_str,
                                        explicit_hydrogen = True,
                                        reinterpret_aromatic = True,
                                        strict = False)

        # add hcount to hydro representation
        base_hcount = nx.get_node_attributes(molecule, "hcount")
        nx.set_node_attributes(molecule_hydro, base_hcount, name = "hcount")

        # get existing atom map
        atom_mapping = nx.get_node_attributes(molecule, "class")
        mapped_atoms = list(atom_mapping.keys())

        # rename hydrogens
        hydro_mapping = dict()
        for v in list(molecule_hydro.nodes()):
            if((molecule_hydro.nodes[v]["element"] == "H") and (not v in mapped_atoms)):

                # count hydrogens for testing
                hydro_count_R = hydro_count_R + 1

                # set new name for hydrogen
                hydro_index = hydro_index + 1
                hydro_mapping[v] = "h." + str(hydro_index)

        # assign new names for hydrogens
        molecule_hydro = nx.relabel_nodes(molecule_hydro, hydro_mapping)

        # rename originally matched vertices using atom map
        molecule = nx.relabel_nodes(molecule, atom_mapping)
        molecule_hydro = nx.relabel_nodes(molecule_hydro, atom_mapping)

        # make (disjoint) union over G
        G = deepcopy(nx.union(G, molecule))
        full_hydro_G = deepcopy(nx.union(full_hydro_G, molecule_hydro))

    # get products graph H
    products_list = products_side.split(".")
    H = nx.Graph()
    full_hydro_H = nx.Graph()

    # turn each product into a graph
    for each_product in products_list:

        # read smiles into networkx graph
        try:
            temp_mol = Chem.MolFromSmiles(each_product, sanitize = True)
            temp_mol_str = Chem.MolToSmiles(temp_mol,
                                            canonical = True,
                                            kekuleSmiles = True,
                                            allHsExplicit = True)
        except:
            temp_mol = Chem.MolFromSmiles(each_product, sanitize = False)
            temp_mol_str = Chem.MolToSmiles(temp_mol,
                                            canonical = True,
                                            kekuleSmiles = True,
                                            allHsExplicit = True)

        # get molecule with implicit hydrogens
        molecule = ps.read_smiles(temp_mol_str,
                                  explicit_hydrogen = False,
                                  reinterpret_aromatic = True,
                                  strict = False)

        # get molecule with explicit hydrogens
        molecule_hydro = ps.read_smiles(temp_mol_str,
                                        explicit_hydrogen = True,
                                        reinterpret_aromatic = True,
                                        strict = False)

        # add hcount to hydro representation
        base_hcount = nx.get_node_attributes(molecule, "hcount")
        nx.set_node_attributes(molecule_hydro, base_hcount, name = "hcount")

        # get existing atom map
        atom_mapping = nx.get_node_attributes(molecule, "class")
        mapped_atoms = list(atom_mapping.keys())

        # rename hydrogens
        hydro_mapping = dict()
        for v in list(molecule_hydro.nodes()):
            if((molecule_hydro.nodes[v]["element"] == "H") and (not v in mapped_atoms)):

                # count hydrogens for testing
                hydro_count_P = hydro_count_P + 1

                # set new name for hydrogen
                hydro_index = hydro_index + 1
                hydro_mapping[v] = "h." + str(hydro_index)

        # assign new names for hydrogens
        molecule_hydro = nx.relabel_nodes(molecule_hydro, hydro_mapping)

        # rename originally matched vertices using atom map
        molecule = nx.relabel_nodes(molecule, atom_mapping)
        molecule_hydro = nx.relabel_nodes(molecule_hydro, atom_mapping)

        # make (disjoint) union over H
        H = deepcopy(nx.union(H, molecule))
        full_hydro_H = deepcopy(nx.union(full_hydro_H, molecule_hydro))

    # remove reaction (later) if not balanced up to hydrogens
    if(not hydro_count_R == hydro_count_P):
        hydro_balanced = False

    # remove unnecessary attributes of G
    for (v, info) in G.nodes(data = True):
        if("_atom_str" in G.nodes[v]):
            del G.nodes[v]["_atom_str"]
        if("_pos" in G.nodes[v]):
            del G.nodes[v]["_pos"]
        if("class" in G.nodes[v]):
            del G.nodes[v]["class"]
    for (u, v, info) in G.edges(data = True):
        if("_bond_str" in G[u][v]):
            del G[u][v]["_bond_str"]
        if("_pos" in G[u][v]):
            del G[u][v]["_pos"]

    # remove unnecessary attributes of H
    for (v, info) in H.nodes(data = True):
        if("_atom_str" in H.nodes[v]):
            del H.nodes[v]["_atom_str"]
        if("_pos" in H.nodes[v]):
            del H.nodes[v]["_pos"]
        if("class" in H.nodes[v]):
            del H.nodes[v]["class"]
    for (u, v, info) in H.edges(data = True):
        if("_bond_str" in H[u][v]):
            del H[u][v]["_bond_str"]
        if("_pos" in H[u][v]):
            del H[u][v]["_pos"]

    # remove unnecessary attributes of full_hydro_G
    for (v, info) in full_hydro_G.nodes(data = True):
        if("_atom_str" in full_hydro_G.nodes[v]):
            del full_hydro_G.nodes[v]["_atom_str"]
        if("_pos" in full_hydro_G.nodes[v]):
            del full_hydro_G.nodes[v]["_pos"]
        if("class" in full_hydro_G.nodes[v]):
            del full_hydro_G.nodes[v]["class"]
    for (u, v, info) in full_hydro_G.edges(data = True):
        if("_bond_str" in full_hydro_G[u][v]):
            del full_hydro_G[u][v]["_bond_str"]
        if("_pos" in full_hydro_G[u][v]):
            del full_hydro_G[u][v]["_pos"]

    # remove unnecessary attributes of full_hydro_H
    for (v, info) in full_hydro_H.nodes(data = True):
        if("_atom_str" in full_hydro_H.nodes[v]):
            del full_hydro_H.nodes[v]["_atom_str"]
        if("_pos" in full_hydro_H.nodes[v]):
            del full_hydro_H.nodes[v]["_pos"]
        if("class" in full_hydro_H.nodes[v]):
            del full_hydro_H.nodes[v]["class"]
    for (u, v, info) in full_hydro_H.edges(data = True):
        if("_bond_str" in full_hydro_H[u][v]):
            del full_hydro_H[u][v]["_bond_str"]
        if("_pos" in full_hydro_H[u][v]):
            del full_hydro_H[u][v]["_pos"]

    # end of function
    return(hydro_balanced, G, H, full_hydro_G, full_hydro_H)




# function: evaluate unmatched hydrogens ---------------------------------------
def get_unmatched_hydrogens(base_G, base_H, full_hydro_G, full_hydro_H):

    # local variables
    k = 0
    min_H = 0
    num_HR = 0
    num_HP = 0
    last_index = 0
    all_maps = []
    extra_HR = []
    extra_HP = []
    indices_HR = []
    indices_HP = []
    base_map_H = []
    hydrogens_R = []
    hydrogens_P = []
    non_H_atoms = []
    all_domains = []
    removable_HR = []
    removable_HP = []
    new_extension = []
    the_extended_maps = []
    names_HR = dict()
    names_HP = dict()
    v = None
    h = None
    each_map = None
    unmatched_hydro_G = None
    unmatched_hydro_H = None

    # obtain mapped non-H vertices
    non_H_atoms = [v for v in base_G.nodes(data = False) if(not base_G.nodes[v]["element"] == "H")]

    # iterate obtaining unmatched hydrogens
    for v in non_H_atoms:

        # get all hydrogens
        hydrogens_R = [h for h in list(full_hydro_G.neighbors(v)) if("h." in str(h))]
        hydrogens_P = [h for h in list(full_hydro_H.neighbors(v)) if("h." in str(h))]
        num_HR = len(hydrogens_R)
        num_HP = len(hydrogens_P)

        # map hydrogens only if not both zero
        if((num_HR > 0) or (num_HP > 0)):

            # case 1: atom lost hydrogens in products
            if((num_HR > 0) and (num_HP == 0)):
                extra_HR = extra_HR + deepcopy(hydrogens_R)
                continue

            # case 2: atom had no hydrogens in reactants
            if((num_HR == 0) and (num_HP > 0)):
                extra_HP = extra_HP + deepcopy(hydrogens_P)
                continue

            # case 3: positive number of hydrogens in both sides
            if((num_HR > 0) and (num_HP > 0)):

                # get removable (matchable) hydrogens if any
                if(num_HR == num_HP):

                    # all hydrogens are removable if all are matchable
                    removable_HR = removable_HR + deepcopy(hydrogens_R)
                    removable_HP = removable_HP + deepcopy(hydrogens_P)

                else:

                    # get indices of hydrogens adjacent with this atom
                    indices_HR = [int((h.split("."))[1]) for h in hydrogens_R]
                    indices_HP = [int((h.split("."))[1]) for h in hydrogens_P]
                    indices_HR.sort()
                    indices_HP.sort()

                    # get minimum number of matchable hydrogens
                    min_H = min([num_HR, num_HP])

                    # get removable hydrogens
                    removable_HR = removable_HR + ["h." + str(k) for k in indices_HR[:min_H]]
                    removable_HP = removable_HP + ["h." + str(k) for k in indices_HP[:min_H]]

                    # save unmatched hydrogens
                    if(num_HR > min_H):
                        extra_HR = extra_HR + ["h." + str(k) for k in indices_HR[min_H:]]
                    if(num_HP > min_H):
                        extra_HP = extra_HP + ["h." + str(k) for k in indices_HP[min_H:]]

    # get copies and remove matchable hydrogens
    # * by construction at this point len(extra_HR) == len(extra_HP)
    unmatched_hydro_G = deepcopy(full_hydro_G)
    unmatched_hydro_G.remove_nodes_from(removable_HR)
    unmatched_hydro_H = deepcopy(full_hydro_H)
    unmatched_hydro_H.remove_nodes_from(removable_HP)

    # end of function
    return(extra_HR, extra_HP, unmatched_hydro_G, unmatched_hydro_H)




# function: get the ITS graph of a reaction and an AAM -------------------------
def build_ITS_graph(nx_G, nx_H, atom_map):

    # local variables
    null_symbol = "*"
    edges_G = []
    edges_H = []
    new_label = ()
    info_dict = dict()
    forward_map = dict()
    inverse_map = dict()
    u = None
    v = None
    nx_ITS = None

    # get inverse AAM
    for (x, y) in atom_map:
        forward_map[x] = y
        inverse_map[y] = x

    # canonical ITS
    nx_ITS = nx.Graph()

    # add node labels to ITS graph
    for (v, info_dict) in list(nx_G.nodes(data = True)):
        # new node label for ITS graph
        new_label = deepcopy((info_dict, info_dict))
        # add node with compressed attributes
        nx_ITS.add_node(v, GMNL_ITS = new_label)

    # copy original edges
    edges_G = list(nx_G.edges(data = False))
    edges_H = list(nx_H.edges(data = False))

    # add edges from G
    for (u, v) in edges_G:

        # check if also an edge in H and add the respective edge
        if((forward_map[u], forward_map[v]) in edges_H):
            new_label = (nx_G[u][v], nx_H[forward_map[u]][forward_map[v]])
            nx_ITS.add_edge(u, v, GMEL_ITS = deepcopy(new_label))
        else:
            if((forward_map[v], forward_map[u]) in edges_H):
                new_label = (nx_G[u][v], nx_H[forward_map[v]][forward_map[u]])
                nx_ITS.add_edge(u, v, GMEL_ITS = deepcopy(new_label))
            else:
                new_label = (nx_G[u][v], null_symbol)
                nx_ITS.add_edge(u, v, GMEL_ITS = deepcopy(new_label))

    # add edges from H
    for (u, v) in edges_H:

        # check if also an edge in G and add the respective edge
        if((inverse_map[u], inverse_map[v]) in edges_G):
            continue
        else:
            if((inverse_map[v], inverse_map[u]) in edges_G):
                continue
            else:
                new_label = (null_symbol, nx_H[u][v])
                nx_ITS.add_edge(inverse_map[u], inverse_map[v], GMEL_ITS = deepcopy(new_label))

    # end of function
    return(nx_ITS)




# buid copy with inteligent labeling and run isomorphism -----------------------
def inteligent_labeling_isomorphism(G, H, anchor, algorithm):

    # local variables
    are_isomorphic = False
    anchor_counter = 0
    nodes_G = []
    nodes_H = []
    edges_G = []
    edges_H = []
    anchor_G = []
    anchor_H = []
    isomorphisms = []
    reaction_edges_G = []
    reaction_edges_H = []
    empty_dict = dict()
    anchor_labels_G = dict()
    anchor_labels_H = dict()
    rc_G = nx.Graph()
    rc_H = nx.Graph()
    remainder_G = nx.Graph()
    remainder_H = nx.Graph()
    GM = None
    inteligent_node_match = None
    inteligent_edge_match = None

    # label for anchor
    for (a, b) in anchor:
        anchor_counter = anchor_counter + 1
        anchor_labels_G[a] = anchor_counter
        anchor_labels_H[b] = anchor_counter
        anchor_G.append(a)
        anchor_H.append(b)

    # relabel vertices
    nodes_G = list(G.nodes(data = True))
    for (node, info)  in nodes_G:
        if(node in anchor_G):
            remainder_G.add_node(node, inteligent_label = (anchor_labels_G[node], info))
        else:
            remainder_G.add_node(node, inteligent_label = (0, info))
    nodes_H = list(H.nodes(data = True))
    for (node, info)  in nodes_H:
        if(node in anchor_H):
            remainder_H.add_node(node, inteligent_label = (anchor_labels_H[node], info))
        else:
            remainder_H.add_node(node, inteligent_label = (0, info))

    # get edges in reaction center
    rc_G = deepcopy(G)
    rc_H = deepcopy(H)
    rc_G.remove_nodes_from([v for (v, info) in nodes_G if(v not in anchor_G)])
    rc_H.remove_nodes_from([v for (v, info) in nodes_H if(v not in anchor_H)])
    reaction_edges_G = list(rc_G.edges(data = True))
    reaction_edges_H = list(rc_H.edges(data = True))

    # add edges
    remainder_G.add_edges_from(list(G.edges(data = True)))
    remainder_H.add_edges_from(list(H.edges(data = True)))
    remainder_G.remove_edges_from(reaction_edges_G)
    remainder_H.remove_edges_from(reaction_edges_H)

    # run specified algorithm
    if(algorithm == "gm_iso"):
        isomorphisms, are_isomorphic = gm.search_isomorphisms(remainder_G, remainder_H, node_labels = True, edge_labels = True, all_isomorphisms = False)

    if(algorithm == "nx_iso"):
        inteligent_node_match = generic_node_match("inteligent_label", (0, empty_dict), eq)
        inteligent_edge_match = generic_edge_match("bond_type", "b0", eq)
        are_isomorphic   = nx.is_isomorphic(remainder_G, remainder_H, node_match = inteligent_node_match, edge_match = inteligent_edge_match)
        isomorphisms = []

    # end of function
    return(isomorphisms, are_isomorphic)




# analysis #####################################################################




# output data holders
times_method_1 = []
times_method_2 = []
order_of_base_ITS = []
feasible_reactions = []
time_for_automorphisms = []
number_of_automorphisms = []
number_of_equivalence_classes = []
number_of_unmatched_hydrogens = []




# initialize data holders
feasible = False
uneven_aam = []
not_readable = []
unbalanced_hydro = []
no_unmatched_hydrogens = []




# load all reactions from Long's pre-analysis
list_of_reactions = load_from_pickle("hydrogen.pkl.gz")




# node and edge match for isomorphism tests
GM_node_match = generic_node_match("GMNL_ITS", "e", eq)
GM_edge_match = generic_edge_match("GMEL_ITS", ("*", "*"), eq)




# filter reactions
for reaction_entry in list_of_reactions:

    # reinitialize control
    feasible = True

    # print all atom maps
    if(print_all):
        print("######################################")
        print(reaction_entry["R-id"])
        print(reaction_entry["aam"])
        print("######################################")

    # test eveness and readibility
    try:

        # test eveness: all atoms in the reactants are mapped to all atoms
        # in the products and viceversa, meaning that...
        # - all atoms in the molecules of the base reaction are being mapped
        # - all atom-map indices appear in both sides of the reactions
        # - both sides have exactly the same set of indices
        even_aam = check_eveness(reaction_entry["aam"])

        # save if aam uneven
        if(not even_aam):

            # save as uneven
            feasible = False
            uneven_aam.append(reaction_entry["R-id"])
            # print("--------------------------------------")
            # print("- uneven")
            # print(reaction_entry["R-id"])
            # print(reaction_entry["aam"])
            # print("--------------------------------------")

    except:

        # save if aam is not readable
        feasible = False
        not_readable.append(reaction_entry["R-id"])
        # print("--------------------------------------")
        # print("- unreadable")
        # print(reaction_entry["R-id"])
        # print(reaction_entry["aam"])
        # print("--------------------------------------")

    # skip if unfeasible
    if(not feasible):
        continue

    # get reaction smiles for hydrogen testing
    reaction_smiles = reaction_entry["aam"]

    # get base and full-hydrogen graph representations of the SMILES
    hydro_balanced, G, H, full_hydro_G, full_hydro_H = get_graph_representations(reaction_smiles)

    # analyze reaction
    if(not hydro_balanced):

        # store reaction as not balanced up-to hydrogens
        unbalanced_hydro.append(reaction_entry["R-id"])
        # print("--------------------------------------")
        # print("- unbalanced hydrogens")
        # print(reaction_entry["R-id"])
        # print(reaction_entry["aam"])
        # print("--------------------------------------")

    else:

        # get unmatched hydrogens and reduced hydro-versions
        extra_hydro_G, extra_hydro_H, unmatched_hydro_G, unmatched_hydro_H = get_unmatched_hydrogens(G, H, full_hydro_G, full_hydro_H)

        # from SMILES to graphs with Long's functions
        reactants_side, products_side = rsmi_to_graph(rsmi = reaction_smiles,
                                                      drop_non_aam = True,
                                                      sanitize = True,
                                                      use_index_as_atom_map = True,
                                                      node_attrs = ["element", "aromatic", "hcount", "charge"],
                                                      edge_attrs = ["order"])

        # get hcount for comparison
        hcount_change = check_hcount_change(reactants_side, products_side)
        if((not hcount_change == len(extra_hydro_G)) or (not hcount_change == len(extra_hydro_H))):
            print("######################################")
            print("BIG PROBLEM")
            print("######################################")
        if(hcount_change == 0):
            no_unmatched_hydrogens.append(reaction_entry["R-id"])
            continue
        print(hcount_change, len(extra_hydro_G), len(extra_hydro_H))

        # save number of unmatched_hydrogens
        number_of_unmatched_hydrogens.append(len(extra_hydro_G))




        # report running
        print("- analyzing: ", reaction_entry["R-id"])

        # save reaction as feasible
        feasible_reactions.append(reaction_entry["R-id"])

        # get base atom-map
        base_aam = [(v, v) for v in G.nodes(data = False)]

        # save order of base ITS graph
        order_of_base_ITS.append(len(base_aam))




        # ------------------------------------------------------------------------
        # METHOD 1: the usual way by testing the isomorphism of the ITS graphs
        # ------------------------------------------------------------------------
        last_class = 0
        all_classes_method_1 = []
        representatives_method_1 = dict()
        equivalence_classes_method_1 = dict()

        # initialize timer
        initial_time = time.time()

        # classify hydrogen extensions
        for each_permutation in permutations(extra_hydro_G):

            # get extended aam
            extended_aam = deepcopy(base_aam) + list(zip(each_permutation, extra_hydro_H))

            # get hydro ITS graph
            hydro_ITS = build_ITS_graph(unmatched_hydro_G, unmatched_hydro_H, extended_aam)

            # restart control variable
            found_equivalent = False

            # compare with classified ITSs (if any yet)
            for each_class in all_classes_method_1:

                # test isomorphism
                are_isomorphic = nx.is_isomorphic(hydro_ITS,
                                                  representatives_method_1[each_class][0],
                                                  node_match = GM_node_match,
                                                  edge_match = GM_edge_match)

                # update class
                if(are_isomorphic):
                    # store permutation in the corresponding class
                    equivalence_classes_method_1[each_class].add(deepcopy(each_permutation))
                    # update and finish loop
                    found_equivalent = True
                    break

            # if equivalent was not found create new equivalence class
            if(not found_equivalent):
                # set class number and save it
                last_class = last_class + 1
                all_classes_method_1.append(last_class)
                # save (unique) representative and the first permutation
                representatives_method_1[last_class] = [deepcopy(hydro_ITS)]
                equivalence_classes_method_1[last_class] = {deepcopy(each_permutation)}

        # end timer
        final_time = time.time()

        # save running time - method 1
        times_method_1.append(final_time - initial_time)




        # ------------------------------------------------------------------------
        # METHOD 2: using the automorphisms of the base ITS graph and co-extension
        # ------------------------------------------------------------------------
        last_class = 0
        base_automorphisms = []
        all_classes_method_2 = []
        representatives_method_2 = dict()
        equivalence_classes_method_2 = dict()

        # initialize timer
        initial_time = time.time()

        # get base ITS graph
        base_ITS = build_ITS_graph(G, H, base_aam)

        # initialize timer automorphisms
        initial_time_auto = time.time()

        # get automorphisms of base ITS graph
        GM = isomorphism.GraphMatcher(base_ITS,
                                      base_ITS,
                                      node_match = GM_node_match,
                                      edge_match = GM_edge_match)
        base_automorphisms_dicts = list(GM.isomorphisms_iter())
        base_automorphisms = []
        for each_automorphism in base_automorphisms_dicts:
            base_automorphisms.append(list(each_automorphism.items()))

        # end timer automorphisms
        final_time_auto = time.time()

        # save running time of automorphism search
        time_for_automorphisms.append(final_time_auto - initial_time_auto)

        # classify hydrogen extensions
        for each_permutation in permutations(extra_hydro_G):

            # get extended aam
            extended_aam = deepcopy(base_aam) + list(zip(each_permutation, extra_hydro_H))

            # get hydro ITS graph
            hydro_ITS = build_ITS_graph(unmatched_hydro_G, unmatched_hydro_H, extended_aam)

            # restart control variable
            found_equivalent = False

            # compare with classified ITSs (if any yet)
            for each_class in all_classes_method_2:

                # test isomorphism via co-extension
                for each_automorphism in base_automorphisms:

                    # test isomorphism
                    isomorphisms, are_isomorphic = inteligent_labeling_isomorphism(hydro_ITS,
                                                                                   representatives_method_2[each_class][0],
                                                                                   each_automorphism,
                                                                                   "nx_iso")

                    # update class
                    if(are_isomorphic):
                        # store permutation in the corresponding class
                        equivalence_classes_method_2[each_class].add(deepcopy(each_permutation))
                        # update and finish loop
                        found_equivalent = True
                        break

                # break external loop
                if(found_equivalent):
                    break

            # if equivalent was not found create new equivalence class
            if(not found_equivalent):
                # set class number and save it
                last_class = last_class + 1
                all_classes_method_2.append(last_class)
                # save (unique) representative and the first permutation
                representatives_method_2[last_class] = [deepcopy(hydro_ITS)]
                equivalence_classes_method_2[last_class] = {deepcopy(each_permutation)}

        # end timer
        final_time = time.time()

        # save running time and automorphisms - method 2
        times_method_2.append(final_time - initial_time)
        number_of_automorphisms.append(len(base_automorphisms))




        # ------------------------------------------------------------------------
        # Compare equivalence classes obtained with both methods for santiy check
        # ------------------------------------------------------------------------

        # check that the same amount of classes was obtained
        if(not len(all_classes_method_1) == len(all_classes_method_2)):
            print("* something is wrong...")

        # initialize data holders
        found = False
        matched_class = 0
        correspondence = []
        unmatched_method_1 = []
        unmatched_method_2 = []
        classes_method_1 = deepcopy(all_classes_method_1)
        classes_method_2 = deepcopy(all_classes_method_2)

        # compare classes
        for each_class_1 in classes_method_1:

            # reinitialize control
            found = False

            # loop through (remaining) classes
            for each_class_2 in classes_method_2:
                if(equivalence_classes_method_1[each_class_1] == equivalence_classes_method_2[each_class_2]):
                    matched_class = each_class_2
                    correspondence.append((each_class_1, each_class_2))
                    found = True
                    break

            # remove matched class
            if(found):
                classes_method_2.remove(matched_class)

        # report unmatched classes
        unmatched_method_1 = deepcopy(all_classes_method_1)
        unmatched_method_2 = deepcopy(all_classes_method_2)
        for (each_class_1, each_class_2) in correspondence:
            unmatched_method_1.remove(each_class_1)
            unmatched_method_2.remove(each_class_2)
        if(len(unmatched_method_1) > 0):
            print("- unmatched classes from method 1: ", unmatched_method_1)
        if(len(unmatched_method_2) > 0):
            print("- unmatched classes from method 2: ", unmatched_method_2)
        if((len(unmatched_method_1) == 0) and (len(unmatched_method_2) == 0)):
            print("- all equivalence classes of extensions where correctly matched")
            print(correspondence)
        print("---------------------------------------------------------")

        # save number of classes
        number_of_equivalence_classes.append(len(all_classes_method_1))




# ------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------

# initialize data holders
X = []
Y = []
C = []
E = []
Z = []
W = []

df_dict = dict()
df_dict["reaction"] = []
df_dict["hydrogens"] = []
df_dict["time"] = []
df_dict["method"] = []

time_dict = dict()
time_dict["time_proportion"] = []

average_m1 = dict()
average_m1[2] = []
average_m1[3] = []
average_m1[4] = []
average_m1[5] = []

average_m2 = dict()
average_m2[2] = []
average_m2[3] = []
average_m2[4] = []
average_m2[5] = []

# get data to plot
for i in range(len(feasible_reactions)):

    # proportion for X axis
    X.append(math.log(times_method_1[i] * 1000, 10))
    Y.append(math.log(times_method_2[i] * 1000, 10))
    C.append(math.log(number_of_automorphisms[i], 10))
    E.append(number_of_equivalence_classes[i])
    Z.append(math.log(number_of_automorphisms[i], 10))
    W.append(math.log(number_of_equivalence_classes[i], 10))

    # print results for preview
    print("- reaction: ", feasible_reactions[i])
    print("- order: ", order_of_base_ITS[i])
    print("- unmatched hydrogens: ", number_of_unmatched_hydrogens[i])
    print("- hydro proportion: ", number_of_unmatched_hydrogens[i] * 100 / (order_of_base_ITS[i] + number_of_unmatched_hydrogens[i]))
    print("- equivalence classes: ", number_of_equivalence_classes[i])
    print("- automorphisms: ", number_of_automorphisms[i])
    print("- times 1: ", times_method_1[i], ", times 2: ", times_method_2[i])
    print("- times auto: ", time_for_automorphisms[i] * 100 / times_method_2[i])
    print("---------------------------------------------------------")

    # save data frame info for method 1
    df_dict["reaction"].append("entry-m1:" + feasible_reactions[i])
    df_dict["hydrogens"].append(number_of_unmatched_hydrogens[i])
    df_dict["time"].append(math.log(times_method_1[i]*1000, 10))
    df_dict["method"].append("m1")

    # save data frame info for method 2
    df_dict["reaction"].append("entry-m2:" +feasible_reactions[i])
    df_dict["hydrogens"].append(number_of_unmatched_hydrogens[i])
    df_dict["time"].append(math.log(times_method_2[i]*1000, 10))
    df_dict["method"].append("m2")

    # save time taken for automorpisms as proportion of time for method 2
    time_dict["time_proportion"].append(time_for_automorphisms[i] * 100 / times_method_2[i])

    # get least for avergae
    average_m1[number_of_unmatched_hydrogens[i]].append(times_method_1[i]*1000)
    average_m2[number_of_unmatched_hydrogens[i]].append(times_method_2[i]*1000)




# print total of feasible reactions and with each amount of hydrogens
print("\n")
print("- Total of feasible reactions: ", len(feasible_reactions))
print("\n")
print("- number mA, 2H: ", len(average_m1[2]))
print("- number mA, 3H: ", len(average_m1[3]))
print("- number mA, 4H: ", len(average_m1[4]))
print("- number mA, 5H: ", len(average_m1[5]))
print("---------------------------------")
print("- number mB, 2H: ", len(average_m2[2]))
print("- number mB, 3H: ", len(average_m2[3]))
print("- number mB, 4H: ", len(average_m2[4]))
print("- number mB, 5H: ", len(average_m2[5]))
print("\n")

# print average running time
print("- average and std mA, 2H: ", np.mean(average_m1[2]), np.std(average_m1[2]))
print("- average and std mA, 3H: ", np.mean(average_m1[3]), np.std(average_m1[3]))
print("- average and std mA, 4H: ", np.mean(average_m1[4]), np.std(average_m1[4]))
print("- average and std mA, 5H: ", np.mean(average_m1[5]), np.std(average_m1[5]))
print("---------------------------------")
print("- average and std mB, 2H: ", np.mean(average_m2[2]), np.std(average_m2[2]))
print("- average and std mB, 3H: ", np.mean(average_m2[3]), np.std(average_m2[3]))
print("- average and std mB, 4H: ", np.mean(average_m2[4]), np.std(average_m2[4]))
print("- average and std mB, 5H: ", np.mean(average_m2[5]), np.std(average_m2[5]))
print("\n")




# dump data
output_file = open("3_results.pickle", "wb")
pickle.dump([times_method_1,
             times_method_2,
             order_of_base_ITS,
             feasible_reactions,
             time_for_automorphisms,
             number_of_automorphisms,
             number_of_equivalence_classes,
             number_of_unmatched_hydrogens],
            output_file)
output_file.close()




# make scatter plot
scp = plt.scatter(X, Y, c = C, alpha = 0.75, cmap = cmr.get_sub_cmap("viridis", 0.15, 0.95))
plt.plot([0, 4.5], [0, 4.5], color = "grey", linestyle = "-.", linewidth = 1)

# figure attributes
plt.xlabel(r"Log$_{10}$ of running time of method A [ms]", color = "w", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method B [ms]", color = "w", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)

# plot color bar of scatter plots
cbar = plt.colorbar(scp)
cbar.solids.set(alpha = 1)
cbar.ax.set_ylabel(r"Log$_{10}$ of number of automorphisms", labelpad = 25, rotation = 270, size = 10)
cbar.ax.tick_params(labelsize = 10)
plt.tight_layout()

# save figure
file_name = "3_colorbar.pdf"
plt.savefig(file_name)
plt.close()




# make scatter plot
scp = plt.scatter(X, Y, c = C, alpha = 0.75, cmap = cmr.get_sub_cmap("viridis", 0.15, 0.95))
plt.plot([0, 4.5], [0, 4.5], color = "grey", linestyle = "-.", linewidth = 1)

# figure attributes
plt.xlabel(r"Log$_{10}$ of running time of method A [ms]", color = "k", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method B [ms]", color = "k", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)

# plot color bar of scatter plots
cbar = plt.colorbar(scp)
cbar.solids.set(alpha = 1)
cbar.ax.set_ylabel(r"Log$_{10}$ of number of automorphisms", labelpad = 25, rotation = 270, size = 10)
cbar.ax.tick_params(labelsize = 10)
plt.tight_layout()

# save figure
file_name = "3_anchored_iso_gm.pdf"
plt.savefig(file_name)
plt.close()




# make scatter plot
scp = plt.scatter(E, X, c = "tab:blue", alpha = 0.5)

# figure attributes
plt.xlabel("Number of equivalence classes", color = "k", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method A [ms]", color = "k", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylim(-0.5, 4.5)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
plt.tight_layout()

# save figure
file_name = "3_time_vs_classes_A.pdf"
plt.savefig(file_name)
plt.close()




# make scatter plot
scp = plt.scatter(E, Y, c = "tab:red", alpha = 0.5)

# figure attributes
plt.xlabel("Number of equivalence classes", color = "w", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method B [ms]", color = "k", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylim(-0.5, 4.5)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)
plt.tight_layout()

# save figure
file_name = "3_time_vs_classes_B.pdf"
plt.savefig(file_name)
plt.close()




# make scatter plot
scp = plt.scatter(X, Y, c = E, alpha = 0.6, cmap = cmr.get_sub_cmap("plasma", 0.15, 0.95))
plt.plot([0, 4.5], [0, 4.5], color = "grey", linestyle = "-.", linewidth = 1)

# figure attributes
plt.xlabel(r"Log$_{10}$ of running time of method A [ms]", color = "w", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method B [ms]", color = "w", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)

# plot color bar of scatter plots
cbar = plt.colorbar(scp)
cbar.solids.set(alpha = 1)
cbar.ax.set_ylabel("Number of equivalence classes of ITS graphs", labelpad = 25, rotation = 270, size = 10)
cbar.ax.tick_params(labelsize = 10)
plt.tight_layout()

# save figure
file_name = "3_classes_without_labels.pdf"
plt.savefig(file_name)
plt.close()




# make scatter plot
scp = plt.scatter(X, Y, c = E, alpha = 0.6, cmap = cmr.get_sub_cmap("plasma", 0.15, 0.95))
plt.plot([0, 4.5], [0, 4.5], color = "grey", linestyle = "-.", linewidth = 1)

# figure attributes
plt.xlabel(r"Log$_{10}$ of running time of method A [ms]", color = "k", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time of method B [ms]", color = "k", labelpad = 10, size = 10)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.4)

# plot color bar of scatter plots
cbar = plt.colorbar(scp)
cbar.solids.set(alpha = 1)
cbar.ax.set_ylabel("Number of equivalence classes of ITS graphs", labelpad = 25, rotation = 270, size = 10)
cbar.ax.tick_params(labelsize = 10)
plt.tight_layout()

# save figure
file_name = "3_classes_with_labels.pdf"
plt.savefig(file_name)
plt.close()




# load the dataset
data_frame = pd.DataFrame(df_dict)

# create figure and axes
fig, ax = plt.subplots()

# draw violins for each sex
sns.violinplot(data = data_frame, x = "hydrogens", y = "time", hue = "method", split = True, gap = 0.1, inner = "quart", ax = ax)

# set attributes
for violin in ax.collections:
    violin.set_alpha(0.75)
plt.legend([],[], frameon = False)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.5)
plt.xlabel("Number of unmatched hydrogens", color = "k", labelpad = 10, size = 10)
plt.ylabel(r"Log$_{10}$ of running time [ms]", color = "k", labelpad = 10, size = 10)
plt.tight_layout()

# save figure
file_name = "3_violins.pdf"
plt.savefig(file_name)
plt.close()



# distribution of time for calculating automorphisms
time_data_frame = pd.DataFrame(time_dict)
sns.displot(time_data_frame, x = "time_proportion")
plt.legend([],[], frameon = False)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.5)
plt.xlabel("Proportion of time for eval. of automorphism group [%]", color = "k", labelpad = 10, size = 10)
plt.ylabel("Count of reactions", color = "k", labelpad = 10, size = 10)
plt.xlim(-3, 50)
plt.ylim(0, 60)
plt.tight_layout()

# save figure
file_name = "3_distribution.pdf"
plt.savefig(file_name)
plt.close()




################################################################################
################################################################################
