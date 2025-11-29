################################################################################
#                                                                              #
#  - Analysis of Hydrogen Extensions                                           #
#                                                                              #
#  - Random ITS graphs                                                         #
#                                                                              #
#  - Made by Marcos Laffitte - Github @MarcosLaffitte                          #
#                                                                              #
#  - 2 - Anchor isomorphism with Relabeling and GranMapache                    #
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




# base graphs for ITS with 1 automorphism --------------------------------------

# reactants graph and nodes
G = nx.Graph()
for k in range(1, 26+1):
    G.add_node(k, element = "C", aromatic = False, hcount = 0, charge = 0)
# reaction center edges
G.add_edge(1, 24, order = 1.0)
G.add_edge(24, 25, order = 1.0)
G.add_edge(21, 22, order = 1.0)
G.add_edge(23, 26, order = 1.0)
G.add_edge(26, 11, order = 1.0)
# left path
for k in range(1, 9+1):
    G.add_edge(k, k+1, order = 1.0)
# right path
for k in range(11, 19+1):
    G.add_edge(k, k+1, order = 1.0)

# reactants graph and nodes
H = nx.Graph()
for k in range(1, 26+1):
    H.add_node(k, element = "C", aromatic = False, hcount = 0, charge = 0)
# reaction center edges
H.add_edge(1, 24, order = 1.0)
H.add_edge(24, 21, order = 1.0)
H.add_edge(22, 23, order = 1.0)
H.add_edge(25, 26, order = 1.0)
H.add_edge(26, 11, order = 1.0)
# left path
for k in range(1, 9+1):
    H.add_edge(k, k+1, order = 1.0)
# right path
for k in range(11, 19+1):
    H.add_edge(k, k+1, order = 1.0)




# base atom-map for ITS with 1 automorphism ------------------------------------
base_aam = [(v, v) for v in G.nodes(data = False)]




# functions ####################################################################




# function: evaluate unmatched hydrogens ---------------------------------------
def get_unmatched_hydrogens(hydro_in_reactants, hydro_in_products):

    # global variables
    global G
    global H

    # local variables
    k = 0
    hydrogen_name = ""
    extra_HR = []
    extra_HP = []
    unmatched_hydro_G = None
    unmatched_hydro_H = None

    # copy global base graphs
    unmatched_hydro_G = deepcopy(G)
    unmatched_hydro_H = deepcopy(H)

    # add hydrogens to G and H
    for k in range(len(hydro_in_reactants)):

        # make new name
        hydrogen_name = "h." + str(k+1)

        # add hydrogen
        unmatched_hydro_G.add_node(hydrogen_name, element = "H", aromatic = False, hcount = 0, charge = 0)
        unmatched_hydro_H.add_node(hydrogen_name, element = "H", aromatic = False, hcount = 0, charge = 0)

        # add bonds
        unmatched_hydro_G.add_edge(hydro_in_reactants[k], hydrogen_name, order = 1.0)
        unmatched_hydro_H.add_edge(hydro_in_products[k], hydrogen_name, order = 1.0)

        # save hydrogens
        extra_HR.append(hydrogen_name)
        extra_HP.append(hydrogen_name)

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
        empty_dict = dict()
        inteligent_node_match = generic_node_match("inteligent_label", (0, empty_dict), eq)
        inteligent_edge_match = generic_edge_match("bond_type", "b0", eq)
        GM = isomorphism.GraphMatcher(remainder_G, remainder_H, node_match = inteligent_node_match, edge_match = inteligent_edge_match)
        iter_iso = GM.isomorphisms_iter()
        iter_iso_first =  next(iter_iso)
        isomorphisms = [iter_iso_first]

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




# load all reactions from Long's pre-analysis
input_file = open("random_reactions.pkl", "rb")
list_of_reactions = pickle.load(input_file)
input_file.close()




# filter reactions
for (hydro_in_reactants, hydro_in_products) in list_of_reactions:

    # get unmatched hydrogens and reduced hydro-versions
    extra_hydro_G, extra_hydro_H, unmatched_hydro_G, unmatched_hydro_H = get_unmatched_hydrogens(hydro_in_reactants, hydro_in_products)

    # save number of unmatched_hydrogens
    number_of_unmatched_hydrogens.append(len(extra_hydro_G))

    # report running
    print("- analyzing: ", hydro_in_reactants, hydro_in_products)
    print(len(extra_hydro_G), len(extra_hydro_H))

    # save reaction as feasible
    feasible_reactions.append((hydro_in_reactants, hydro_in_products))

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
            isomorphisms, are_isomorphic = gm.search_isomorphisms(nx_G = hydro_ITS,
                                                                  nx_H = representatives_method_1[each_class][0],
                                                                  node_labels = True,
                                                                  edge_labels = True,
                                                                  all_isomorphisms = False)

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
    base_automorphisms, dummy = gm.search_isomorphisms(nx_G = base_ITS,
                                                       nx_H = base_ITS,
                                                       node_labels = True,
                                                       edge_labels = True,
                                                       all_isomorphisms = True)

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
                                                                               "gm_iso")

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

df_dict_A = dict()
df_dict_A["reaction"] = []
df_dict_A["hydrogens"] = []
df_dict_A["time"] = []
df_dict_A["method"] = []

df_dict_B = dict()
df_dict_B["reaction"] = []
df_dict_B["hydrogens"] = []
df_dict_B["time"] = []
df_dict_B["method"] = []

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

    # save time taken for automorpisms as proportion of time for method 2
    time_dict["time_proportion"].append(time_for_automorphisms[i] * 100 / times_method_2[i])

    # identifiers
    h_str_G = [str(num) for num in feasible_reactions[i][0]]
    h_str_H = [str(num) for num in feasible_reactions[i][1]]

    # save entries for smaller and bigger hue values
    bigger = [4, 5]
    # bigger = [3, 4]
    if(number_of_unmatched_hydrogens[i] in bigger):

        # save data frame info for method 1
        df_dict_B["reaction"].append("entry-m1:" + ",".join(h_str_G) + ">" + ",".join(h_str_H))
        df_dict_B["hydrogens"].append(number_of_unmatched_hydrogens[i])
        df_dict_B["time"].append(math.log(times_method_1[i]*1000, 10))
        df_dict_B["method"].append("m1")

        # save data frame info for method 2
        df_dict_B["reaction"].append("entry-m2:" + ",".join(h_str_G) + ">" + ",".join(h_str_H))
        df_dict_B["hydrogens"].append(number_of_unmatched_hydrogens[i])
        df_dict_B["time"].append(math.log(times_method_2[i]*1000, 10))
        df_dict_B["method"].append("m2")

    else:

        # save data frame info for method 1
        df_dict_A["reaction"].append("entry-m1:" + ",".join(h_str_G) + ">" + ",".join(h_str_H))
        df_dict_A["hydrogens"].append(number_of_unmatched_hydrogens[i])
        df_dict_A["time"].append(math.log(times_method_1[i]*1000, 10))
        df_dict_A["method"].append("m1")

        # save data frame info for method 2
        df_dict_A["reaction"].append("entry-m2:" + ",".join(h_str_G) + ">" + ",".join(h_str_H))
        df_dict_A["hydrogens"].append(number_of_unmatched_hydrogens[i])
        df_dict_A["time"].append(math.log(times_method_2[i]*1000, 10))
        df_dict_A["method"].append("m2")

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
output_file = open("2_results_random.pickle", "wb")
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




# Load the dataset
data_frame_A = pd.DataFrame(df_dict_A)

# create figure and axes
fig, ax = plt.subplots()

# draw violins for each sex
sns.violinplot(data = data_frame_A, x = "hydrogens", y = "time", hue = "method", split = True, gap = 0.1, inner = "quart", ax = ax)

# set attributes
for violin in ax.collections:
    violin.set_alpha(0.75)
plt.legend([],[], frameon = False)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.5)
plt.xlabel("Number of unmatched hydrogens", color = "k", labelpad = 10, size = 12)
plt.ylabel(r"Log$_{10}$ of running time [ms]", color = "k", labelpad = 10, size = 12)
plt.ylim(-0.05, 2.25)
plt.tight_layout()

# save figure
file_name = "2_violins_A.pdf"
plt.savefig(file_name)
plt.close()




# Load the dataset
data_frame_B = pd.DataFrame(df_dict_B)

# create figure and axes
fig, ax = plt.subplots()

# draw violins for each sex
sns.violinplot(data = data_frame_B, x = "hydrogens", y = "time", hue = "method", split = True, gap = 0.1, inner = "quart", ax = ax)

# set attributes
for violin in ax.collections:
    violin.set_alpha(0.75)
plt.legend([],[], frameon = False)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.5)
plt.xlabel("Number of unmatched hydrogens", color = "k", labelpad = 10, size = 12)
plt.ylabel(r"Log$_{10}$ of running time [ms]", color = "k", labelpad = 10, size = 12)
plt.ylim(2.25, 4.30)
plt.tight_layout()

# save figure
file_name = "2_violins_B.pdf"
plt.savefig(file_name)
plt.close()




# distribution of time for calculating automorphisms
time_data_frame = pd.DataFrame(time_dict)
sns.displot(time_data_frame, x = "time_proportion")
plt.legend([],[], frameon = False)
plt.grid(visible = True, axis = "both", color = "grey", linestyle = "--", linewidth = 0.5)
plt.xlabel("Proportion of time for eval. of automorphism group [%]", color = "k", labelpad = 10, size = 10)
plt.ylabel("Count of reactions", color = "k", labelpad = 10, size = 12)
plt.xlim(-3, 70)
plt.ylim(0, 250)
plt.tight_layout()

# save figure
file_name = "2_distribution_random.pdf"
plt.savefig(file_name)
plt.close()




################################################################################
################################################################################
