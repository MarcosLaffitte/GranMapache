################################################################################
#                                                                              #
#  - Analysis of Random Graphs for Extension of Partial Atonm-to-Atom Maps     #
#                                                                              #
#  - Extensions                                                                #
#                                                                              #
#  - for WABI 2025                                                             #
#                                                                              #
#  - Made by Marcos Laffitte - Github @MarcosLaffitte                          #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import math
import time
import random
import pickle
from copy import deepcopy


# additional dependencies ------------------------------------------------------
import numpy
import cython
import networkx as nx
import matplotlib.pyplot as plt


# networkx isomorphism with node and edge match --------------------------------
from operator import eq
from networkx.algorithms import isomorphism
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match


# custom dependencies ----------------------------------------------------------
import gmapache as gm


# parameters ###################################################################


# input ------------------------------------------------------------------------
in_G = None
in_H = None
in_anchor = None


# parameters of remainder graphs -----------------------------------------------
nodes_remainders = [100, 125, 150, 175, 200]
# NOTE: density of at least 2 is needed for connectedness with 100 nodes,
# but with density exactly 2 and 125 nodes it seems way too hard to find a random connected graph,
# and density of at most 97 is also needed for possibility of construction with 100 vertices,
# thus the 3-97 range for density was selected, which seems to work well
densities_remainders = list(range(3, 97 + 1))
density_breakpoints = [10, 20, 30, 40, 50, 60, 70, 80, 90, 97]
density_counter = 0


# general variables ------------------------------------------------------------
tail = "".join(10*[" "])
reaction_counter = 0
reactions_by_order = []
each_list_of_reactions_with_fixed_density = []


# output -----------------------------------------------------------------------
out_time_gm_extender = 0
out_time_gm_isomorphism = 0
out_time_nx_isomorphism = 0
out_time_nx_ismags_algo = 0
times_by_order = []
times_by_fixed_density = []


# functions ####################################################################


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

    # get endges in reaction center
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
        are_isomorphic = nx.is_isomorphic(remainder_G, remainder_H, node_match = inteligent_node_match, edge_match = inteligent_edge_match)
    if(algorithm == "ismags"):
        empty_dict = dict()
        inteligent_node_match = generic_node_match("inteligent_label", (0, empty_dict), eq)
        inteligent_edge_match = generic_edge_match("bond_type", "b0", eq)
        ismags = nx.isomorphism.ISMAGS(remainder_G, remainder_H, node_match = inteligent_node_match, edge_match = inteligent_edge_match)
        are_isomorphic = ismags.is_isomorphic(symmetry = True)

    # end of function
    return(isomorphisms, are_isomorphic)


# analysis #####################################################################


# task message
print("\n")
print("######################################################################")
print("> Analyzing random graphs")
print("######################################################################")
print("\n")


# get number of nodes
for each_order in nodes_remainders:

    # reinitialize variables
    times_by_order = []
    density_counter = 0

    # get density breakpoint
    for each_density_breakpoint in density_breakpoints:

        # open the file with graphs of each order
        file_name = "Random_Graphs_" + str(each_order) + "_nodes_and_upper_density_" + str(each_density_breakpoint) + ".pkl"
        print("* Loading graphs with " + str(each_order) + " nodes and upper density " + str(each_density_breakpoint))
        with open(file_name, "rb") as handle:
            reactions_by_order = pickle.load(handle)
        print("\n")

        # open each sublist of graphs by density
        for each_list_of_reactions_with_fixed_density in reactions_by_order:

            # task message
            print("----------------------------------------------------------------------")
            print("* extending maps in reactions with density of " + str(densities_remainders[density_counter]) + " %")
            print("----------------------------------------------------------------------")
            print("\n")

            # reinitialize variables
            reaction_counter = 0
            times_by_fixed_density = []

            # get each reaction of a given density
            for (in_G, in_H, in_anchor) in each_list_of_reactions_with_fixed_density:

                # update reaction counter
                reaction_counter = reaction_counter + 1

                # task message
                print("- reaction number: " + str(reaction_counter)
                      + " - density: " + str(densities_remainders[density_counter])
                      + " - order: " + str(each_order))

                # analysis and time with gm_extender
                print("|| running gm_extender" + tail)
                initial_time = time.time()
                all_extensions, good_map = gm.search_stable_extension(in_G, in_H,
                                                                      input_anchor = in_anchor,
                                                                      node_labels = True, edge_labels = True,
                                                                      all_extensions = False)
                final_time = time.time()
                out_time_gm_extender = final_time - initial_time

                # analysis and time gm_isomorphism
                print("|| running gm_isomorphism" + tail)
                initial_time = time.time()
                all_extensions, good_map = inteligent_labeling_isomorphism(in_G, in_H, in_anchor, "gm_iso")
                final_time = time.time()
                out_time_gm_isomorphism = final_time - initial_time

                # analysis and time nx_isomorphism
                print("|| running nx_isomorphism" + tail)
                initial_time = time.time()
                all_extensions, good_map = inteligent_labeling_isomorphism(in_G, in_H, in_anchor, "nx_iso")
                final_time = time.time()
                out_time_nx_isomorphism = final_time - initial_time

                # analysis and time nx_ismags
                print("|| running nx_ismags" + tail)
                initial_time = time.time()
                all_extensions, good_map = inteligent_labeling_isomorphism(in_G, in_H, in_anchor, "ismags")
                final_time = time.time()
                out_time_nx_ismags_algo = final_time - initial_time

                # save times of reaction
                times_by_fixed_density.append((out_time_gm_extender, out_time_gm_isomorphism, out_time_nx_isomorphism, out_time_nx_ismags_algo))
                print("\n")

            # update density counter
            density_counter = density_counter + 1

            # save times of all reactions with same density
            times_by_order.append(times_by_fixed_density)

    # dump times for graphs with each order
    file_name = "Times_Random_Graphs_" + str(each_order) + "_nodes.pkl"
    with open(file_name, "wb") as handle:
        pickle.dump(times_by_order, handle)
    print("######################################################################")
    print("\n")


# line jump message
print("\n")
print("> Finished")
print("\n")


################################################################################
################################################################################
