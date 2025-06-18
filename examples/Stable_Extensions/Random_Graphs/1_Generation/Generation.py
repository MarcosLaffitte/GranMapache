################################################################################
#                                                                              #
#  - Generate Random Graphs for Extension of Partial Atonm-to-Atom Maps        #
#                                                                              #
#  - Generation                                                                #
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


# custom dependencies ----------------------------------------------------------
import gmapache as gm


# parameters ###################################################################


# number of graphs -------------------------------------------------------------
reactions_per_density = 10


# parameters of remainder graphs -----------------------------------------------
nodes_remainders = [200, 175, 150, 125, 100]
# NOTE: density of at least 2 is needed for connectedness with 100 nodes,
# but with density exactly 2 and 125 nodes it seems way too hard to find a random connected graph,
# and density of at most 97 is also needed for possibility of construction with 100 vertices,
# thus the 3-97 range for density was selected, which seems to work well
densities_remainders = list(range(3, 97 + 1))
density_breakpoints = [10, 20, 30, 40, 50, 60, 70, 80, 90, 97]


# parameters of reaction centers -----------------------------------------------
nodes_reaction_center = 15
edges_reaction_center = 20
anchor = [(k, k) for k in range(1, nodes_reaction_center + 1)]


# label sets -------------------------------------------------------------------
node_label_set = ["a1", "a2", "a3", "a4", "a5", "a6"]
edge_label_set = ["b1", "b2", "b3"]
edge_origin_set = ["G", "H", "GH"]


# data holders -----------------------------------------------------------------
# all underlying graphs
done_underlying_ITS_graphs = []
# produce connected underlyging graph of remainder with false reaction center
underlying_remainder = None
# produce connected underlyging graph of true reacion center
underlying_center = None
# remove false reaction edges from underlyging graph of remainder
reduced_remainder = None
false_center = None
# add true reaction center to remainder and cover edge-defficit to match initial density
underlying_ITS = None
remainder_edges = []
center_edges = []
# obtain labeled graphs
original_G = None
original_H = None
# aleatorize node names and positions
randomized_G = None
randomized_H = None
randomized_anchor = []


# output -----------------------------------------------------------------------
out_G = None
out_H = None
out_anchor = None
reactions_by_order = []
reactions_with_fixed_density = []


# functions ####################################################################


# buid unlabeled remainder graph -----------------------------------------------
def generate_unlabeled_connected_remainder(remainder_order, remainder_density):

    # local variables
    remainder = None
    connected = False
    remainder_size = 0
    rename_nodes = dict()

    # get size from ceiling of density
    remainder_size = math.ceil((remainder_density/100)*((remainder_order*(remainder_order-1))/2))

    # generate random graph
    while(not connected):
        # generate arbitrary random graph or tree if necessary
        if(remainder_size == remainder_order-1):
            remainder = nx.random_unlabeled_tree(remainder_order)
        else:
            remainder = nx.gnm_random_graph(remainder_order, remainder_size, directed = False)
        # check connectednes in general case
        connected = nx.is_connected(remainder)

    # rename nodes
    rename_nodes = {k:k+1 for k in list(remainder.nodes())}
    remainder = nx.relabel_nodes(remainder, rename_nodes, copy = True)

    # end of function
    return(remainder)


# buid unlabeled connected reaction center -------------------------------------
def generate_unlabeled_connected_reaction_center(rc_order, rc_size):

    # local variables
    center = None
    connected = False
    rename_nodes = dict()

    # generate true reaction center
    while(not connected):
        center = nx.gnm_random_graph(rc_order, rc_size, directed = False)
        connected = nx.is_connected(center)

    # rename nodes
    rename_nodes = {k:k+1 for k in list(center.nodes())}
    center = nx.relabel_nodes(center, rename_nodes, copy = True)

    # end of function
    return(center)


# remove false reaction edges --------------------------------------------------
def remove_false_reaction_center(remainder_graph):

    # global variables
    global nodes_reaction_center

    # local variables
    false_center = None
    reduced_remainder = None
    false_center_size = 0
    false_center_edges = []
    remainder_graph_nodes = list(range(nodes_reaction_center + 1, remainder_graph.order() + 1))

    # induce false center
    false_center = deepcopy(remainder_graph)
    false_center.remove_nodes_from(remainder_graph_nodes)
    false_center_edges = list(false_center.edges())
    false_center_size = false_center.size()

    # get reduced remainder
    reduced_remainder = deepcopy(remainder_graph)
    reduced_remainder.remove_edges_from(false_center_edges)

    # end of function
    return(reduced_remainder, false_center_size)


# add true reaction edges and recover density ----------------------------------
def recover_true_reaction_center_and_density(reduced_remainder_graph, reaction_center_graph, false_center_size):

    # local variables
    u = 0
    v = 0
    center_edges = []
    remainder_edges = []
    underlying_ITS = None
    valid = False

    # add true reaction center
    underlying_ITS = deepcopy(reduced_remainder_graph)
    underlying_ITS.add_edges_from(list(reaction_center_graph.edges()))

    # recover density of remainder
    for i in range(false_center_size):
        # get new edge at random
        valid = False
        while(not valid):
            # sample edge
            candidate = list(random.sample(list((nx.complement(underlying_ITS)).edges()), 1))[0]
            u = candidate[0]
            v = candidate[1]
            # test validity of candidate
            if(not ((u in list(reaction_center_graph.nodes())) and (v in list(reaction_center_graph.nodes())))):
                remainder_edges.append((u, v))
                underlying_ITS.add_edge(u, v)
                valid = True

    # save edges
    center_edges = list(reaction_center_graph.edges())
    remainder_edges = remainder_edges + list(reduced_remainder_graph.edges())

    # end of function
    return(underlying_ITS, remainder_edges, center_edges)


# get labeled reactans and products sides --------------------------------------
def get_labeled_reaction(underlying_ITS, remainder_edges, center_edges):

    # global variables
    global node_label_set # = ["a1", "a2", "a3", "a4", "a5", "a6"]
    global edge_label_set # = ["b1", "b2", "b3"]
    global edge_origin_set # = ["G", "H", "GH"]

    # local variables
    G = None
    H = None
    u = 0
    v = 0
    atom_label = ""
    bond_label = ""
    edge_origin = ""
    bond_tuple = ()

    # add labeled vertices
    G = nx.Graph()
    H = nx.Graph()
    for v in list(underlying_ITS.nodes()):
        atom_label = list(random.sample(node_label_set, 1))[0]
        G.add_node(v, atom_type = atom_label)
        H.add_node(v, atom_type = atom_label)

    # add remainder labeld edges
    for (u, v) in remainder_edges:
        bond_label = list(random.sample(edge_label_set, 1))[0]
        G.add_edge(u, v, bond_type = bond_label)
        H.add_edge(u, v, bond_type = bond_label)

    # add reaction edges
    for (u, v) in center_edges:
        # choose to which graph belongs the edge
        edge_origin = list(random.sample(edge_origin_set, 1))[0]
        # add reaction edge only to G
        if(edge_origin == "G"):
            bond_label = list(random.sample(edge_label_set, 1))[0]
            G.add_edge(u, v, bond_type = bond_label)
        # add reaction edge only to H
        if(edge_origin == "H"):
            bond_label = list(random.sample(edge_label_set, 1))[0]
            H.add_edge(u, v, bond_type = bond_label)
        # add reaction edge to both graphs
        if(edge_origin == "GH"):
            bond_tuple = list(random.sample(edge_label_set, 2))
            G.add_edge(u, v, bond_type = bond_tuple[0])
            H.add_edge(u, v, bond_type = bond_tuple[1])

    # end of function
    return(G, H)


# do anchor preserving randomization -------------------------------------------
def make_anchor_preserving_randomization(G, H):

    # global variables
    global anchor

    # local variables
    new_name = 0
    rand_G = None
    rand_H = None
    renamed_G = None
    renamed_H = None
    rand_anchor = []
    old_names_G = []
    old_names_H = []
    nodes_G_renamed = []
    nodes_H_renamed = []
    edges_G_renamed = []
    edges_H_renamed = []
    new_names_G = dict()
    new_names_H = dict()

    # copy input graphs
    renamed_G = deepcopy(G)
    renamed_H = deepcopy(H)

    # rename vertices
    old_names_G = list(renamed_G.nodes())
    old_names_H = list(renamed_H.nodes())
    for each_node in list(renamed_G.nodes()):
        new_name = random.choice(old_names_G)
        new_names_G[each_node] = "x_" + str(new_name)
        old_names_G.remove(new_name)
    for each_node in list(renamed_H.nodes()):
        new_name = random.choice(old_names_H)
        new_names_H[each_node] = "y_" + str(new_name)
        old_names_H.remove(new_name)
    nx.relabel_nodes(renamed_G, new_names_G, copy = False)
    nx.relabel_nodes(renamed_H, new_names_H, copy = False)

    # relabel anchor
    rand_anchor = []
    for (a, b) in anchor:
        rand_anchor.append((new_names_G[a], new_names_H[b]))

    # aleatorized default order of nodes and edges
    nodes_G_renamed = list(renamed_G.nodes(data = True))
    nodes_H_renamed = list(renamed_H.nodes(data = True))
    edges_G_renamed = list(renamed_G.edges(data = True))
    edges_H_renamed = list(renamed_H.edges(data = True))
    random.shuffle(nodes_G_renamed)
    random.shuffle(nodes_H_renamed)
    random.shuffle(edges_G_renamed)
    random.shuffle(edges_H_renamed)

    # make aleatorized copies
    rand_G = nx.Graph()
    rand_G.add_nodes_from(nodes_G_renamed)
    rand_G.add_edges_from(edges_G_renamed)
    rand_H = nx.Graph()
    rand_H.add_nodes_from(nodes_H_renamed)
    rand_H.add_edges_from(edges_H_renamed)

    # end of function
    return(rand_G, rand_H, rand_anchor)


# analysis #####################################################################


# task message
print("\n")
print("> Generating random graphs")
print("\n")
print("----------------------------------------------------------------------")


# fix number of vertices of remainder graphs
for each_order in nodes_remainders:

    # reinitialize variables
    reactions_by_order = []

    # fix density of remainder graphs
    for each_density in densities_remainders:

        # reinitialize variables
        done_underlying_ITS_graphs = []
        reactions_with_fixed_density = []

        # task message
        print("* generating random ITS graphs with - order: ", each_order, " - remainder density: ", each_density)

        # start time for reference
        initial_time = time.time()

        # generate ITS graphs whose remainder has fixed density
        while(len(done_underlying_ITS_graphs) < reactions_per_density):

            # produce connected underlyging graph of remainder with false reaction center
            underlying_remainder = generate_unlabeled_connected_remainder(each_order, each_density)

            # produce connected underlyging graph of true reacion center
            underlying_center = generate_unlabeled_connected_reaction_center(nodes_reaction_center, edges_reaction_center)

            # remove false reaction edges from underlyging graph of remainder
            reduced_remainder, false_center = remove_false_reaction_center(underlying_remainder)

            # add true reaction center to remainder and cover edge-defficit to match initial density
            underlying_ITS, remainder_edges, center_edges = recover_true_reaction_center_and_density(reduced_remainder, underlying_center, false_center)

            # test connectendness of ITS just in case (all should be conected by construction, this is just a safe-check)
            connected = nx.is_connected(underlying_ITS)
            if(not connected):
                print("- NOT COOL")

            # continue if connected
            if(connected):

                # check underlying ITS graph is not isomorphic to the saved ones
                isomorphic = False
                for each_graph in done_underlying_ITS_graphs:

                    # here networkx VF2-isomorphism uses the construction-order of the reaction center,
                    # which should not be done in experiments with randomize data for the sake of a fair comparison
                    isomorphic = nx.is_isomorphic(underlying_ITS, each_graph)

                    # discard if isomorphic
                    if(isomorphic):
                        break

                # continue if underlying is new
                if(not isomorphic):

                    # if underlying is not isomorphic then save
                    done_underlying_ITS_graphs.append(underlying_ITS)

                    # obtain labeled graphs
                    original_G, original_H  = get_labeled_reaction(underlying_ITS, remainder_edges, center_edges)

                    # aleatorize node names and positions to prevent NetworkX VF2-isomorphism from using the construction-order for G and H,
                    # which should not be done in experiments with randomize data for the sake of a fair comparison
                    randomized_G, randomized_H, randomized_anchor = make_anchor_preserving_randomization(original_G, original_H)

                    # store new reaction
                    out_G = deepcopy(randomized_G)
                    out_H = deepcopy(randomized_H)
                    out_anchor = deepcopy(randomized_anchor)
                    reactions_with_fixed_density.append((out_G, out_H, out_anchor))

        # reinitialize variables
        reactions_by_order.append(reactions_with_fixed_density)

        # print time taken for generation
        final_time = time.time()
        print("|| total time [s]: ", final_time - initial_time)

        # dump graphs with fixed order in pickle file
        if(each_density in density_breakpoints):
            # make pickle file
            file_name = "Random_Graphs_" + str(each_order) + "_nodes_and_upper_density_" + str(each_density) + ".pkl"
            print("|| Saving file: " + file_name)
            with open(file_name, "wb") as handle:
                pickle.dump(reactions_by_order, handle)
            # reinitialize variables
            reactions_by_order = []
            print("----------------------------------------------------------------------")


# line jump message
print("\n")
print("> Finished")
print("\n")


################################################################################
################################################################################
