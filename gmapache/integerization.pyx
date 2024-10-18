################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: integerization                                                     #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import time
import random
from copy import deepcopy
from math import modf, sqrt
from operator import eq, itemgetter
from itertools import product, combinations
from sys import argv, exit, getrecursionlimit, setrecursionlimit


# not in python ----------------------------------------------------------------
import cython


# functions ####################################################################


# function: homogenizes dicts of labels of list of graphs into integers --------
def encode_graphs(input_graphs = []):
    # description
    """
    - Receives a list of inputs graphs and turns each dictionary of
    node_labels and edge_labels into integers, considering their
    repetitions across the list.
    - Nodes are also renamed into integers starting from 1.
    - The function returns a list with the copies of the input_graph
    preserving their order, together with dictionaries from integers
    into the original node_labels, edge_labels, and node_names.
    """
    # cython variables
    cdef int i = 0
    # local variables
    new_node_label = 1
    new_edge_label = 1
    all_nodes = []
    all_node_labels = []
    all_edge_labels = []
    encoded_graphs = []
    node_name_encoding = dict()        # from ints to node names
    node_label_encoding = dict()       # from ints to node label-dicts
    edge_label_encoding = dict()       # from ints to edge label-dicts
    node_name_encoding_inv = dict()    # from node names to ints
    node_label_encoding_inv = [None]   # from node label-dicts to ints (indices)
    edge_label_encoding_inv = [None]   # from edge label-dicts to ints (indices)
    # homogenize node names across input graphs
    for eachGraph in input_graphs:
        all_nodes = list(set(all_nodes + list(eachGraph.nodes())))
    for i in range(len(all_nodes)):
        node_name_encoding[int(i+1)] = deepcopy(all_nodes[i])
        node_name_encoding_inv[deepcopy(all_nodes[i])] = int(i+1)
    # homogenize node labels across input graphs
    for eachGraph in input_graphs:
        for (v, nodeInfo) in list(eachGraph.nodes(data = True)):
            if(nodeInfo not in node_label_encoding_inv):
                node_label_encoding[new_node_label] = deepcopy(nodeInfo)
                node_label_encoding_inv.append(deepcopy(nodeInfo))
                new_node_label = new_node_label + 1                
    # homogenize edge labels across input graphs
    for eachGraph in input_graphs:
        for (u, v, edgeInfo) in list(eachGraph.edges(data = True)):
            if(edgeInfo not in edge_label_encoding_inv):
                edge_label_encoding[new_edge_label] = deepcopy(edgeInfo)
                edge_label_encoding_inv.append(deepcopy(edgeInfo))
                new_edge_label = new_edge_label + 1

    
    # end of function
    return([])







# # function: condense labels into one dict for comparison -----------------------
# def condensedLabel(someG):
#     # local variables
#     uniformG = None
#     # get graph of corresponding type
#     if(someG.is_directed()):
#         uniformG = nx.DiGraph()
#     else:
#         uniformG = nx.Graph()
#     # iterate over nodes condensing their labels
#     for (v, nodeInfo) in list(someG.nodes(data = True)):
#         uniformG.add_node(v, nodeLabel = nodeInfo)
#     # iterate over edges condensing their labels
#     for (u, v, edgeInfo) in list(someG.edges(data = True)):
#         uniformG.add_edge(u, v, edgeLabel = edgeInfo)
#     # end of function
#     return(uniformG)


# # function: extend condensed labels inside dict into single labels -------------
# def extendedLabel(someG):
#     # local variables
#     extendedG = None
#     nodeAttributes = dict()
#     edgeAttributes = dict()
#     # get graph of corresponding type
#     if(someG.is_directed()):
#         extendedG = nx.DiGraph()
#     else:
#         extendedG = nx.Graph()
#     # get attributes of nodes
#     for (v, nodeInfo) in list(someG.nodes(data = True)):
#         extendedG.add_node(v)
#         nodeAttributes[v] = deepcopy(nodeInfo["nodeLabel"])
#     # get attributes of edges
#     for (u, v, edgeInfo) in list(someG.edges(data = True)):
#         extendedG.add_edge(u, v)
#         edgeAttributes[(u, v)] = deepcopy(edgeInfo["edgeLabel"])
#     # asign new attributes
#     nx.set_node_attributes(extendedG, nodeAttributes)
#     nx.set_edge_attributes(extendedG, edgeAttributes)
#     # end of function
#     return(extendedG)


# # function: condense labels into one single string -----------------------------
# def condensedLabelToStr(someG):
#     # local variables
#     strLabels = []
#     strLabelFinal = ""
#     uniformG = None
#     nodeInfoTuple = None
#     edgeInfoTuple = None
#     # get graph of corresponding type
#     if(someG.is_directed()):
#         uniformG = nx.DiGraph()
#     else:
#         uniformG = nx.Graph()
#     # iterate over nodes condensing their labels into tuple
#     for (v, nodeInfo) in list(someG.nodes(data = True)):
#         nodeInfoTuple = [(labelName, nodeInfo[labelName]) for labelName in list(nodeInfo.keys())]
#         nodeInfoTuple.sort()
#         nodeInfoTuple = tuple(nodeInfoTuple)
#         if(len(nodeInfoTuple) > 0):
#             strLabels = []
#             for (labelName, labelValue) in nodeInfoTuple:
#                 strLabels.append("(" + str(labelName) + ", " + str(labelValue) + ")")
#             strLabelFinal = ", ".join(strLabels)
#         else:
#             strLabelFinal = "empty"
#         uniformG.add_node(v, nodeLabelStr = strLabelFinal)
#     # iterate over edges condensing their labels into tuple
#     for (u, v, edgeInfo) in list(someG.edges(data = True)):
#         edgeInfoTuple = [(labelName, edgeInfo[labelName]) for labelName in list(edgeInfo.keys())]
#         edgeInfoTuple.sort()
#         edgeInfoTuple = tuple(edgeInfoTuple)
#         if(len(edgeInfoTuple) > 0):
#             strLabels = []
#             for (labelName, labelValue) in edgeInfoTuple:
#                 strLabels.append("(" + str(labelName) + ", " + str(labelValue) + ")")
#             strLabelFinal = ", ".join(strLabels)
#         else:
#             strLabelFinal = "empty"
#         uniformG.add_edge(u, v, edgeLabelStr = strLabelFinal)
#     # end of function
#     return(uniformG)


################################################################################
################################################################################
