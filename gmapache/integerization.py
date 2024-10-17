############################################################
#                                                          #
# - GranMapache: GRAphs-and-Networks MAPping Applications  #
#   with Cython and HEuristics                             #
#                                                          #
# - Module: integerization                                 #
#                                                          #
############################################################


# dependencies #############################################


# already in python ----------------------------------------
import time
import random
from copy import deepcopy
from math import modf, sqrt
from operator import eq, itemgetter
from itertools import product, combinations
from sys import argv, exit, getrecursionlimit, setrecursionlimit


# not in python --------------------------------------------
import cython


# functions ################################################


# function: transform labeled graphs to integers -----------
def encoding(input_graphs = []):
    # local variables



# function: condense labels into one dict for comparison -----------------------
def condensedLabel(someG):
    # local variables
    uniformG = None
    # get graph of corresponding type
    if(someG.is_directed()):
        uniformG = nx.DiGraph()
    else:
        uniformG = nx.Graph()
    # iterate over nodes condensing their labels
    for (v, nodeInfo) in list(someG.nodes(data = True)):
        uniformG.add_node(v, nodeLabel = nodeInfo)
    # iterate over edges condensing their labels
    for (u, v, edgeInfo) in list(someG.edges(data = True)):
        uniformG.add_edge(u, v, edgeLabel = edgeInfo)
    # end of function
    return(uniformG)


# function: extend condensed labels inside dict into single labels -------------
def extendedLabel(someG):
    # local variables
    extendedG = None
    nodeAttributes = dict()
    edgeAttributes = dict()
    # get graph of corresponding type
    if(someG.is_directed()):
        extendedG = nx.DiGraph()
    else:
        extendedG = nx.Graph()
    # get attributes of nodes
    for (v, nodeInfo) in list(someG.nodes(data = True)):
        extendedG.add_node(v)
        nodeAttributes[v] = deepcopy(nodeInfo["nodeLabel"])
    # get attributes of edges
    for (u, v, edgeInfo) in list(someG.edges(data = True)):
        extendedG.add_edge(u, v)
        edgeAttributes[(u, v)] = deepcopy(edgeInfo["edgeLabel"])
    # asign new attributes
    nx.set_node_attributes(extendedG, nodeAttributes)
    nx.set_edge_attributes(extendedG, edgeAttributes)
    # end of function
    return(extendedG)


# function: condense labels into one single string -----------------------------
def condensedLabelToStr(someG):
    # local variables
    strLabels = []
    strLabelFinal = ""
    uniformG = None
    nodeInfoTuple = None
    edgeInfoTuple = None
    # get graph of corresponding type
    if(someG.is_directed()):
        uniformG = nx.DiGraph()
    else:
        uniformG = nx.Graph()
    # iterate over nodes condensing their labels into tuple
    for (v, nodeInfo) in list(someG.nodes(data = True)):
        nodeInfoTuple = [(labelName, nodeInfo[labelName]) for labelName in list(nodeInfo.keys())]
        nodeInfoTuple.sort()
        nodeInfoTuple = tuple(nodeInfoTuple)
        if(len(nodeInfoTuple) > 0):
            strLabels = []
            for (labelName, labelValue) in nodeInfoTuple:
                strLabels.append("(" + str(labelName) + ", " + str(labelValue) + ")")
            strLabelFinal = ", ".join(strLabels)
        else:
            strLabelFinal = "empty"
        uniformG.add_node(v, nodeLabelStr = strLabelFinal)
    # iterate over edges condensing their labels into tuple
    for (u, v, edgeInfo) in list(someG.edges(data = True)):
        edgeInfoTuple = [(labelName, edgeInfo[labelName]) for labelName in list(edgeInfo.keys())]
        edgeInfoTuple.sort()
        edgeInfoTuple = tuple(edgeInfoTuple)
        if(len(edgeInfoTuple) > 0):
            strLabels = []
            for (labelName, labelValue) in edgeInfoTuple:
                strLabels.append("(" + str(labelName) + ", " + str(labelValue) + ")")
            strLabelFinal = ", ".join(strLabels)
        else:
            strLabelFinal = "empty"
        uniformG.add_edge(u, v, edgeLabelStr = strLabelFinal)
    # end of function
    return(uniformG)


############################################################
############################################################
