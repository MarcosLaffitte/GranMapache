################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: partial_maps                                                       #
#                                                                              #
# - Description: analysis of properties of partial maps, like maximum          #
#   connected extensions, overlaps, consistency, and others.                   #
#                                                                              #
################################################################################



# dependencies #################################################################



# already in python ------------------------------------------------------------
from copy import deepcopy



# not in python ----------------------------------------------------------------
import cython
import networkx as nx



# numpy in cython --------------------------------------------------------------
cimport numpy as cnp



# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match



# algorithms ###################################################################



# functions - maximum connected extension - wrapper ############################



# function: callable wrapper for the maximum connected extension ---------------
def maximum_connected_extension(G = nx.Graph(),   # can still receive DiGraph
                                H = nx.Graph(),   # can still receive DiGraph
                                anchor = []):
    # description
    """
    > description: receives two graphs G and H, and a match between them (here called
    anchor), and uses a VF2-like approach to obtain a maximum extension of the anchor
    producing a connected common subgraph (not necessarily maximum itslef). The anchor
    alone also produces a subgraph, which may not be an induced common subgraph, but
    the subgraph produced by the extension after removing the achor is always induced.

    > input:
    * G - first networkx (di)graph being matched.
    * H - second networkx (di)graph being matched.
    * anchor - inyective map as list of 2-tuples (x, y) of nodes x from G and y from H.

    > output:
    * extension - injective map as a list of 2-tuples (x, y) of nodes x from G and y
    from H representing the maximum connected extension of the anchor (contains the
    anchor as a sublist).
    * good_extension - boolean value indicating if the extension covers all nodes of
    G and H, i.e., if it is a bijection between G and H. If so, the anchor is what we
    have refered to as a "good partial atom map", and equivalenteĺy the match obtained
    by removing the anchhor from the extension is a graph-isomorphism between the graphs
    it induces from G and H.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.encode_match
    * .integerization.decode_match
    *
    """
    # output holders
    extension = []
    good_extension = False
    # cython variables
    # local variables
    encoded_graphs = []
    encoded_node_names = dict()
    encoded_node_label = dict()
    encoded_edge_label = dict()
    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([G, H])
    # encode match
    encoded_anchor = encode_match(anchor, encoded_node_names)
    # get cython-numpy structures for analysis
    # - nodes
    # - edges
    # - node labels
    # - edge labels
    # get total order for VF2-like analysis
    # get maximum extension
    # decode maximum extension
    # end of function
    return(anchor)



# functions - maximum connected extension - undirected #########################



# # function: recursive MATCH for undir MCS search -------------------------------
# def undirRecursiveExpansionMCS(someG1, someG2,
#                                existenceG1, existenceG2,
#                                containedG1, containedG2,
#                                someMatch = [], allMatches = [],
#                                ambiguous1 = [], ambiguous2 = [],
#                                score = "order",
#                                vLabels = True,
#                                eLabels = True,
#                                printProgressMCS = True,
#                                ambiguousPairsCheck = False,
#                                totOrder = dict()):
#     # local variables
#     progress = 0
#     expOrder = 0
#     scoreNewMatch = 0
#     scoreOldMatch = 0
#     newMatch = []
#     currMatch1 = []
#     currMatch2 = []
#     candidatePairs = []
#     forMatch = dict()
#     invMatch = dict()
#     allMatchesSet = set()
#     ansSiFy = False
#     ansSeFy = False
#     foundSubIso = False
#     progressReport = False
#     foundMaxColoredMatch = False
#     # get expected order for decision making
#     expOrder = min([someG1.order(), someG2.order()])
#     # define total order if not yet defined
#     if(len(totOrder) == 0):
#         totOrder = {v: i for (i, v) in list(enumerate(list(someG2.nodes()), start = 1))}
#         progressReport = True
#     # test initial alignment and improvement in alignment
#     if(len(allMatches) == 0):
#         if(len(someMatch) > 0):
#             # save match
#             if(len(someMatch) < expOrder):
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                         foundMaxColoredMatch = True
#                         return(allMatches, foundSubIso, foundMaxColoredMatch)
#             else:
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(ambiguousPairsCheck):
#                         foundSubIso = isSubIso(someG1, someG2, someMatch)
#                         if(foundSubIso):
#                             foundMaxColoredMatch = True
#                     else:
#                         foundSubIso = True
#                         foundMaxColoredMatch = True
#     else:
#         # pick score based on arguments
#         scoreNewMatch = matchScore(someG1, someG2, someMatch, existenceG1, existenceG2, containedG1, containedG2, score = score)
#         scoreOldMatch = matchScore(someG1, someG2, allMatches[0], existenceG1, existenceG2, containedG1, containedG2, score = score)
#         # save to MCS list if gets the same score of alignment
#         if(scoreNewMatch == scoreOldMatch):
#             allMatchesSet = [set(eachMatch) for eachMatch in allMatches]
#             if(not set(someMatch) in allMatchesSet):
#                 # append match
#                 if(len(someMatch) < expOrder):
#                     allMatches = allMatches + [someMatch]
#                     if(vLabels and eLabels):
#                         if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                             foundMaxColoredMatch = True
#                             return(allMatches, foundSubIso, foundMaxColoredMatch)
#                 else:
#                     allMatches = allMatches + [someMatch]
#                     if(vLabels and eLabels):
#                         if(ambiguousPairsCheck):
#                             foundSubIso = isSubIso(someG1, someG2, someMatch)
#                             if(foundSubIso):
#                                 foundMaxColoredMatch = True
#                         else:
#                             foundSubIso = True
#                             foundMaxColoredMatch = True
#         # overwrite MCS list if there is inprovement in alignment
#         if(scoreNewMatch > scoreOldMatch):
#             if(len(someMatch) < expOrder):
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                         foundMaxColoredMatch = True
#                         return(allMatches, foundSubIso, foundMaxColoredMatch)
#             else:
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(ambiguousPairsCheck):
#                         foundSubIso = isSubIso(someG1, someG2, someMatch)
#                         if(foundSubIso):
#                             foundMaxColoredMatch = True
#                     else:
#                         foundSubIso = True
#                         foundMaxColoredMatch = True
#     # pre-evaluate available pairs
#     if(len(someMatch) < expOrder):
#         # generate auxiliary structures
#         currMatch1 = [x for (x, y) in someMatch]
#         currMatch2 = [y for (x, y) in someMatch]
#         forMatch = {x:y for (x, y) in someMatch}
#         invMatch = {y:x for (x, y) in someMatch}
#         # get candidate pairs (if any)
#         candidatePairs = undirCandidatesMCS(someMatch, currMatch1, currMatch2, someG1, someG2, totOrder)
#         # evaluate candidate pairs
#         for (n1, n2) in candidatePairs:
#             # print progress only in first call
#             if(progressReport and printProgressMCS):
#                 progress = progress + 1
#                 printProgress(round(progress*100/len(candidatePairs), 2), progressIn = "in case", reportCase = False)
#             # evaluate sintactic feasibility
#             ansSiFy = undirSintacticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, invMatch, someG1, someG2,
#                                                    ambiguousPairsCheck, ambiguous1, ambiguous2)
#             if(ansSiFy):
#                 # evaluate semantic feasibility
#                 ansSeFy = undirSemanticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, someG1, someG2,
#                                                       vLabels = vLabels, eLabels = eLabels)
#                 if(ansSeFy):
#                     # DFS over feasible pairs
#                     newMatch = someMatch + [(n1, n2)]
#                     allMatches, foundSubIso, foundMaxColoredMatch = undirRecursiveExpansionMCS(someG1, someG2,
#                                                                                                existenceG1, existenceG2,
#                                                                                                containedG1, containedG2,
#                                                                                                newMatch, allMatches,
#                                                                                                ambiguous1 = ambiguous1, ambiguous2 = ambiguous2,
#                                                                                                score = score,
#                                                                                                vLabels = vLabels,
#                                                                                                eLabels = eLabels,
#                                                                                                ambiguousPairsCheck = ambiguousPairsCheck,
#                                                                                                totOrder = totOrder)
#                     # stop serach if one (sub)graph isomorphism was found (when preseving all labels)
#                     if(vLabels and eLabels and foundSubIso):
#                         break
#                     # stop search if one maximal colored match was found (when preserving all labels)
#                     if(vLabels and eLabels and foundMaxColoredMatch):
#                         break
#     # end of function
#     return(allMatches, foundSubIso, foundMaxColoredMatch)


# # function: get candidate pairs for undir MCS search ---------------------------
# def undirCandidatesMCS(someMatch, currMatch1, currMatch2, someG1, someG2, totOrder):
#     # local variables
#     maxMatchedIndex = 0
#     P = []
#     valid1 = []
#     valid2 = []
#     # get candidates preserving order (if no previous match just take everything)
#     if(len(someMatch) > 0):
#         maxMatchedIndex = max([totOrder[n2] for (n1, n2) in someMatch])
#     # get candidate pairs
#     valid1 = [x for x in list(someG1.nodes()) if(not x in currMatch1)]
#     valid2 = [y for y in list(someG2.nodes()) if((not y in currMatch2) and (totOrder[y] > maxMatchedIndex))]
#     P = list(product(valid1, valid2))
#     # end of function
#     return(P)



# # function: evaluate the sintactic feasability of mapping n1 to n2 -------------
# def undirSintacticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, invMatch, someG1, someG2,
#                                  ambiguousCheck, ambiguousG1, ambiguousG2):
#     # local variables
#     neigh1 = []
#     neigh2 = []
#     ambNeigh1 = []
#     ambNeigh2 = []
#     matchNeigh1 = []
#     matchNeigh2 = []
#     # get neighbors of n1 and n2
#     neigh1 = list(someG1.neighbors(n1))
#     neigh2 = list(someG2.neighbors(n2))
#     # loop-consistency-test
#     if((n1 in neigh1) and (not n2 in neigh2)):
#         return(False)
#     if((not n1 in neigh1) and (n2 in neigh2)):
#         return(False)
#     # look ahead 0: consistency of neighbors in match
#     matchNeigh1 = [x for x in neigh1 if(x in currMatch1)]
#     matchNeigh2 = [y for y in neigh2 if(y in currMatch2)]
#     # get ambiguous neighbors if requested
#     if(ambiguousCheck):
#         ambNeigh1 = list(set([u for (u, v) in ambiguousG1 if(v == n1)] + [v for (u, v) in ambiguousG1 if(u == n1)]))
#         ambNeigh2 = list(set([u for (u, v) in ambiguousG2 if(v == n2)] + [v for (u, v) in ambiguousG2 if(u == n2)]))
#     # compare neighborhoods
#     for v1 in matchNeigh1:
#         # check if either true or ambiguous neighbor
#         if(not forMatch[v1] in (matchNeigh2 + ambNeigh2)):
#             return(False)
#     for v2 in matchNeigh2:
#         # check if either true or ambiguous neighbor
#         if(not invMatch[v2] in (matchNeigh1 + ambNeigh1)):
#             return(False)
#     # end of function
#     return(True)



# # function: evaluate the semantic feasability of mapping n1 to n2 --------------
# def undirSemanticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, someG1, someG2,
#                                 vLabels = True, eLabels = True):
#     # local variables
#     neigh1 = []
#     neigh2 = []
#     matchNeigh1 = []
#     matchNeigh2 = []
#     # compare vertex-labels
#     if(vLabels):
#         if(not someG1.nodes[n1]["nodeLabel"] == someG2.nodes[n2]["nodeLabel"]):
#             return(False)
#     # compare edge labels
#     if(eLabels):
#         # get neighborhoods
#         neigh1 = list(someG1.neighbors(n1))
#         # compare loop-label (if any)
#         if(n1 in neigh1):
#             if(not someG1[n1][n1]["edgeLabel"] == someG2[n2][n2]["edgeLabel"]):
#                 return(False)
#         # compare edge-labels of true-edges and ignore ambiguous neighbors
#         neigh2 = list(someG2.neighbors(n2))
#         matchNeigh1 = [v for v in neigh1 if(v in currMatch1)]
#         matchNeigh2 = [v for v in neigh2 if(v in currMatch2)]
#         for v in matchNeigh1:
#             # only true or ambiguous neighbors at this point (this is the intersection)
#             if(forMatch[v] in matchNeigh2):
#                 if(not someG1[n1][v]["edgeLabel"] == someG2[n2][forMatch[v]]["edgeLabel"]):
#                     return(False)
#     # end of function
#     return(True)



# # functions - maximum connected extension - directed ###########################



# # function: recursive MATCH for dir MCS search ---------------------------------
# def dirRecursiveExpansionMCS(someG1, someG2,
#                              existenceG1, existenceG2,
#                              containedG1, containedG2,
#                              someMatch = [], allMatches = [],
#                              ambiguous1 = [], ambiguous2 = [],
#                              score = "order",
#                              vLabels = True,
#                              eLabels = True,
#                              printProgressMCS = True,
#                              ambiguousPairsCheck = False,
#                              totOrder = dict()):
#     # local variables
#     expOrder = 0
#     progress = 0
#     scoreNewMatch = 0
#     scoreOldMatch = 0
#     newMatch = []
#     currMatch1 = []
#     currMatch2 = []
#     candidatePairs = []
#     forMatch = dict()
#     invMatch = dict()
#     allMatchesSet = set()
#     ansSiFy = False
#     ansSeFy = False
#     foundSubIso = False
#     progressReport = False
#     foundMaxColoredMatch = False
#     # get expected order for decision making
#     expOrder = min([someG1.order(), someG2.order()])
#     # define total order if not yet defined
#     if(len(totOrder) == 0):
#         totOrder = {v: i for (i, v) in list(enumerate(list(someG2.nodes()), start = 1))}
#         progressReport = True
#     # test initial alignment and improvement in alignment
#     if(len(allMatches) == 0):
#         if(len(someMatch) > 0):
#             # save match
#             if(len(someMatch) < expOrder):
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                         foundMaxColoredMatch = True
#                         return(allMatches, foundSubIso, foundMaxColoredMatch)
#             else:
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(ambiguousPairsCheck):
#                         foundSubIso = isSubIso(someG1, someG2, someMatch)
#                         if(foundSubIso):
#                             foundMaxColoredMatch = True
#                     else:
#                         foundSubIso = True
#                         foundMaxColoredMatch = True
#     else:
#         # pick score based on arguments
#         scoreNewMatch = matchScore(someG1, someG2, someMatch, existenceG1, existenceG2, containedG1, containedG2, score = score)
#         scoreOldMatch = matchScore(someG1, someG2, allMatches[0], existenceG1, existenceG2, containedG1, containedG2, score = score)
#         # save to MCS list if gets the same score of alignment
#         if(scoreNewMatch == scoreOldMatch):
#             allMatchesSet = [set(eachMatch) for eachMatch in allMatches]
#             if(not set(someMatch) in allMatchesSet):
#                 # append match
#                 if(len(someMatch) < expOrder):
#                     allMatches = allMatches + [someMatch]
#                     if(vLabels and eLabels):
#                         if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                             foundMaxColoredMatch = True
#                             return(allMatches, foundSubIso, foundMaxColoredMatch)
#                 else:
#                     allMatches = allMatches + [someMatch]
#                     if(vLabels and eLabels):
#                         if(ambiguousPairsCheck):
#                             foundSubIso = isSubIso(someG1, someG2, someMatch)
#                             if(foundSubIso):
#                                 foundMaxColoredMatch = True
#                         else:
#                             foundSubIso = True
#                             foundMaxColoredMatch = True
#         # overwrite MCS list if there is inprovement in alignment
#         if(scoreNewMatch > scoreOldMatch):
#             if(len(someMatch) < expOrder):
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(isMaxColoredMatch(someG1, someG2, someMatch)):
#                         foundMaxColoredMatch = True
#                         return(allMatches, foundSubIso, foundMaxColoredMatch)
#             else:
#                 allMatches = [someMatch]
#                 if(vLabels and eLabels):
#                     if(ambiguousPairsCheck):
#                         foundSubIso = isSubIso(someG1, someG2, someMatch)
#                         if(foundSubIso):
#                             foundMaxColoredMatch = True
#                     else:
#                         foundSubIso = True
#                         foundMaxColoredMatch = True
#     # pre-evaluate available pairs
#     if(len(someMatch) < expOrder):
#         # generate auxiliary structures
#         currMatch1 = [x for (x, y) in someMatch]
#         currMatch2 = [y for (x, y) in someMatch]
#         forMatch = {x:y for (x, y) in someMatch}
#         invMatch = {y:x for (x, y) in someMatch}
#         # get candidate pairs (if any)
#         candidatePairs = dirCandidatesMCS(someMatch, currMatch1, currMatch2, someG1, someG2, totOrder)
#         # evaluate candidate pairs
#         for (n1, n2) in candidatePairs:
#             # print progress only in first call
#             if(progressReport and printProgressMCS):
#                 progress = progress + 1
#                 printProgress(round(progress*100/len(candidatePairs), 2), progressIn = "in case", reportCase = False)
#             # evaluate sintactic feasibility
#             ansSiFy = dirSintacticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, invMatch, someG1, someG2,
#                                                  ambiguousPairsCheck, ambiguous1, ambiguous2)
#             if(ansSiFy):
#                 # evaluate semantic feasibility
#                 ansSeFy = dirSemanticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, someG1, someG2,
#                                                     vLabels = vLabels, eLabels = eLabels)
#                 if(ansSeFy):
#                     # DFS over feasible pairs
#                     newMatch = someMatch + [(n1, n2)]
#                     allMatches, foundSubIso, foundMaxColoredMatch = dirRecursiveExpansionMCS(someG1, someG2,
#                                                                                              existenceG1, existenceG2,
#                                                                                              containedG1, containedG2,
#                                                                                              newMatch, allMatches,
#                                                                                              ambiguous1 = ambiguous1, ambiguous2 = ambiguous2,
#                                                                                              score = score,
#                                                                                              vLabels = vLabels,
#                                                                                              eLabels = eLabels,
#                                                                                              ambiguousPairsCheck = ambiguousPairsCheck,
#                                                                                              totOrder = totOrder)
#                     # stop serach is one (sub)graph isomorphism was found (when preseving all labels)
#                     if(vLabels and eLabels and foundSubIso):
#                         break
#                     # stop search if one maximal colored match was found (when preseving all labels)
#                     if(vLabels and eLabels and foundMaxColoredMatch):
#                         break
#     # end of function
#     return(allMatches, foundSubIso, foundMaxColoredMatch)



# # function: get candidate pairs for dir MCS search ------------------------------
# def dirCandidatesMCS(someMatch, currMatch1, currMatch2, someG1, someG2, totOrder):
#     # local variables
#     maxMatchedIndex = 0
#     P = []
#     valid1 = []
#     valid2 = []
#     # get candidates preserving order (if no previous match just take everything)
#     if(len(someMatch) > 0):
#         maxMatchedIndex = max([totOrder[n2] for (n1, n2) in someMatch])
#     # alternatively try pairing all unpaired vertices
#     valid1 = [x for x in list(someG1.nodes()) if(not x in currMatch1)]
#     valid2 = [y for y in list(someG2.nodes()) if((not y in currMatch2) and (totOrder[y] > maxMatchedIndex))]
#     P = list(product(valid1, valid2))
#     # end of function
#     return(P)



# # function: evaluate the sintactic feasability of mapping n1 to n2 -------------
# def dirSintacticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, invMatch, someG1, someG2,
#                                ambiguousCheck, ambiguousG1, ambiguousG2):
#     # local variables
#     inNeigh1 = []
#     inNeigh2 = []
#     outNeigh1 = []
#     outNeigh2 = []
#     ambNeigh1 = []
#     ambNeigh2 = []
#     inMatchNeigh1 = []
#     inMatchNeigh2 = []
#     outMatchNeigh1 = []
#     outMatchNeigh2 = []
#     # get neighbors of n1 and n2
#     inNeigh1 = list(someG1.predecessors(n1))
#     inNeigh2 = list(someG2.predecessors(n2))
#     outNeigh1 = list(someG1.neighbors(n1))
#     outNeigh2 = list(someG2.neighbors(n2))
#     # loop-consistency-test
#     if((n1 in inNeigh1) and (not n2 in inNeigh2)):
#         return(False)
#     if((not n1 in inNeigh1) and (n2 in inNeigh2)):
#         return(False)
#     # get ambiguous neighbors if requested
#     if(ambiguousCheck):
#         ambNeigh1 = list(set([u for (u, v) in ambiguousG1 if(v == n1)] + [v for (u, v) in ambiguousG1 if(u == n1)]))
#         ambNeigh2 = list(set([u for (u, v) in ambiguousG2 if(v == n2)] + [v for (u, v) in ambiguousG2 if(u == n2)]))
#     # look ahead 0: consistency of neighbors in match
#     inMatchNeigh1 = [a1 for a1 in inNeigh1 if(a1 in currMatch1)]
#     inMatchNeigh2 = [a2 for a2 in inNeigh2 if(a2 in currMatch2)]
#     for a1 in inMatchNeigh1:
#         if(not forMatch[a1] in (inMatchNeigh2 + ambNeigh2)):
#             return(False)
#     for a2 in inMatchNeigh2:
#         if(not invMatch[a2] in (inMatchNeigh1 + ambNeigh1)):
#             return(False)
#     outMatchNeigh1 = [b1 for b1 in outNeigh1 if(b1 in currMatch1)]
#     outMatchNeigh2 = [b2 for b2 in outNeigh2 if(b2 in currMatch2)]
#     for b1 in outMatchNeigh1:
#         if(not forMatch[b1] in (outMatchNeigh2 + ambNeigh2)):
#             return(False)
#     for b2 in outMatchNeigh2:
#         if(not invMatch[b2] in (outMatchNeigh1 + ambNeigh1)):
#             return(False)
#     # end of function
#     return(True)



# # function: evaluate the semantic feasability of mapping n1 to n2 --------------
# def dirSemanticFeasabilityMCS(n1, n2, currMatch1, currMatch2, forMatch, someG1, someG2,
#                               vLabels = True, eLabels = True):
#     # local variables
#     inNeigh1 = []
#     inNeigh2 = []
#     outNeigh1 = []
#     outNeigh2 = []
#     inMatchNeigh1 = []
#     inMatchNeigh2 = []
#     outMatchNeigh1 = []
#     outMatchNeigh2 = []
#     # compare vertex-labels
#     if(vLabels):
#         if(not someG1.nodes[n1]["nodeLabel"] == someG2.nodes[n2]["nodeLabel"]):
#             return(False)
#     # compare edge labels
#     if(eLabels):
#         # get neighborhoods
#         inNeigh1 = list(someG1.predecessors(n1))
#         outNeigh1 = list(someG1.neighbors(n1))
#         # compare loop-label (if any)
#         if(n1 in inNeigh1):
#             if(not someG1[n1][n1]["edgeLabel"] == someG2[n2][n2]["edgeLabel"]):
#                 return(False)
#         # compare edge-labels of true-edges and ignore ambiguous neighbors
#         inNeigh2 = list(someG2.predecessors(n2))
#         outNeigh2 = list(someG2.neighbors(n2))
#         inMatchNeigh1 = [a for a in inNeigh1 if(a in currMatch1)]
#         outMatchNeigh1 = [b for b in outNeigh1 if(b in currMatch1)]
#         inMatchNeigh2 = [a for a in inNeigh2 if(a in currMatch2)]
#         outMatchNeigh2 = [b for b in outNeigh2 if(b in currMatch2)]
#         for a in inMatchNeigh1:
#             # only true or ambiguous neighbors at this point (this is the in-intersection)
#             if(forMatch[a] in inMatchNeigh2):
#                 if(not someG1[a][n1]["edgeLabel"] == someG2[forMatch[a]][n2]["edgeLabel"]):
#                     return(False)
#         for b in outMatchNeigh1:
#             # only true or ambiguous neighbors at this point (this is the out-intersection)
#             if(forMatch[b] in outMatchNeigh2):
#                 if(not someG1[n1][b]["edgeLabel"] == someG2[n2][forMatch[b]]["edgeLabel"]):
#                     return(False)
#     # end of function
#     return(True)



################################################################################
################################################################################
