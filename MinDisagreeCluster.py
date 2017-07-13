#!/usr/bin/env python


import argparse
ap = argparse.ArgumentParser(description="Heuristic for min-disagree clustering")
ap.add_argument("--graph", help="Graph file.", required=True)
ap.add_argument("--radius", help="Radius for repulsion", type=int, default=None)
ap.add_argument("--penalty", help="Repulsion edge penalty.", type=float, default=1.0)
ap.add_argument("--factor", help="Number of repulstion nodes to add relative to neighboring", type=float,default=1)
ap.add_argument("--swap", help="Perform this many iterations of swapping nodes at cut boundaries", type=int, default=0)
ap.add_argument("--layout", help="Apply layout from this graph.", default=None)
ap.add_argument("--out", help="Output graph.")
ap.add_argument("--plot", help="Plot iterations", default=None)
ap.add_argument("--plotRepulsion", help="Plot repulsion", default=False,action='store_true')
ap.add_argument("--niter", help="Number of clustering iterations.", type=int,default=10)
ap.add_argument("--embed", help="Stop in ipython shell.", action='store_true', default=False)
ap.add_argument("--cuts", help="Write cuts to this file.", default=None)

ap.add_argument("--sites", help="Write sites of cuts to this file.", default=None)
args = ap.parse_args()


import networkx as nx
import ABPUtils
import itertools
import numpy as np
import pdb
from sets import Set
import copy
import random
from multiprocessing import Pool
import operator



g = ABPUtils.ReadGraph(args.graph)


sites = [ (int(g.node[n]['pos']), n) for n in g.node ]
sites = sorted(sites)

i = 0
nRepulsion = 0
nOK = 0


def StoreRepulsionAdjacency(g, n, pos, radius=None):
    repulsion = Set()
    if radius is None:
        neighbors = g[n].keys()        
        start = min([pos[node] for node in neighbors])
        end   = max([pos[node] for node in neighbors])
        
    for i in g.nodes():
        addRepulsion = False
        if radius is None:
            if pos[i] > start and  pos[i] < end and i not in g[n] and i != n:
                addRepulsion = True
        else:

            dist = abs(int(pos[i]) - int( pos[n]))
            if dist < radius and i not in g[n]:
                if i != n:
                    addRepulsion = True
        if addRepulsion:
            repulsion.update([i])
    return repulsion


def StoreRepulsion(g, radius=None):
    pos = {i:g.node[i]['pos'] for i in g.nodes()}
    repulsion = { i: StoreRepulsionAdjacency(g,i,pos,radius) for i in g.nodes() }
    nUpdated = 0
    # Make sure the repulsion sets are reflexive
    for src in repulsion:
        for dest in repulsion[src]:
            if src not in repulsion[dest]:
                repulsion[dest].update([src])
                nUpdated+=1
    return repulsion


def SubsampleRepulsion(gAdjList, repulsion, factor=1.0):
    nReduced = 0
    for n in gAdjList.keys():
        if len(gAdjList[n]) *factor < len(repulsion[n]):
            toRemove = int(len(repulsion[n]) - len(gAdjList[n])*factor)
            repList = sorted(list(repulsion[n]))
            repSub  = Set([])
            fraction= len(repulsion[n])/(len(gAdjList[n])*factor)
            for i in range(0,len(repulsion[n])):
                if len(repSub) * fraction < i:
                    repSub.update(Set([repList[i]]))
                else:
                    nReduced+=1
            repulsion[n] = repSub

    return nReduced

def NeighborSimilarity(gAdjList,src,dest):
    return len(gAdjList[src].intersection(gAdjList[dest]))


def GetCost(g, i, j, r, m):
    if j in g[i]:
        return 1
    else:
        if j in m[i]:
            return -1
        else:
            return 0
    
allNodeSet = Set(g.nodes())



def CountRepulsion(cut, nodes, repulsion):
    n=0
    for node in nodes:
        n+= len(repulsion[node].intersection(cut))
    return n


def CountEdgesOutside(gadj, cut):
    #
    # return the number of edges between nodes in a and b
    #
    
    n=0
    for src in cut:
        n += len(gadj[src].difference(cut))
    return n
        

def ScoreCut(gadj, cut, repulsion):
    #
    # Count repulsion edges inside the cut.
    #
    nNeg = CountRepulsion(cut, cut, repulsion)/2

    #
    # Count edges from the cut to other nodes.
    #
    nPos = CountEdgesOutside(gadj, cut)
    return nNeg*args.penalty + nPos

def Neighboring(g, subgraph, n):
    neighbors = Set(sorted(g[n].keys()))
    return list(neighbors.difference(subgraph))

def GrowSubgraph(g, subgraph, neighbors, n):
    subgraph.update(n)
    for i in n:
        neighbors.update(g[i].keys())
    neighbors.difference_update(subgraph)
    neighbors.difference_update(n)


def GetAdjacent(gadj, subgraph, exclude=None):
    adjacent = Set()
    for n in subgraph:
        #
        # Possibly faster method to check for neighbors.
        #
        adjacent.update(gadj[n])

    if exclude is not None:
        adjacent.difference_update(exclude)

    # Make sure adjacent is not self referencing
    adjacent.difference_update(subgraph)
    return adjacent

def CountUnion(a,b):
    if (len(a) < len(b)):
        smaller = a
        larger = b
    else:
        smaller = b
        larger = a
    n = 0
    for i in smaller:
        if i in larger:
            n+=1
    return n


def ScoreCutExpansion(g, gadj, cut, expansion, repulsion):

    #
    # Thrifty method to compute new score of cut when nodes from
    # expansion are added to cut.
    #
    numEdgesCrossingCut = 0

    for src in expansion:
        neighbors = Set(g[src].keys())
        neighbors.difference_update(cut)
        neighbors.difference_update(expansion)
        numEdgesCrossingCut += len(neighbors)

    numEdgesInExpansion = 0
    for src in expansion:
        numEdgesInExpansion += len(gadj[src].intersection(expansion))

    # 
    # 'black' edges that were previously from the cut to the
    # expansion.  The score will decrease by this.
    #

    numEdgesBetweenExpansionAndCut = 0
    for src in expansion:
        for dest in gadj[src]:
            if dest in cut:
                numEdgesBetweenExpansionAndCut += 1

    #
    # 'red' edges that were outside the cut but are now inside.
    #
    numNewIncludedRepulsion = CountRepulsion(cut, expansion, repulsion)
#    print "  Score cut expansion exp-tot {}\texp-repl {}\texp-cut {}\texp-exp{}".format(numEdgesCrossingCut , numNewIncludedRepulsion ,  numEdgesBetweenExpansionAndCut , numEdgesInExpansion)

    return numEdgesCrossingCut + numNewIncludedRepulsion*args.penalty - numEdgesBetweenExpansionAndCut - numEdgesInExpansion

def IsolatedCutScore(gadj, testCut):
    n=0
    for node in testCut:
        for dest in gadj[node]:
            if dest > node:
                n+=1
    return n


class GraphOutput:
    def __init__(self, fileName):
        self.fileName = fileName
        self.vertexColor = 1
        self.iter = 0

def TestCutExpansion(g, gAdjList, cut, repulsion, expandedCut):
    curScore = ScoreCut(gAdjList, expandedCut, repulsion)
    #
    # Find the cost of all of the neighbors in what will be expanded.
    #
    isolatedTestCutScore = ScoreCut(gAdjList, expandedCut, repulsion)

    #
    # Try and come up with a thrifty method to compute expansion.
    #
    cutExpansion = ScoreCutExpansion(g, gAdjList, cut, expandedCut, repulsion)
    newCutScore = curScore + cutExpansion
    return curScore, newCutScore

    

def GrowCut(g, cut, repulsion, minNeighborSimilarity=3):
    neighbors = GetAdjacent(g, cut)
    grewCut = True
    iter = 0
    nodes = g.nodes()
    adjList = g.adjacency_list()
    gAdjList = { nodes[i]: Set(adjList[i]) for i in range(0,len(nodes)) }
    it = 0

    curScore = ScoreCut(gAdjList, cut, repulsion)
    while grewCut:
        grewCut = False
        nSearched = 0
        searched = []
        neighbors = Set([])
        cutNodes = list(cut)
        random.shuffle(cutNodes)
        newCutScore = 0
        iter +=1
        for cutNode in cutNodes:
            cutNeighbors = gAdjList[cutNode].difference(cut)
            for n in cutNeighbors:

                if NeighborSimilarity(gAdjList, cutNode, n) < minNeighborSimilarity:
                    continue
                # 
                # Define the nodes that will be expanded
                #
                expandedCut = Set([n])# #O(n)
            
                #
                # When the expansion is a search, e.g. Set([n]+ gAdjList[n]), remove overlap with current cut
                #
                expandedCut.difference_update(cut) # O(n)
        
                #
                # Find the cost of all of the neighbors in what will be expanded.
                #
                isolatedTestCutScore =  ScoreCut(gAdjList, expandedCut, repulsion)

                #
                # Try and come up with a thrifty method to compute expansion.
                #
                cutExpansion = ScoreCutExpansion(g, gAdjList, cut, expandedCut, repulsion)

            
                #
                # Get the score of the entire new test cut.  Ideally this is equal to the cut expansion
                #
                #newCutScore = ScoreCut(gAdjList, testCut, repulsion) 
                newCutScore = curScore + cutExpansion

                #
                # The benefit of adding this cut is equal to the new test cut minus
                # all the nodes in isolation. I'm not sure if this is a good thing --
                # maybe the cost of the cut should be equal to the assumption that 
                # all other nodes are in the same cut.
                #
                newCutScoreDiff = newCutScore - isolatedTestCutScore

                it += 1
#                print str((newCutScore, isolatedTestCutScore, newCutScore - isolatedTestCutScore , curScore))
                if  newCutScore - isolatedTestCutScore < curScore:
                    GrowSubgraph(g, cut, neighbors, expandedCut)
                    grewCut = True
                    curScore = newCutScore
                    break
                nSearched+=1
            searched.append(str(n))
    return cut

class Input:
    def __init__(self,g,it,repulsion):
        self.g = g
        self.it = it
        self.repulsion = repulsion



def SampleCuts(g, nIter, repulsion):
    clustered = Set([])
    cuts = []
    it = 0
    g.graph['NumVertexColors'] = nIter+1
    unclustered = None
    while (len(clustered) < len(g.nodes()) and it < nIter):
        it +=1
        unclustered = list(Set(g.nodes()).difference(clustered))
        if len(unclustered) > 0:
            seed = unclustered[random.randint(0,len(unclustered)-1)]
            cut = Set([seed])
            GrowCut(g, cut, repulsion)
            clustered.update(cut)
            cuts.append(cut)
            print "iter: " + str(it) + "\t" + str(len(cuts[-1])) + "\t" + str(len(clustered)) + "/" + str(len(g.nodes()))
        else:
            break
    return (cuts, unclustered)
        
def SampleCutsIt(cl):
    return SampleCuts(cl.g, cl.it, cl.repulsion)

def ReciprocalOverlap(a,b,f):
    if f == 0:
        return True
    ovp = len(a.intersection(b))
    if ovp == 0:
        return False
    else:
        return float(ovp) / len(a) > f and float(ovp) / len(b) > f
    

def MergeCuts(cuts, ratio=0.9):
    #
    # Do a greedy merge of cuts.
    #
    i = 0
    cutMerged = True
    while cutMerged:
        cutMerged = False
        while i < len(cuts) - 1:
            j = i + 1
            while j < len(cuts):
                if ReciprocalOverlap(cuts[i], cuts[j], ratio):
                    cuts[i].update(cuts[j])
                    del cuts[j]
                    cutMerged = True
                else:
                    j+=1
            i+=1


# Now try swapping nodes that border two cuts

def GetCut(n, cuts):
    for i in range(0,len(cuts)):
        if n in cuts[i]:
            return i
    return None

#
# Check overlap between a set of nodes and all other cuts
#
def CountCutMembership(nodeSet, cuts):
    membership = { i: len(nodeSet.intersection(cuts[i])) for i in range(0,len(cuts)) }
    for cutIndex in membership.keys():
        if membership[cutIndex] == 0:
            del membership[cutIndex]
    return membership



def GetCutIndices(cuts):
    return  { j : i for i in range(0,len(cuts)) for j in cuts[i] }

def TestCutDeletion(adjList, repulsion, cuts, cutIndex, node):
    cuts[cutIndex].difference_update(Set([node]))            
    score = ScoreCut(adjList, cuts[cutIndex], repulsion)
    cuts[cutIndex].update(Set([node]))
    return score

def PickCutDeletion(gAdjList, repulsion, cuts, cutIndices, node):

    scores = [TestCutDeletion(gAdjList, repulsion, cuts, i, node) for i in cutIndices]
    lowestScore = min(scores)
    for minScoreIndex in range(0,len(scores)):
        if scores[minScoreIndex] == lowestScore:
            break
    for i in range(0,minScoreIndex):
        cuts[cutIndices.keys()[i]].difference_update(Set([node]))
    for i in range(minScoreIndex+1, len(cutIndices)):
        cuts[cutIndices.keys()[i]].difference_update(Set([node]))


def AssignNodesToUniqueClusters(gAdjList, repulsion, cuts):
    for node in gAdjList.keys():
        membership = CountCutMembership(Set([node]), cuts)
        if len(membership) > 1:
            PickCutDeletion(gAdjList, repulsion, cuts, membership, node)

def TestCutSwap(adjList, repulsion, cuts, cutScores, a, b, swap):
    cuts[a].difference_update(swap)
    cuts[b].update(swap)
    scoreA = ScoreCut(adjList, cuts[a], repulsion)
    scoreB = ScoreCut(adjList, cuts[b], repulsion)
    #
    # Put things back where they came from
    #
    cuts[a].update(swap)
    cuts[b].difference_update(swap)
    return (cutScores[a] + cutScores[b] - scoreA - scoreB, scoreA, scoreB)

def OptimizeBySwappingNodes(gAdjList, cuts, repulsion, nIter=0, minGain=10, neighborSimilarityCutoff=3):

    cutScores = [ScoreCut(gAdjList, cut, repulsion) for cut in cuts]    
    cutIndices = GetCutIndices(cuts)
    gain = 0
    swapMade = True
    iterIndex = 0
    while  (nIter == 0 or iterIndex < nIter) and swapMade:
        swapMade = False
        iterIndex+=1        
        for cutIt in range(0,len(cuts)):
            cut = cuts[cutIt]
            neighbors = ABPUtils.GetCutNeighbors(gAdjList, cut)
            maxNeighbor = None
            maxGain     = 0
            maxNeighborCut = 0
        
            for node in neighbors:
                if node not in cutIndices:
                    continue
                maxNeighborSimilarity = 0
                for cutNode in cut:
                    if node in gAdjList[cutNode]:
                        neighborSimilarity= NeighborSimilarity(gAdjList,cutNode, node)
                        maxNeighborSimilarity=max( neighborSimilarity , maxNeighborSimilarity)
                if maxNeighborSimilarity < neighborSimilarityCutoff:
                    continue
                            
                neighborCut = cutIndices[node]
                (swapScore, score0, score1) = TestCutSwap(gAdjList,
                                                          repulsion,
                                                          cuts,
                                                          cutScores,
                                                          neighborCut,
                                                          cutIt,
                                                          Set([node]))
                if swapScore > minGain:
                    print "Swapping " + str(node) + "\t" + str(cutIt) + "\t" + str(neighborCut) + "\t" + str(swapScore)
                    swapMade = True
                    cutScores[neighborCut]  = score0
                    cutScores[cutIt]        = score1
            
                    cuts[cutIt].update(Set([node]))
                    cuts[neighborCut].difference_update(Set([node]))
                    gain += swapScore

    return gain

def FilterConflictingNodes(gAdjList, repulsion, cuts, nodes, cutoff):
    cutScores = [ScoreCut(gAdjList, cut, repulsion) for cut in cuts]    
    cutIndices = GetCutIndices(cuts)
    gain = 0
    for node in nodes:
        if node not in cutIndices:
            continue
        neighborMembers = CountCutMembership(gAdjList[node], cuts)
        if len(neighborMembers) > 1:
            if min(neighborMembers.values()) > cutoff:
                print "Filtering " + str(node) + " " + str(min(neighborMembers.values()))
                cuts[cutIndices[node]].difference_update(Set([node]))

def RemoveSmallCuts(cuts, minCutSize=3):
    i = 0
    while i < len(cuts):
        if len(cuts[i]) < 2:
            cuts.remove(cuts[i])
        else:
            i+=1
    

def AddRepulsionEdges(g, repulsion):
    for src in repulsion.keys():
        for dest in repulsion[src]:
            g.add_edge(src, dest, weight=0.1, color=1)

adjList = g.adjacency_list()
nodes = list(g.nodes())
gAdjList = { n: Set(g[n].keys()) for n in g.nodes() }

repulsion = StoreRepulsion(g, args.radius)

nReduced = SubsampleRepulsion(gAdjList, repulsion, factor=args.factor)
print nReduced

step = Input(g, args.niter, repulsion)
(cuts, unclustered) = SampleCutsIt(step)
if args.layout:
    ABPUtils.ApplyLayout(g,args.layout)

MergeCuts(cuts,0.75)    
AssignNodesToUniqueClusters(gAdjList, repulsion, cuts)


cutLengths = [len(c) for c in cuts]
RemoveSmallCuts(cuts)

#
# Search for components not added to any set
#

def MaxOverlap(query, cuts):
    maxOverlap = 0
    for cut in cuts:
        maxOverlap = max(maxOverlap, len(query.intersection(cut)))
    return maxOverlap
                        
components = ABPUtils.GetComponents(g,4)

componentSets = [Set(c) for c in components]
for cc in componentSets:
    maxOverlap = MaxOverlap(cc, cuts)
    if maxOverlap < 4:
        print "adding component of size: " + str(len(cc))
        cuts.append(cc)


if args.swap > 0:
    OptimizeBySwappingNodes(gAdjList, cuts, repulsion, args.swap)


ABPUtils.ColorGraphByCut(g,cuts)
    
if args.plot is not None:
    g.graph['NumVertexColors'] = len(cuts)+1

    labels = {n: "" for n in g.nodes()}
    if args.plotRepulsion:
        repGraph = g.copy()
        AddRepulsionEdges(repGraph, repulsion)
        ABPUtils.DrawGraph(repGraph, args.plot, labels=labels)
    else:
        ABPUtils.DrawGraph(g, args.plot, labels=labels)

if args.cuts is not None:
    cutsFile = open(args.cuts,'w')

    for cut in cuts:
        cutsFile.write("\t".join([str(g.node[c]["index"]) for c in sorted(list(cut))]) + "\n")
    cutsFile.close()

if args.sites is not None:
    sitesFile = open(args.sites, 'w')
    for cut in cuts:
	sites = sorted([g.node[n]['pos'] for n in cut])
        sitesFile.write("\t".join(str(s) for s in sites ) + "\n")
    
if args.out is not None:
    allCuts = Set()
    for c in cuts:
        allCuts.update(c)
    notCut = Set(g.nodes()) - allCuts
    notCutList = sorted(list(notCut))
    g.remove_nodes_from(notCutList)
    for n in g.nodes():
        g.node[n]['color'] = -1
        
    for i in range(0,len(cuts)):
        for c in cuts[i]:
            g.node[c]['color']= i
    ABPUtils.WriteGraph(g, args.out)


if args.embed:
    IPython.embed()

#AddRepulsionEdges(g,repulsion)
#ABPUtils.DrawGraph(g, "test.repl.png")

