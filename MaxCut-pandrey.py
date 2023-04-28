from QCutting import distance, printClustering, createGraph, file_to_list
from networkx import get_edge_attributes
import sys
from maxcut import MaxCutSDP
G=createGraph(file_to_list(sys.argv[1]))
sol=MaxCutSDP(G)
sol.solve(verbose=False)
printClustering(G,[bool(x+1) for x in sol.get_results('cut')])
