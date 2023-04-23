from QCutting import distance, printClustering, createGraph
from networkx import get_edge_attributes
import sys
from maxcut import MaxCutSDP
from ast import literal_eval
with open(sys.argv[1]) as f:
  points=literal_eval(f.read())
G=createGraph(points)
sol=MaxCutSDP(G)
sol.solve(verbose=False)
printClustering(G,[bool(x+1) for x in sol.get_results('cut')])
