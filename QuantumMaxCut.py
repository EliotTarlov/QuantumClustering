from QCutting import get_expectation,createQAOACirc,file_to_list,createGraph,printGraphWithWeights,plot_counts,printClustering
import sys
from networkx import from_numpy_array
from numpy import array
from scipy.optimize import minimize
from qiskit import Aer
import matplotlib as mpl
import matplotlib.pyplot as plt
save_to_latex=True



G=createGraph(file_to_list(sys.argv[1]))
G.remove_edge(0,2)
#printGraphWithWeights(G,"000")
#G=from_numpy_array(array(file_to_list(sys.argv[1])))
res = minimize(get_expectation(G),[1.0, 1.0],method='COBYLA')
qc_res = createQAOACirc(G, res.x)
backend = Aer.get_backend('qasm_simulator')
backend.shots = 512
counts =backend.run(qc_res, seed_simulator=10).result().get_counts()
#plot_counts(counts)

sorted_counts = sorted(counts.items())
# extract the keys and values from the sorted dictionary
keys = [item[0] for item in sorted_counts]
values = [item[1] for item in sorted_counts]
# create a bar chart from the sorted dictionary
if save_to_latex: #save figure for latex
   mpl.use("pgf")
   mpl.rcParams.update({
       "pgf.texsystem": "pdflatex",
       'font.family': 'serif',
       'text.usetex': True,
       'pgf.rcfonts': False,
   })
#plt.bar(keys, values, align='center')
#plt.savefig('histogram.pgf')
# printGraphWithWeights(G,bestSolution)
   

#show maximum cut
bestSolution= max(counts, key=counts.get)
printClustering(G,bestSolution)
plt.savefig("cut.pgf")
