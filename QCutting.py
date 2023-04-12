import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter
from random import random
import networkx as nx
from itertools import combinations_with_replacement

G = nx.Graph()
#random scattering of points to cluster
num_points=10
max_dimensions=100 #exclusive
points=[]
def distance(origin,destination): 
   """returns distance between two points of arbitrary dimension"""
   accumulator=0
   for dimension in zip(origin,destination):
      accumulator+=(dimension[0]-dimension[1])**2 #order doesn't matter
   return accumulator**.5
      
for i in range(num_points):
   points.append((int(random()*max_dimensions),int(random()*max_dimensions)))
#plt.scatter([x[0] for x in points],[x[1] for x in points])
G.add_nodes_from(list(range(num_points)))#each node corresponds to a point
for origin,destination in combinations_with_replacement(range(num_points), 2):
   if origin!=destination:
      G.add_edge(origin,destination,weight=distance(points[origin],points[destination])) #connect each vertex to each other vertex with an edge where the weight is proportional to the distance between them. This is a relatively expensive operation (O(n^2))
#nx.draw(G, with_labels=True, alpha=0.8, node_size=500)
#plt.show()
# Adjacency is essentially a matrix which tells you which nodes are
# connected. This matrix is given as a sparse matrix, so we need to
# convert it to a dense matrix
adjacency = nx.adjacency_matrix(G).todense()
nqubits = 4
if(False):
   beta = Parameter("$\\beta$")
   qc_mix = QuantumCircuit(nqubits)
   for i in range(0, nqubits):
       qc_mix.rx(2 * beta, i)
       
   qc_mix.draw()
   def create_qaoa_circ(G, theta):
       
       """
       Creates a parametrized qaoa circuit
       
       Args:  
           G: networkx graph
           theta: list
                  unitary parameters
                        
       Returns:
           qc: qiskit circuit
       """
       
       nqubits = len(G.nodes())
       p = len(theta)//2  # number of alternating unitaries
       qc = QuantumCircuit(nqubits)
       
       beta = theta[:p]
       gamma = theta[p:]
       
       # initial_state
       for i in range(0, nqubits):
           qc.h(i)
       
       for irep in range(0, p):
           
           # problem unitary
           for pair in list(G.edges()):
               qc.rzz(2 * gamma[irep], pair[0], pair[1])

           # mixer unitary
           for i in range(0, nqubits):
               qc.rx(2 * beta[irep], i)
               
       qc.measure_all()
           
       return qc

   # Finally we write a function that executes the circuit on the chosen backend
   def get_expectation(G, shots=512):
       
       """
       Runs parametrized circuit
       
       Args:
           G: networkx graph
           p: int,
              Number of repetitions of unitaries
       """
       
       backend = Aer.get_backend('qasm_simulator')
       backend.shots = shots
       
       def execute_circ(theta):
           
           qc = create_qaoa_circ(G, theta)
           counts = backend.run(qc, seed_simulator=10, 
                                nshots=512).result().get_counts()
           
           return compute_expectation(counts, G)
       
       return execute_circ
       
       from qiskit.visualization import plot_histogram

   backend = Aer.get_backend('aer_simulator')
   backend.shots = 512

   qc_res = create_qaoa_circ(G, res.x)

   counts = backend.run(qc_res, seed_simulator=10).result().get_counts()

   plot_histogram(counts)
