import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter
from random import random
import networkx as nx
from itertools import combinations_with_replacement
import argparse
from scipy.optimize import minimize

default_num_verticies=4

def createRandomScattering(num_points,scale=1,x_trans=0,y_trans=0):
   '''random scattering of num_points scaled by scale, and translated by x_trans and y_trans'''
   return [(random()*scale+x_trans,random()*scale+y_trans) for _ in range(num_points)]
   
def createTestClustering(clusters,max_num_points=9, max_dist_btw_clusters=4):
   out=[]
   for _ in range(clusters):
      out+=createRandomScattering(int(random()*max_num_points)+1,x_trans=max_dist_btw_clusters*random(),y_trans=max_dist_btw_clusters*random()) 
   return out
   
   
def distance(point1,point2):
   '''distance between two points (tuples) of arbitrary dimension'''
   if len(point1) != len(point2):
       raise ValueError("The two points must have the same number of dimensions.")     
   sum_of_squares = 0
   for i in range(len(point1)):
       diff = point1[i] - point2[i] #order doesn't matter
       sum_of_squares += diff ** 2
   return sum_of_squares**.5

#handle command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--points',type=list,help='List of points. To be given as a list of tuples eg [(1,2),(2,3),...,(2.1,.01)]. If empty list is given, random ')
args = parser.parse_args()

# Checking if the random option was given
if not args.points:
    points = createRandomScattering(default_num_verticies)#4 random points distributed across 0,1
else:
    points = args.points
      
#plt.scatter([x[0] for x in points],[x[1] for x in points])

#create graph
G = nx.Graph()
G.add_nodes_from(list(range(len(points))))#each node corresponds to a point
for origin,destination in combinations_with_replacement(range(len(points)),2): #create the complete graph
   if origin!=destination: #no self loops
      G.add_edge(origin,destination,weight=distance(points[origin],points[destination])) #connect each vertex to each other vertex with an edge where the weight is proportional to the distance between them. This is a relatively expensive operation (O(n^2))

# Adjacency is essentially a matrix which tells you which nodes are connected. This matrix is given as a sparse matrix, so we need to
# convert it to a dense matrix
adjacency = nx.adjacency_matrix(G).todense()
nqubits = len(points) #each qubit corresponds to a vertex in the graph and each boolean assignment for each point corresponds with what side of the cut the vertex is in. 

def mincut_obj(solution, graph):
    """Given a bit string as a solution, this function returns
    the sum of the weights of the edges between the sections.
       solution: solution bit string
        graph: networkx graph
    """
    if len(solution)!=len(graph.verticies()):
      raise ValueError("solution bit string must have the same length as the number of verticies in the graph.")
    obj = 0
    for i, j in graph.edges():
        if solution[i] != solution[j]:
            obj += graph.get_edge_data(i,j)["weight"]
    return obj

def compute_expectation(counts, graph):
    """Computes expectation value based on measurement results
    Args:
        counts: (dict) key as bit string, val as count
        graph: networkx graph
    Returns the expected value
    """
    avg = 0
    sum_count = 0
    for bit_string, count in counts.items():
        obj = mincut_obj(bit_string, graph)
        avg += obj * count
        sum_count += count
    return avg/sum_count


def createQAOACirc(G, theta):
    """
    Creates a parametrized qaoa circuit
    Args:  
        G: networkx graph
        theta: list unitary parameters
    Returns:
        qc: qiskit circuit
    """
    nqubits = len(G.nodes())
    p = len(theta)//2  # number of alternating unitaries
    qc = QuantumCircuit(nqubits)
    
    beta = theta[:p] #list of parameters for the problem layer
    gamma = theta[p:] #list of parameters for the mixing layer
    
    # initial state; uniform superposition
    for i in range(0, nqubits):
        qc.h(i)
    
    for irep in range(0, p):
        # problem unitary
        for pair in list(G.edges()): 
        #since we are always dealing with the complete graph, 
        #this will apply RZZ(\gamma) to each pair of qubits
            qc.rzz(2 * gamma[irep], pair[0], pair[1])
        qc.barrier()
        # mixer unitary
        for i in range(0, nqubits):
            qc.rx(2 * beta[irep], i)
        qc.barrier()  
    qc.measure_all()    
    return qc


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

#backend = Aer.get_backend('aer_simulator')
#backend.shots = 512

#qc_res = create_qaoa_circ(G, res.x)

#counts = backend.run(qc_res, seed_simulator=10).result().get_counts()

#plot_histogram(counts)