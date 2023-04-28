import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, Aer, execute
#!from qiskit.circuit import Parameter
from random import random
import networkx as nx
from itertools import combinations_with_replacement
import argparse
from scipy.optimize import minimize
from numpy import set_printoptions
set_printoptions(linewidth=500)#so matricies don't wrap as often

def createRandomScattering(num_points,scale=1,x_trans=0,y_trans=0):
   '''random scattering of num_points scaled by scale, and translated by x_trans and y_trans'''
   return [(random()*scale+x_trans,random()*scale+y_trans) for _ in range(num_points)]
   

def createTestClustering(clusters,max_num_points=4, max_dist_btw_clusters=4):
   '''create a set of points with a number of clustered together scatterings'''
   out=[]
   for _ in range(clusters):
      #this only looks complicated because of the beefy variable names, the +1 is to stop the cluster from being of size 0
      out+=createRandomScattering(int(random()*max_num_points)+1,x_trans=max_dist_btw_clusters*random(),y_trans=max_dist_btw_clusters*random()) 
   f=open("lastRandDist.txt", "w")
   f.write(str(out))
   f.close()
   return out
   
def createGraph(points):
   G = nx.Graph()
   for i in range(len(points)):
      G.add_node(i,pos=points[i])
   for origin,destination in combinations_with_replacement(range(len(points)),2): #create the complete graph
      if origin!=destination: #no self loops
         G.add_edge(origin,destination,weight=distance(points[origin],points[destination])) #connect each vertex to each 
         #other vertex with an edge where the weight is proportional to the distance between them.
   return G
def distance(point1,point2):
   '''distance between two points (tuples) of arbitrary dimension'''
   if len(point1) != len(point2):
      raise ValueError("The two points must have the same number of dimensions.")     
   sum_of_squares = 0
   for i in range(len(point1)):
      diff = point1[i] - point2[i] #order doesn't matter
      sum_of_squares += diff ** 2
   return sum_of_squares**.5





def printClustering(G,solution):
   '''prints a networkxgraph with associated positions colored according to a solution bitstring'''
   colors=[]
   for i in solution:
      if int(i):
         colors.append("blue")
      else:
         colors.append("pink")
   nx.draw_networkx(G,nx.get_node_attributes(G,'pos'), node_color=colors,with_labels=True,edgelist=[])
   
def printGraphWithWeights(G,solution):
   labels=nx.get_edge_attributes(G,'weight')
   labels= {k:str(v)[:4] for k, v in labels.items()}
   nx.draw_networkx_edge_labels(G,nx.get_node_attributes(G,'pos'),edge_labels=labels,bbox=dict(alpha=0))
   printClustering(G,solution)


def mincut_obj(solution, graph):
    """Given a bit string as a solution, this function returns
    the sum of the weights of the edges between the sections.
       solution: solution bit string
       graph: networkx graph
    """
    if len(solution)!=len(graph.nodes()):
       raise ValueError("solution bit string must have the same length as the number of verticies in the graph.")
    obj = 0
    for i, j in graph.edges():
       if solution[i] != solution[j]:
          obj -= graph.get_edge_data(i,j)["weight"] #Q_ij=-W_ij
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
    nqubits = len(G.nodes()) #each qubit corresponds to a vertex in the graph and each boolean assignment for each point corresponds with 
   #what side of the cut the vertex is in. 
    p = len(theta)//2  # number of alternating unitaries
    qc = QuantumCircuit(nqubits)
    
    beta = theta[:p] #list of parameters for the problem layers
    gamma = theta[p:] #list of parameters for the mixing layers
    
    # initial state; uniform superposition
    for i in range(0, nqubits):
        qc.h(i)
    
    for irep in range(0, p):
        # problem unitary
        for i,j in list(G.edges()): 
        #since we are always dealing with the complete graph, 
        #this will apply RZZ(\gamma) to each pair of qubits
            qc.rzz(G.get_edge_data(i,j)["weight"] * gamma[irep], i,j)
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
        
        qc = createQAOACirc(G, theta)
        counts = backend.run(qc, seed_simulator=10,nshots=512).result().get_counts()
        
        return compute_expectation(counts, G)
    
    return execute_circ
def file_to_list(f): #takes a textiowrapper aka the thing returned by open("filename.txt")
    from ast import literal_eval
    with open(f) as f:
         return literal_eval(f.read())
def plot_counts(counts):
      #we need to sort to get the bar chart in decent order
   sorted_counts = sorted(counts.items())
   # extract the keys and values from the sorted dictionary
   keys = [item[0] for item in sorted_counts]
   values = [item[1] for item in sorted_counts]
   # create a bar chart from the sorted dictionary
   plt.bar(keys, values, align='center')
   plt.show()

if __name__=="__main__":
   #handle command line arguments
   parser = argparse.ArgumentParser()
   parser.add_argument('-p', '--points',type=str,help='List of points. To be given as the path to a file which contains a list of tuples eg [(1,2),(2,3),...,(2.1,.01)]. If no points are given, random')
   parser.add_argument('-l', '--layers',type=str,help="Number of layers to apply to the Quantum Circuit. More is more expensive but more accurate.")
   args = parser.parse_args()

   # Checking if the random option was given
   if not args.points:
      points = createTestClustering(2) #args.clusters)#2 random clusters
   else:
      points=file_to_list(args.points)
   if not args.layers:
      layers=1
   else:
      layers=1*layers

   #create graph  
   G=createGraph(points)   
   res = minimize(get_expectation(G),[1.0, 1.0]*(layers),method='COBYLA') #get_expectation returns the function 
   #execute_circ which takes in a 1D vector of even dimension. 
   backend = Aer.get_backend('qasm_simulator')
   backend.shots = 512
   qc_res = createQAOACirc(G, res.x)
   qc_res.draw(output="mpl")
   plt.show()
   
   counts =backend.run(qc_res, seed_simulator=10).result().get_counts()
   plot_counts(counts)
   #show maximum cut
   bestSolution= max(counts, key=counts.get)
   printClustering(G,bestSolution)
   # printGraphWithWeights(G,bestSolution)
   plt.show()



