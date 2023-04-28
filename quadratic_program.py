import numpy as np
import qiskit as qk
np.set_printoptions(linewidth=np.inf)#so matricies don't wrap as often
from QCutting import createGraph
"""this program converts an adjacency matrix A into a quadratic optimization 
problem with linear part compute_c and quadratic matrix compute_Q"""
def compute_c(A):
    n = len(A)
    c = [0] * n   # initialize c as a list of n zeros
    for i in range(n):
        c[i] = sum(A[i])
    return np.array([c]).T
def compute_Q(A):
    n=len(A)
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = -A[i,j]
    return Q
def cost(x,Q,c):
   '''cost of a binary assignment x given quadratic coefficient matrix Q and linear coefficient vector c'''
   return (np.matmul(np.matmul(x.T, Q), x) + np.matmul(c.T, x))
if __name__=="__main__":
    from sys import argv
    from random import randint
    if len(argv)==1:
        A=np.array([[0,1,0],[1,0,1],[0,1,0]])
    else:
        from ast import literal_eval
        from networkx import adjacency_matrix
        with open(argv[1]) as f:
            points=literal_eval(f.read())
            G=createGraph(points)
            A=adjacency_matrix(G).todense()
    x=[randint(0,1) for _ in range(len(A))]
    print(f"x={x}")
    print(f"Q= {compute_Q(A)}")
    print(f"c={compute_c(A)}")
    print(cost(np.array([x]).T,compute_Q(A),compute_c(A))[0][0],end="\n\n\n")

