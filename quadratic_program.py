import numpy as np
import qiskit as qk
np.set_printoptions(linewidth=np.inf)#so matricies don't wrap as often
def compute_c(A):
    n = len(A)
    c = [0] * n   # initialize c as a list of n zeros
    for i in range(n):
        c[i] = sum(A[i]) / 2
    return np.array(c)
def compute_Q(A):
   #nah chief, this ain't right
    n = A.shape[0]
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = (A[i,j] - A[j,i]) / 4
    return Q
def cost(x,Q,c):
   return (np.matmul(np.matmul(x.T, Q), x) + np.matmul(c.T, x))
   
   
def info_matrix(qc, conv=complex):
   '''gives the matrix representation of a quantum circuit with the function "conv" applied to each of its elements'''
   print([conv(y) for y in [z for z in np.array(qk.quantum_info.Operator(qc))]])

qc=qk.QuantumCircuit(2)
qc2=qk.QuantumCircuit(2)
qc2.h(0)

qc.h(0)
qc.cx(0,1)
print(qk.quantum_info.Operator(qc))
qc.measure_all()
#print(info_matrix(qc))
backend = qk.Aer.get_backend('qasm_simulator')
# Execute the circuit on the simulator and get the result
result = qk.execute(qc, backend=backend).result()
# Print the counts of each measurement outcome
print(result.get_counts(qc))
print(qc)
if False:
   A=np.array([[0,1,0],[1,0,1],[0,1,0]])
   print(cost(np.array([1,1,0]),compute_Q(A),compute_c(A)),end="\n\n\n")
   print(compute_Q(A))
   print(compute_c(A))
   print(np.matmul(np.array((1,1,1)),compute_c(A).T))
