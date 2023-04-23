import numpy as np
from latex2sympy2 import latex2sympy as l2s
import sympy
class A(np.ndarray):
  def __new__(cls, input_array):
    obj = np.asarray(input_array).view(cls)
    return obj
  def __mod__(self, other):
    return np.kron(self, other)
  def __eq__(self,other):
    return np.array_equal(self,other)
  def __repr__(self):
    if type(self[0][0])==type(np.complex128(0)):
      pass
    return np.ndarray.__repr__(self)
   
def bmatrix(A:A):
  """formats to LaTeX bmatrix"""
  if len(A.shape) > 2:
    raise ValueError('bmatrix can at most display two dimensions')
  lines = str(A).replace('[', '').replace(']', '').splitlines()
  rv = [r'\begin{bmatrix}']
  rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
  rv +=  [r'\end{bmatrix}']
  return '\n'.join(rv)
def texify(s:str,replace_gates=True):
  if replace_gates:
    gates={\
      "Z":A([[1,0],[0,-1]]),
      "Y":A([[0,-j],[j,0]]),
      "X":A([[0,1],[1,0]]),
      "I":A([[1,0],[0,1]]),
      "H":(1/(2**.5))*A([[1,1],[1,-1]]),
      "CX":(zero@zero.T)%I+(one@one.T)%X,
      "CY":(zero@zero.T)%I+(one@one.T)%Y,
      "CZ":(zero@zero.T)%I+(one@one.T)%Z}
    for gate in gates.keys():
      s=s.replace(gate,bmatrix(gates[gate]))
  vects={
    "zero.T%zero.T":r"\bra{00}",
    "zero.T%one.T":r"\bra{01}",
    "one.T%zero.T":r"\bra{10}",
    "one.T%one.T":r"\bra{11}",
    "one%one":r"\ket{11}",
    "zero%one":r"\ket{01}",
    "one%zero":r"\ket{10}",
    "zero%zero":r"\ket{00}",
    "one.T":r"\bra{1}",
    "zero.T":r"\bra{0}",
    "one":r"\ket{1}",
    "zero":r"\ket{0}",
    "plus":r"\ket{+}",
    "minus":r"\ket{-}"}    
  for vect in vects.keys():
    s=s.replace(vect,vects[vect])
  s=s.replace("%",r"\otimes ")
  s=s.replace("@","")
  return s
#def tex2math(wires:int, tex:str):
#  tex=tex.replace("Z_i", bmatrix(Zi(wires,i))
#@=@ numpy.matmul. Whenever you want to multiply matrices, use this 
#whenever you want to kronecker product things use %
#basis vectors
j=complex("0+1j")
one=A([[0,1]]).T
zero=A([[1,0]]).T
def basis(bit_string:str):
  text=[]
  for i in bit_string:
    if i=="1":
      text.append(one)
    else:
      text.append(zero)
  for i in range(len(text)):
    if i==0:
      out=text[i]
    else:
      out=out%text[i]
  return out
zz=basis("00")
zo=basis("01")
oz=basis("10")
oo=basis("11")
#pauli x y and z gates
Z=A([[1,0],[0,-1]])
Y=A([[0,-j],[j,0]])
X=A([[0,1],[1,0]])
I=A([[1,0],[0,1]])
H=(1/(2**.5))*A([[1,1],[1,-1]])
plus=H@one
minus=H@zero
CX=(zero@zero.T)%I+(one@one.T)%X
CY=(zero@zero.T)%I+(one@one.T)%Y
CZ=(zero@zero.T)%I+(one@one.T)%Z
CCX=((zero%zero)@(zero.T%zero.T))%I+(zero%one)@(zero.T%one.T)%I+(one%zero)@(one.T%zero.T)%I+((one%one)@(one.T%one.T))%X
def gateOnWire(i:int,wires:int,gate:A):
  if 0>=i<wires:
    raise(ValueError("invalid qubit index"))
  for wire in range(wires):
    if wire==0:
      if wire==i:
        out=gate
      else:
        out=I
    elif wire==i:
      out=out%gate
    else:
      out=out%I
  return out
def Zi(i,wires):
  return gateOnWire(i,wires,Z)
def Yi(i,wires):
  return gateOnWire(i,wires,Y)
def Xi(i,wires):
  return gateOnWire(i,wires,X)
def Hi(i,wires):
  return gateOnWire(i,wires,H)



