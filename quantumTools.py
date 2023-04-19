import numpy as np
class A(np.ndarray):
  def __new__(cls, input_array):
    obj = np.asarray(input_array).view(cls)
    return obj
  def __mod__(self, other):
    return np.kron(self, other)
  def __repr__(self):
    if type(self[0][0])==type(np.complex128(0)):
      pass
    return np.ndarray.__repr__(self)
      
#@=@ numpy.matmul. Whenever you want to multiply matrices, use this 
#whenever you want to kronecker product things use %
j=complex("0+1j")
one=A([[0,1]]).T
zero=A([[1,0]]).T
#pauli x y and z gates
Z=A([[1,0],[0,-1]])
Y=A([[0,-j],[j,0]])
X=A([[0,1],[1,0]])
I=A([[1,0],[0,1]])
H=(1/(2**.5))*A([[1,1],[1,-1]])
CX=(zero@zero.T)%I+(one@one.T)%X
CY=(zero@zero.T)%I+(one@one.T)%Y
CZ=(zero@zero.T)%I+(one@one.T)%Z
#basis vectors

#zerozero=zero%zero
plus=H@one
minus=H@zero


