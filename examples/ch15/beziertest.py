import numpy as np
import matplotlib.pyplot as plt

C = np.zeros(shape=[4, 4, 4])

C[0, 0, 0] = 1.0
C[0, 1, 1] = 1.0
C[0, 1, 2] = 0.5
C[0, 1, 3] = 0.25
C[0, 2, 2] = 0.5
C[0, 2, 3] = 0.5
C[0, 3, 3] = 0.25

C[1, 0, 0] = 0.25
C[1, 1, 0] = 0.5
C[1, 1, 1] = 0.5
C[1, 2, 0] = 0.25
C[1, 2, 1] = 0.5
C[1, 2, 2] = 1.0
C[1, 2, 3] = 0.5
C[1, 3, 3] = 0.5

C[2, 0, 0] = 0.5
C[2, 1, 0] = 0.5
C[2, 1, 1] = 1.0
C[2, 2, 2] = 1.0
C[2, 3, 3] = 1.0

C[3, 0, 0] = 1.0
C[3, 1, 1] = 1.0
C[3, 2, 2] = 1.0
C[3, 3, 3] = 1.0

coords = np.zeros(shape=(10, 2))

coords[0, :] = [0.0, 0.0]
coords[1, :] = [1.0, 1.0]
coords[2, :] = [1.0, 3.0]
coords[3, :] = [2.0, 3.0]
coords[4, :] = [2.5, 1.5]
coords[5, :] = [1.5, 0.5]
coords[6, :] = [3.0, 0.0]
coords[7, :] = [3.2, 2.0]
coords[8, :] = [3.8, 2.5]
coords[9, :] = [4.0, 0.0]

elems = np.zeros(shape=(4, 4), dtype=int)

elems[0, :] = [0, 1, 2, 3]
elems[1, :] = [1, 2, 3, 4]
elems[2, :] = [3, 4, 5, 6]
elems[3, :] = [6, 7, 8, 9]

output = []
length = 0.

from pyfem.util.BezierShapeFunctions import getElemBezierData

for elemNodes, Celem in zip(elems, C):

    sdata = getElemBezierData(coords[elemNodes, :], Celem, order=100, elemType='Line4')

    for idata in sdata:

        x = np.dot(idata.h, coords[elemNodes, :])

        output.append(x)

        length += idata.weight

print("The length of the curve is", length)

plt.plot([x[0] for x in output], [x[1] for x in output], '-')
plt.plot([x[0] for x in coords], [x[1] for x in coords], 'ro-')

plt.show()
