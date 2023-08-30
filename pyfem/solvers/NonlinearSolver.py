############################################################################
#  This Python file is part of PyFEM-1.0, released on Aug. 29, 2012.       #
#  The PyFEM code accompanies the book:                                    #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#  Comments and suggestions can be sent to:                                #
#     PyFEM-support@tue.nl                                                 #
#                                                                          #
#  The latest version can be downloaded from the web-site:                 #                                                                          
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

from pyfem.util.BaseModule import BaseModule
import numpy as np
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness

import sys

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class NonlinearSolver(BaseModule):

    def __init__(self, props, globdat):
        self.tol = 1.0e-3
        self.iterMax = 10
        self.maxCycle = sys.maxsize
        self.maxLam = 1.0e20

        super().__init__(props)

        globdat.lam = 0.0

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def run(self, props, globdat):
        globdat.cycle += 1

        globdat.lam = 1.0 * globdat.cycle

        dofCount = len(globdat.dofs)

        a = globdat.state
        Da = globdat.Dstate
        fhat = globdat.fhat

        Da[:] = np.zeros(dofCount)
        fint = np.zeros(dofCount)
        fext = globdat.lam * fhat

        print('=================================')
        print(' Load step %i' % globdat.cycle)
        print('=================================')
        print('  NR iter : L2-norm residual')

        globdat.iiter = 0

        K, fint = assembleTangentStiffness(props, globdat)

        error = 1.

        while error > self.tol:
            globdat.iiter += 1

            da = globdat.dofs.solve(K, fext - fint)

            Da[:] += da[:]
            a[:] += da[:]

            K, fint = assembleTangentStiffness(props, globdat)

            norm = np.linalg.norm(fext)

            if norm < 1.0e-16:
                error = np.linalg.norm(fext - fint)
            else:
                error = np.linalg.norm(fext - fint) / norm

            print('  Iter', globdat.iiter, ':', error)

            if globdat.iiter == self.iterMax:
                raise RuntimeError('Newton-Raphson iterations did not converge!')

        # Converged

        globdat.elements.commitHistory()

        Da[:] = np.zeros(len(globdat.dofs))

        globdat.fint = fint

        if globdat.cycle == self.maxCycle or globdat.lam > self.maxLam:
            globdat.active = False
