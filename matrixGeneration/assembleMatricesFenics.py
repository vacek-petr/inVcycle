import sys
import petsc4py
import time
petsc4py.init(sys.argv)
from petsc4py import PETSc
print = PETSc.Sys.Print
opts = PETSc.Options()

import numpy as np
from scipy import sparse
import scipy.io
from dolfin import *
from mshr import *
from matplotlib import pyplot
import dolfin.common.plotting as fenicsplot


# peak 2D BJR1995
uex = Expression("10*sin(2*pi*x[0])*sin(pi*x[1])*exp(pow((x[0]-frac34),2) + pow((x[1]-frac34),2))", degree=4, pi=np.pi, frac34 = 0.75)
grad_uex = Expression( ["10*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(2*pi*x[0])*sin(pi*x[1])*(2*x[0] -frac32) + 20*pi*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*cos(2*pi*x[0])*sin(pi*x[1])","10*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(pi*x[1])*sin(2*pi*x[0])*(2*x[1] -frac32) + 10*pi*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*cos(pi*x[1])*sin(2*pi*x[0])"], degree=4, pi=np.pi, frac34=0.75,frac32 = 1.5)
f = Expression("50*pi*pi*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(pi*x[1])*sin(2*pi*x[0]) - 10*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(pi*x[1])*sin(2*pi*x[0])*pow((2*x[1] -frac32),2) - 10*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(pi*x[1])*sin(2*pi*x[0])*pow((2*x[0] -frac32),2) - 40*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*sin(pi*x[1])*sin(2*pi*x[0]) - 20*pi*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*cos(pi*x[1])*sin(2*pi*x[0])*(2*x[1] -frac32) - 40*pi*exp(pow((x[0] - frac34),2) + pow((x[1] - frac34),2))*cos(2*pi*x[0])*sin(pi*x[1])*(2*x[0] -frac32)",degree=4, pi=np.pi, frac34=0.75,frac32=1.5)



nlevel = 8 # number of levels
n = 11 # size coarsest level (points on a line)

matricesA = np.empty((nlevel,), dtype=np.object)
vectorsF = np.empty((nlevel,), dtype=np.object)
matricesP = np.empty((nlevel-1,), dtype=np.object)
vectorsBoundaryNodes = np.empty((nlevel,), dtype=np.object)

def boundary(x, on_boundary):
    return on_boundary

def operator(u,v):
    a = inner(grad(u), grad(v))*dx
    return a

def rightside(f,v):
    b = f*v*dx
    return b

    
tic = time.time()

for j in range(0,nlevel):
    mesh = UnitSquareMesh(n*(2**(j)),n*(2**(j)))
    
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    VFine = FunctionSpace(mesh, P1)
    print("level " + str(j) + ", number of DoF " + str(VFine.dim()))
    u = TrialFunction(VFine)
    v = TestFunction(VFine)
    bc = DirichletBC(VFine, 0.0, boundary)

    boundaryNodes = np.fromiter(bc.get_boundary_values().keys(), dtype=int)
    boundaryNodes = boundaryNodes + np.ones(boundaryNodes.shape)
    boundaryNodes.sort()
    vectorsBoundaryNodes[j] = boundaryNodes

    A, F = assemble_system(operator(u,v),rightside(f,v),bc)
    
    tic2 = time.time()
    u = Function(VFine)
    solve(A, u.vector(),F,"cg")

    discretizetionErrorEnergyNorm = sqrt(assemble(inner(grad(u)-grad_uex,grad(u)-grad_uex)*dx(mesh)))
    print("Energy norm of discretization error = {0:16.8e}".format(discretizetionErrorEnergyNorm))
    print("Solution time " + str(time.time()-tic2))
    print("-------------------------------------------")

    APetsc = as_backend_type(A).mat()
    indptr, indices, data = APetsc.getValuesCSR()
    AScipy = sparse.csr_matrix((data,indices,indptr))
    matricesA[j] = AScipy
    
    FPetsc = as_backend_type(F).vec()
    vectorsF[j] = FPetsc 

    if (j != 0):
        P = PETScDMCollection.create_transfer_matrix(VCoarse, VFine)
        PPetsc = P.mat()
        indptr, indices, data = PPetsc.getValuesCSR()
        PScipy = sparse.csr_matrix((data,indices,indptr))
        matricesP[j-1] = PScipy

    VCoarse = VFine

scipy.io.savemat("A.mat", {"A": matricesA})
scipy.io.savemat("F.mat", {"F": vectorsF})
scipy.io.savemat("P.mat", {"P": matricesP})
scipy.io.savemat("BN.mat", {"BN": vectorsBoundaryNodes})

toc = time.time()

print("Total time", toc-tic)