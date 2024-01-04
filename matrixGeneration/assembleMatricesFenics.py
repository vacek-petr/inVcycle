import sys
import petsc4py
import time
petsc4py.init(sys.argv)
from petsc4py import PETSc
opts = PETSc.Options()

import numpy as np
from scipy import sparse
import scipy.io
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

# right hand side
f = Expression("1.0", degree=1)

nlevel = 7 # number of levels
n = 20 # size coarsest level (points on a line)

matricesA = np.empty((nlevel,), dtype=object)
vectorsF = np.empty((nlevel,), dtype=object)
matricesP = np.empty((nlevel-1,), dtype=object)
vectorsBoundaryNodes = np.empty((nlevel,), dtype=object)

def boundary(x, on_boundary):
    return on_boundary


def operator(u,v):
    # Poisson problem 
    a = inner(Constant(1.0)*grad(u), grad(v))*dx(1) + inner(Constant(1.0)*grad(u), grad(v))*dx(3) + inner(Constant(1.0)*grad(u), grad(v))*dx(2) + inner(Constant(1.0)*grad(u), grad(v))*dx(4)
    
    # jump-1024 problem
    # a = inner(Constant(1024.0)*grad(u), grad(v))*dx(1) + inner(Constant(1024.0)*grad(u), grad(v))*dx(3) + inner(Constant(1.0)*grad(u), grad(v))*dx(2) + inner(Constant(1.0)*grad(u), grad(v))*dx(4)
    return a

def rightside(f,v):
    b = f*v*dx
    return b

    
tic = time.time()

for j in range(0,nlevel):

    mesh = UnitSquareMesh(n*(2**(j)),n*(2**(j)))

    leftBottomSubdomain = AutoSubDomain(lambda x, on_exterior: (x[0] <= 0.5)&(x[1] <= 0.5))
    leftTopSubdomain = AutoSubDomain(lambda x, on_exterior: (x[0] <= 0.5)&(x[1] >= 0.5))
    rightBottomSubdomain = AutoSubDomain(lambda x, on_exterior: (x[0] >= 0.5)&(x[1] <= 0.5))
    rightTopSubdomain = AutoSubDomain(lambda x, on_exterior: (x[0] >= 0.5)&(x[1] >= 0.5))
    cf = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    rightTopSubdomain.mark(cf, 1)
    leftTopSubdomain.mark(cf, 2)
    leftBottomSubdomain.mark(cf, 3)
    rightBottomSubdomain.mark(cf, 4)

    dx = Measure('dx', domain=mesh, subdomain_data=cf)
    
    plot(cf)
    plt.show()
    fileName = "subdomains_" + str(j) + ".png"
    plt.savefig(fileName)
    

    # plot meshes
    plot(mesh)
    plt.show()
    fileName = "mesh_" + str(j) + ".png"
    plt.savefig(fileName)
    
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
    
    # solve
    tic2 = time.time()
    u = Function(VFine)
    solve(A, u.vector(),F,"cg")
    plot(u)
    plt.show()
    fileName = "solutution_" + str(j) + ".png"
    plt.savefig(fileName)


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
print("Sucess!")