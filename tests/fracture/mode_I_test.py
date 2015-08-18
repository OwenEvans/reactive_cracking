from dolfin import *

# some parameters
E  = Constant(1.0) # Young's Modulus
nu = Constant(0.0) # Poisson's ratio

# define mesh
mesh = RectangleMesh(0.0, 0.0, 1.0, 2.0, 100, 2)

# define function space
U = VectorFunctionSpace(mesh, "CG", 2) # displacement

# define boundary for Dirichlet condition
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 2.0 - DOLFIN_EPS and on_boundary

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

top    = Top()
bottom = Bottom()

# Dirichlet BC ... this is time dependent, need to think
u_top    = Expression(("0.0", "t + 1.0"), t=0.0)
u_bottom = Constant((0.0, 0.0))

# Stress computation
def sigma(v):
    return E*sym(grad(v))

# Define variational problem
du  = TrialFunction(U)
u_t = TestFunction(U)
u   = Function(U)

t     = 0.0
t_max = 0.1
dt    = 0.01

# export
file = File("displacement.pvd")

while (t < t_max):
    
    # update Dirichlet bc
    u_top = Expression(("0.0", "t + 1.0"), t=t)
    bc_top    = DirichletBC(U, u_top, top)
    bc_bottom = DirichletBC(U, u_bottom, bottom)
    bcs = [bc_top, bc_bottom]

    # weak form (maybe put gravity in if this doesn't work straight off)
    r = inner(sym(grad(u_t)), sigma(u))*dx

    # compute jacobian
    J = derivative(r, u, du)

    # Solve
    solve(r == 0, u, J=J, bcs=bcs, solver_parameters={"newton_solver": {'absolute_tolerance': 1E-8, 'relative_tolerance': 1E-8, 'maximum_iterations': 25, 'relaxation_parameter': 1.0}})

    print (u.vector().array())

    file << u
    t += dt