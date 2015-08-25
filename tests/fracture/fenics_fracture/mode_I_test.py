from dolfin import *

# some parameters
E     = Constant(100000.0) # Young's Modulus
nu    = Constant(0.0) # Poisson's ratio
gamma = 1.0 # Penalty factor

# define mesh
mesh = RectangleMesh(0.0, 0.0, 1.0, 2.0, 100, 2)

# define function space
U = VectorFunctionSpace(mesh, "DG", 2) # displacement

# define boundaries
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 2.0 - DOLFIN_EPS and on_boundary

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

# Mark facets of the mesh                                                                                                                                                   
#boundaries = FacetFunction('size_t', mesh, 0)
#Top().mark(boundaries, 1)

# mark the subdomains                                                                                                                                                                                       
top = Top()
bottom = Bottom()
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
top.mark(sub_domains, 1)
bottom.mark(sub_domains, 0)


# Define outer surface measure aware of Dirichlet and Neumann boundaries                                                                                                                                      
#ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
#top, bottom = Top(), Bottom()

# boundary conditions
u_top    = Expression(("0.0", "100.0"))
u_bottom = Constant((0.0, 0.0))

# stress computation
def sigma(v):
    return E*sym(grad(v))

# define variational problem
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
    u_top     = Constant(("0.0", "100.0"))
    bc_top    = DirichletBC(U, u_top, top)
    bc_bottom = DirichletBC(U, u_bottom, bottom)
    #bcs = [bc_top, bc_bottom]
    
    # cell normal and size
    n = FacetNormal(mesh)
    h = CellSize(mesh)
    h_avg = (h('+') + h('-'))/2

    # weak form 
    t_star = dot(avg(sigma(u)), n('+')) - (gamma/h_avg)*jump(u)
    r = inner(sigma(u), grad(u_t))*dx - inner(jump(u), avg(sigma(u_t))*n('+'))*dS - inner(avg(sigma(u))*n('+'), jump(u_t))*dS + (E*gamma/avg(h))*inner(jump(u_t), jump(u))*dS - inner(u, sigma(u_t)*n)*ds - inner(sigma(u)*n, u_t)*ds + (E*gamma/h)*inner(u, u_t)*ds
    
    
    
    #+ inner(u_top, inner(sym(grad(u_t)), n))*ds(1) - inner(u_t, inner(sym(grad(u)), n))*ds(1) + 


    #(gamma/h)*dot(u_t, u)*ds(1)#- inner(jump(u_t), n('+'))*0.5*(abs(inner(-t_star, n('+'))) - inner(-t_star, n('+')))*D*dS 
    #L = -inner(u_top, dot(sym(grad(u_t)), n))*ds(1) 

    # compute jacobian
    J = derivative(r, u, du)

    # Solve
    solve(r == 0, u, J=J, solver_parameters={"newton_solver": {'absolute_tolerance': 1E-8, 'relative_tolerance': 1E-8, 'maximum_iterations': 25, 'relaxation_parameter': 1.0}})

    file << u
    t += dt
