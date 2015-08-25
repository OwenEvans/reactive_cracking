from dolfin import *

# Define class marking boundaries
class TopDirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near((1 - x[1]), 0)

class NeumanBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0]*(1 - x[0]), 0)

class BottomDirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[1], 0)

# Create mesh and define function space
mesh = UnitSquareMesh(50, 50)
V = VectorFunctionSpace(mesh, 'DG', 2)

# Define test and trial functions
du = TrialFunction(V)
u_t = TestFunction(V)
u = Function(V)

# Define normal vector and mesh size
#n = FacetNormal(mesh)
#h = CellSize(mesh)
#h_avg = (h('+') + h('-'))/2

# Dirichlet bottom boundary condition (constant in time)
#u0 = Expression(('0.0', '5.0'))
#u1 = Expression(('0.0', '0.0'))

# Mark facets of the mesh
boundaries = FacetFunction('size_t', mesh, 0)
NeumanBoundary().mark(boundaries, 2)
TopDirichletBoundary().mark(boundaries, 1)
BottomDirichletBoundary().mark(boundaries, 3)

# Define outer surface measure
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

theta = 0.0
u_init = Constant(('0.0', '0.0'))
u_prev = project(u_init, V)
u_theta = (1 - theta)*u + theta*u_prev

# Define penalty parameters
alpha = 10.0
gamma = 5.0
sigma_c = Constant(1.0)
delta_c = 0.5
G = Constant(1.0)
Kon2G = Constant(1.0)
s = Constant(4.0)

# Stress computation
epsilon = sym(grad(u_theta))
epsilond = epsilon - tr(epsilon)/3.*Identity(2)
sigma = G*epsilond + Kon2G*tr(epsilon)*Identity(2)

# Time-stepping
t     = 0.0
t_max = 0.1
dt    = 0.01

file = File("displacement.pvd")

while (t < t_max):

  u0 = Constant(0.1)
  u1 = Constant(0.0)

  n = FacetNormal(mesh)
  h = CellSize(mesh)

  epsilon_t = sym(grad(u_t))
  epsilond_t = epsilon_t - tr(epsilon_t)/3.*Identity(2)
  sigma_t = G*epsilond_t + Kon2G*tr(epsilon_t)*Identity(2)
  sigma_t_n = dot(avg(sigma_t), n('+'))
  sigma_t_n_d = dot(sigma_t, n)

  sigma_n = dot(avg(sigma), n('+'))
  sigma_n_d = dot(sigma, n)

  epsilon_stab = outer(jump(u_theta), n('+'))
  epsilond_stab = epsilon_stab - tr(epsilon_stab)/3.*Identity(2)
  sigma_stab = G('+')*epsilon_stab + Kon2G('+')*tr(epsilon_stab)*Identity(2)

  # Define variational problem
  overstress = inner(n('+'), dot(avg(sigma), n('+'))) - sigma_c('+')
  D_trial = 0.5*(1. + tanh(overstress/0.5)) # do you need the half factor here?

  r = inner(epsilon_t, sigma)*dx \
      - (1. - D_trial)*inner(sigma_t_n, jump(u_theta))*dS \
      - (1. - D_trial)*inner(jump(u_t), sigma_n)*dS \
      + (1. - D_trial)*(s('+')/avg(h))*inner(outer(jump(u_t), n('+')), sigma_stab)*dS \
      - D_trial*inner(jump(u_t), n('+'))*sigma_c('+')*dS \
      - sigma_t_n_d[1]*u_theta[1]*ds(1) \
      - u_t[1]*sigma_n_d[1]*ds(1) \
      + (s/h)*u_t[1]*(2*G*u_theta[1]/3. + Kon2G*u_theta[1]/3.)*ds(1) \
      - sigma_t_n_d[1]*u_theta[1]*ds(3) \
      - u_t[1]*sigma_n_d[1]*ds(3) \
      + (s/h)*u_t[1]*(2*G*u_theta[1]/3. + Kon2G*u_theta[1]/3.)*ds(3) \
      - (s/h)*u0*u_t[1]*ds(1) \
      -(s/h)*u1*u_t[1]*ds(3)
  
  #r = inner(sigma(u), grad(v))*dx \
  #    - (1 - D_trial)*inner(t_star, jump(v))*dS \
  #    - dot(sigma(u)*n, v)*ds(1) \
  #    + (gamma/h)*inner(u, v)*ds(1) \
  #    - dot(sigma(u)*n, v)*ds(3) \
  #    + (gamma/h)*inner(u, v)*ds(3) \
  #    - (gamma/h)*inner(u0, v)*ds(1) - (gamma/h)*inner(u1, v)*ds(3) \
  #    -D_trial*sigma_c('+')*(1 - (1./delta_c))*inner(jump(u), n('+'))*inner(jump(v), n('+'))*dS

  # Compute jacobian
  J = derivative(r, u, du)

  # Solve
  #solve(r == 0, u, J=J, solver_parameters={"snes": {'linear_solver': "mumps", 'maximum_iterations': 25}})

  # Define the solver parameters
  snes_solver_parameters = {"nonlinear_solver": "snes",
                            "snes_solver"     : { "linear_solver"   : "lu",
                                                  "maximum_iterations": 20,
                                                  "report": True,
                                                  "error_on_nonconvergence": False,
                                                  }}

  # Set up the non-linear problem
  problem = NonlinearVariationalProblem(r, u, J=J)

  # Set up the non-linear solver
  solver  = NonlinearVariationalSolver(problem)
  solver.parameters.update(snes_solver_parameters)
  info(solver.parameters, True)
  
  solver.solve()
#print("Solution vector norm (0): {!r}".format(u.vector().norm("l2")))

  u_prev.assign(u)

  #print unorm
  # Plot solution
  #plot(u, interactive=True)
  
  file << u
  t += dt
  
