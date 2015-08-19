from dolfin import *

# FIXME: Make mesh ghosted
parameters["ghost_mode"] = "shared_facet"

# Define class marking Dirichlet boundary (x = 0 or x = 1)
class TopDirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near((1 - x[1]), 0)

# Define class marking Neumann boundary (y = 0 or y = 1)
class NeumanBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[0]*(1 - x[0]), 0)

class BottomDirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary and near(x[1], 0)

# Create mesh and define function space
mesh = UnitSquareMesh(24, 24)
V = VectorFunctionSpace(mesh, 'DG', 2)

# Define test and trial functions
u = TrialFunction(V)
v = TestFunction(V)

# Define normal vector and mesh size
n = FacetNormal(mesh)
h = CellSize(mesh)
h_avg = (h('+') + h('-'))/2

# Define the source term f, Dirichlet term u0 and Neumann term g
f = Expression(('0.0','0.0'))
u0 = Expression(('0.0', '10.0'))
g = Expression(('0.0', '0.0'))
u1 = Expression(('0.0', '0.0'))

# Mark facets of the mesh
boundaries = FacetFunction('size_t', mesh, 0)
NeumanBoundary().mark(boundaries, 2)
TopDirichletBoundary().mark(boundaries, 1)
BottomDirichletBoundary().mark(boundaries, 3)

# Define outer surface measure aware of Dirichlet and Neumann boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define parameters
alpha = 4.0
gamma = 8.0

def sigma(v):
  return sym(grad(v))

# Define variational problem
a = inner(sigma(u), grad(v))*dx \
   - inner(jump(u), avg(sigma(v))*n('+'))*dS \
   - inner(avg(sigma(u))*n('+'), jump(v))*dS \
   + alpha/h_avg*inner(jump(u), jump(v))*dS \
   - inner(u, sigma(v)*n)*ds(1) \
   - dot(sigma(u)*n, v)*ds(1) \
   + (gamma/h)*inner(u, v)*ds(1) \
   - inner(u, sigma(v)*n)*ds(3) \
   - dot(sigma(u)*n, v)*ds(3) \
   + (gamma/h)*inner(u, v)*ds(3)

L = inner(v, f)*dx + inner(g, v)*ds(2) - inner(u0, sigma(v)*n)*ds(1) + (gamma/h)*inner(u0, v)*ds(1) - inner(u1, sigma(v)*n)*ds(3) + (gamma/h)*inner(u1, v)*ds(3) 

# Compute solution
u = Function(V)
solve(a == L, u)
print("Solution vector norm (0): {!r}".format(u.vector().norm("l2")))

# Save solution to file
#file = File("poisson.pvd")
#file << u

# Plot solution
plot(u, interactive=True)
