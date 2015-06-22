from dolfin import *
from math import atan, pi

# some parameters
N        = 1000 # number grid cells
phi_max  = 0.6
step_size = 0.2
theta    = Constant(1.0)
A        = Constant(1000.)
B        = Constant(50.)
n        = Constant(3.)

# define mesh
mesh = UnitIntervalMesh(N)

# define function space
Q = FunctionSpace(mesh, "CG", 1) # pressure space
W = FunctionSpace(mesh, "CG", 1) # porosity space
ME = MixedFunctionSpace([W,Q])

# define boundaries
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > (1.0 - DOLFIN_EPS) and on_boundary

right_phi = Constant(phi_max - step_size)
left_phi  = Constant(phi_max)
p_0 = Constant(0.1)

left  = Left()
right = Right()

# dirichlet BCs
bc0 = DirichletBC(ME.sub(0), right_phi, right)   # specified porosity at right
bc1 = DirichletBC(ME.sub(0), left_phi, left)     # specified porosity at left
bc2 = DirichletBC(ME.sub(1), p_0, right)         # constant pressure at right
bc3 = DirichletBC(ME.sub(1), p_0, left)          # constant pressure at left

bcs = [bc0, bc1, bc2, bc3]

# initial conditions
k = Constant(1.e6) # k should be big

phi_initial = Expression("(-step_size/(pi))*atan(k*(x[0] - 0.5)) + (phi_max - (0.5*step_size))",
                         step_size=step_size, phi_max=phi_max, k=k)

p_initial = p_0

p_prev   = project(p_initial, Q)
phi_prev = project(phi_initial, W)


# define trial and test functions
dv = TrialFunction(ME)
(phi_t, p_t) = TestFunctions(ME)
vfunc = Function(ME)
phi, p = split(vfunc)

file_initial = File("initial_condition.pvd")

file_initial << vfunc.sub(0)

# time-step, initial and final times
dt = 0.0001       # time-step
t  = 0          # initial time
T  = 100*dt      # final time

# export
file = File("1d_dirichlet.pvd")

# compute solution
while (t < T):

    phitheta  = theta*phi + (1. - theta)*phi_prev
    ptheta    = theta*p + (1. - theta)*p_prev

    # weak forms
    r_phi = phi_t*(phi - phi_prev)*dx - phi_t*(1.)*(p - p_prev)*dx
    r_p = p_t*((p - p_prev)/A)*dx + dt*(inner(grad(p_t),(phitheta**n)*grad(ptheta))*dx + (B/A)*inner(grad(p_t),((phitheta**(n - 1.))*grad(phitheta)))*dx) 

    r = r_phi + r_p

    # compute jacobian
    J = derivative(r, vfunc, dv)

    # Solve
    solve(r == 0, vfunc, bcs, J=J, solver_parameters={"newton_solver": {'absolute_tolerance': 1E-8, 'relative_tolerance': 1E-8, 'maximum_iterations': 25, 'relaxation_parameter': 1.0}})

    phi, p = split(vfunc)
    file << vfunc.sub(0)

    assign(phi_prev, vfunc.sub(0))
    assign(p_prev, vfunc.sub(1))

    t += dt
