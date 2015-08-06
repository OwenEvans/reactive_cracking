from dolfin import *
from math import atan, pi

def A_B_sweep(A, B):
    # some parameters
    N        = 2000 # number grid cells
    phi_max  = 0.6
    step_size = 0.2
    theta    = Constant(1.0)
    #A        = Constant(1000.)
    #B        = Constant(50.)
    n        = Constant(3.)

    # define mesh
    mesh = UnitIntervalMesh(N)

    # define function space
    Q = FunctionSpace(mesh, "CG", 1) # pressure space
    W = FunctionSpace(mesh, "CG", 1) # porosity space
    ME = MixedFunctionSpace([W,Q])

    p_0 = Constant(0.1)

    # initial conditions
    phi_min = 0.001
    phi_bnd = 1.
    xlambda = 0.01
    phi_initial = Expression("phi_min + (phi_bnd - phi_min)*exp(-x[0]/xlambda)",
                             phi_min=phi_min, phi_bnd=phi_bnd, xlambda=xlambda)

    phi_prev = project(phi_initial, W)

    #k = Constant(1.e6) # k should be big
    #phi_initial = Expression("(-step_size/(pi))*atan(k*(x[0] - 0.5)) + (phi_max - (0.5*step_size))",
    #                         step_size=step_size, phi_max=phi_max, k=k)
    
    p_initial = p_0
    p_prev   = project(p_initial, Q)
    phi_prev = project(phi_initial, W)

    # define trial and test functions
    dv = TrialFunction(ME)
    (phi_t, p_t) = TestFunctions(ME)
    vfunc = Function(ME)
    phi, p = split(vfunc)

    #file_initial = File("initial_condition_%s_%s.pvd" % (A, B))

    #file_initial << vfunc.sub(0)

    # time-step, initial and final times
    dt = 0.0005       # time-step
    t  = 0          # initial time
    T  = 200*dt      # final time

    # export
    file = File("porosity.pvd")

    # compute solution
    while (t < T):

        phitheta  = theta*phi + (1. - theta)*phi_prev
        ptheta    = theta*p + (1. - theta)*p_prev

        # weak forms
        r_phi = phi_t*(phi - phi_prev)*dx - phi_t*(1.)*(p - p_prev)*dx
        r_p = p_t*((p - p_prev)/A)*dx + dt*(inner(grad(p_t),(phitheta**n)*grad(ptheta))*dx + (B/A)*inner(grad(p_t),((phitheta**n)*(grad(phitheta))))*dx) 

        r = r_phi + r_p

        # compute jacobian
        J = derivative(r, vfunc, dv)

        # Solve
        solve(r == 0, vfunc, J=J, solver_parameters={"newton_solver": {'absolute_tolerance': 1E-8, 'relative_tolerance': 1E-8, 'maximum_iterations': 25, 'relaxation_parameter': 1.0}})

        phi, p = split(vfunc)
        file << vfunc.sub(0)

        assign(phi_prev, vfunc.sub(0))
        assign(p_prev, vfunc.sub(1))

        t += dt

if __name__ == "__main__":
    A, B = Constant(1.), Constant(10.)
    A_B_sweep(A, B)
    
