from dolfin import *
from math import atan, pi

def A_B_sweep(A, B):
    # some parameters
    N        = 1000 # number grid cells
    theta    = Constant(1.0)
    n        = Constant(3.)
    m        = Constant(1.)
    Kb       = Constant(1.0)
    P0       = Constant(0.5)
    C_1       = Constant(0.5)
    C_2       = Constant(0.5)
    C_3       = Constant(-1.0)

    # define mesh
    mesh = UnitIntervalMesh(N)

    # define function space
    Q = FunctionSpace(mesh, "CG", 1) # pressure space
    W = FunctionSpace(mesh, "CG", 1) # porosity space
    C = FunctionSpace(mesh, "CG", 1) # concentration space
    ME = MixedFunctionSpace([W,Q,C])

    # define boundaries
    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] > 1.0 - DOLFIN_EPS and on_boundary
    
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] < DOLFIN_EPS and on_boundary

    left = Left()
    right = Right()

    # mark the subdomains
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    left.mark(sub_domains, 0)
    right.mark(sub_domains, 1)

    # dirichlet bcs
    right_p = Constant(0.)
    left_c = Constant(0.)
    right_c = Constant(0.)
    right_phi = Constant(0.001)
    left_phi = Constant(1.0)
    left_p = Constant(0.0)

    bc0 = DirichletBC(ME.sub(0), right_phi, right)
    bc4 = DirichletBC(ME.sub(0), left_phi, left)
    bc5 = DirichletBC(ME.sub(1), left_p, left)
    bc1 = DirichletBC(ME.sub(1), right_p, right)
    bc2 = DirichletBC(ME.sub(2), right_c, right)
    bc3 = DirichletBC(ME.sub(2), left_c, left)

    bcs = [bc0, bc1, bc2, bc3]

    # initial conditions for p & c
    p_0 = Constant(0.)
    c_0 = Constant(0.)

    p_initial = p_0
    c_initial = c_0

    p_prev = project(p_initial, Q)
    c_prev = project(c_initial, C)

    # initial condition for phi
    phi_min = 0.001 
    phi_bnd = 1.  
    xlambda = 0.01
    phi_initial = Expression("phi_min + (phi_bnd - phi_min)*exp(-x[0]/xlambda)",
                             phi_min=phi_min, phi_bnd=phi_bnd, xlambda=xlambda)

    phi_prev = project(phi_initial, W)

    # define trial and test functions
    dv = TrialFunction(ME)
    (phi_t, p_t, c_t) = TestFunctions(ME)
    vfunc = Function(ME)
    phi, p, c = split(vfunc)

    # time-step, initial and final times
    dt = 0.001       # time-step
    t  = 0          # initial time
    T  = 1000*dt      # final time

    # export
    file1 = File("porosity.pvd")
    #file2 = File("pressure.pvd")
    #file3 = File("conc.pvd")

    j = 0
    # compute solution
    while (t < T):

        Gamma      = phi*(c - 1.)
        Gamma_prev = phi_prev*(c_prev - 1.)
        pdiff      = phi**n*grad(p)
        pdiff_prev = phi_prev**n*grad(p_prev)
        fdiff      = phi**m*grad(phi)
        fdiff_prev = phi_prev**m*grad(phi_prev)

        phitheta  = theta*phi + (1. - theta)*phi_prev
        ptheta    = theta*p + (1. - theta)*p_prev
        ctheta    = theta*c + (1. - theta)*c_prev

        Gamma_theta = theta*Gamma + (1. - theta)*Gamma_prev
        pdiff_theta = theta*pdiff + (1. - theta)*pdiff_prev
        fdiff_theta = theta*fdiff + (1. - theta)*fdiff_prev

        s_p = p_t*Kb*(P0 - ptheta)

        r_phi = phi_t*(phi - phi_prev)*dx - phi_t*(p - p_prev)*dx - dt*phi_t*C_1*Gamma_theta*dx
        
        r_p = p_t*((p - p_prev)/A)*dx + dt*(inner(grad(p_t),pdiff_theta)*dx + (B/A)*inner(grad(p_t),fdiff_theta)*dx - (1./A)*p_t*C_2*Gamma_theta*dx) - (1./A)*s_p*ds(0)
        
        r_c = c_t*(c - c_prev)*dx - dt*c_t*C_3*Gamma_theta*dx

        r = r_phi + r_p + r_c

        # compute jacobian
        J = derivative(r, vfunc, dv)

        # Solve
        solve(r == 0, vfunc, bcs=bcs, J=J, solver_parameters={"newton_solver": {'absolute_tolerance': 1E-8, 'relative_tolerance': 1E-8, 'maximum_iterations': 25, 'relaxation_parameter': 1.0}})

        phi, p, c = split(vfunc)
        
        if j % 50 == 0:
            file1 << vfunc.sub(0)
     #   file2 << vfunc.sub(1)
     #   file3 << vfunc.sub(2)

        assign(phi_prev, vfunc.sub(0))
        assign(p_prev, vfunc.sub(1))
        assign(c_prev, vfunc.sub(2))

        t += dt
        j += 1

if __name__ == "__main__":
    A, B = Constant(100.), Constant(1.)
    A_B_sweep(A, B)
    
