using DifferentialEquations, Plots

# Constants
L = 10.0                # Length of the rod
T = 2.0                 # Total time
Nx = 101                # Number of spatial points
Nt = 1001               # Number of time points
α = 0.1                 # Thermal diffusivity

# Discretization
x = range(0, L, length=Nx)
Δx = x[2] - x[1]
t = range(0, T, length=Nt)
Δt = t[2] - t[1]

# Initial condition (Gaussian distribution)
u0 = exp.(-1.0 .* (x .- L/2).^2)

# Boundary conditions (constant temperature of 0 at both ends)
bcs = Dirichlet0BC(Float64)

# PDE definition
heat_eq = HeatEquation1D(α)

# Discretization using the central difference method
discretization = CenteredDifference(2, 2, Δx, Nx)

# PDE problem
pde_problem = PDEProblem(heat_eq, u0, (0.0, T), bcs, discretization)

# Solve the PDE
sol = solve(pde_problem, Tsit5(), saveat=Δt)

# Plot the solution
heatmap(sol, color=:viridis, xlabel="x", ylabel="t", title="Heat Equation Solution")

savefig("porous.png")

