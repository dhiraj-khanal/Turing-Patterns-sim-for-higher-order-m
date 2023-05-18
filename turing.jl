using Plots

function del2_9(a::Matrix{Float64}, dx=1.0)
    #=
    Returns a finite-difference approximation of the laplacian of the array a,
    with lattice spacing dx, using the nine-point stencil:

                    1/6 2/3 1/6
                    2/3 -10/3 2/3
                    1/6 2/3 1/6
    =#

    del2 = zeros(size(a))
    del2[2:end-1, 2:end-1] .= (4.0 * (a[2:end-1, 3:end] + a[2:end-1, 1:end-2] +
                                      a[3:end, 2:end-1] + a[1:end-2, 2:end-1]) +
                               (a[1:end-2, 1:end-2] + a[1:end-2, 3:end] +
                                a[3:end, 1:end-2] + a[3:end, 3:end]) -
                               20.0 * a[2:end-1, 2:end-1]) ./ (6.0 * dx^2)
    return del2
end
function apply_neumann_conditions!(U::Matrix{Float64}, V::Matrix{Float64})
    for Z in (U, V)
        Z[1, :] = Z[2, :]
        Z[end, :] = Z[end-1, :]
        Z[:, 1] = Z[:, 2]
        Z[:, end] = Z[:, end-1]
    end
end

function simulate_reaction_diffusion(D11, D22, a, b, U, V, dt, dx, steps, fun=1)
    for t in 1:steps
        Uxx = del2_9(U, dx)[2:end-1, 2:end-1]
        Vxx = del2_9(V, dx)[2:end-1, 2:end-1]

        u = U[2:end-1, 2:end-1]
        v = V[2:end-1, 2:end-1]
        if fun == 0
            U[2:end-1, 2:end-1] = u + dt * (D11 * Uxx)
            V[2:end-1, 2:end-1] = v + dt * (D22 * Vxx)
        else
            U[2:end-1, 2:end-1] = u + dt * (D11 * Uxx + (b - 1) * u + a^2 * v)
            V[2:end-1, 2:end-1] = v + dt * (D22 * Vxx - b * u - a^2 * v)

        end

        apply_neumann_conditions!(U, V)
    end
    return U, V
end

concentrations_over_time = []
function verify_simulation(U0::Matrix{Float64}, V0::Matrix{Float64}, U::Matrix{Float64}, V::Matrix{Float64}, tol=1e-3)
    init_concentration = (sum(U0) + sum(V0))/n
    curr_concentration = (sum(U) + sum(V))/n

    append!(concentrations_over_time, curr_concentration)
    return isapprox(init_concentration, curr_concentration, atol=tol)
end

# parameters
D11 = 1.0
D22 = 10.0
a = 4.5
b = 7.5

# grid size 2D domain
n = 100
dx = 1.0

# Time
dt = 0.001
total_time = [0, 0.1, 0.5, 1, 5, 10] #six time periods

steps_per_plot = 1000*total_time

# initialize matrices U and V with random values between 0 and 1
U0 = rand(Float64, n, n)
V0 = rand(Float64, n, n)
#append!(concentrations_over_time, sum(U0) + sum(V0)/n)
# initialize a plots for u and v separately
plt = plot(layout=(2, 3), legend=false, grid=false, xlabel="x", ylabel="y", title=["t=0" "t=0.1" "t=0.5" "t=1" "t=5" "t=10"])
plt1 = plot(layout=(2, 3), legend=false, grid=false, xlabel="x", ylabel="y", title=["t=0" "t=0.1" "t=0.5" "t=1" "t=5" "t=10"])

# simulate the reaction-diffusion system and plot for specified time steps
for (i, time) in enumerate(total_time)
    steps = steps_per_plot[i]
    U, V = simulate_reaction_diffusion(D11, D22, a, b, U0, V0, dt, dx, steps, 1)

    # set f = 0
    #U, V = simulate_reaction_diffusion(D11, D22, a, b, U0, V0, dt, dx, steps, 0)
    # verify the simulation, prints false if initial conditions do not match the final conditions
    is_verified = verify_simulation(U0, V0, U, V)
    if !is_verified
        println("Simulation not verified at t=$time")
    end

    # heatmap of U, V to the corresponding subplots
    heatmap!(plt[i], U)
    heatmap!(plt1[i], V)
end
concentration_plot = plot(total_time, concentrations_over_time, xlabel="Time", ylabel="Concentration (U + V)", title="Concentration over time", legend=false, marker=:circle)
display(concentration_plot)
savefig(concentration_plot, "concentration_over_time.png")

display(plt)
savefig(plt, "U_resultss_with_reaction.png")
display(plt1)
savefig(plt1, "V_resultss_with_reaction.png")

