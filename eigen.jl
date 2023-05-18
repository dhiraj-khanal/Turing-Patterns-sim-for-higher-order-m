using Plots
using LinearAlgebra
using Roots

# determinant of the Brusselator model
function model_determinant(k, d11, d22, a, b, u0, v0, m)
    return (-d11*k^2*u0^m+b-1)*(-d22*k^2*v0^m-a^2) + a^2*b
end

# wave numbers range
k = range(0.3, 1, length = 100)

plot()
m = 1
# find roots of the determinant equation when b=7.5
b_value = 13
a = 4.5
f(k) = model_determinant(k, 1, 10, a, b_value, a, b_value/a, m)
roots = find_zeros(f, 0.6, 2.1)

#b = 7.5
plot!(k, model_determinant.(k, 1, 10, a, b_value, a, b_value/a, m), label="b=13")
b = 7.5

plot!(k, model_determinant.(k, 1, 10, a, b, a, b/a, m), label="b=7.5")
hline!([0], color=:black)
title!("Determinant of matrix")
xlabel!("wave number k")
ylabel!("Determinant")
legend=:topleft

# shade the region where determinant is negative (spatial pattern occurs)
plot!(x->model_determinant(x, 1, 10, a, b_value, a, b_value/a, m), fillrange=0, fillalpha=0.2, fillcolor=:green, label="Pattern region for b=13", xlims=(0.3, 1))

# display roots on the plot
scatter!(roots, [0, 0], color=:red, markersize=5, label="Roots for b=13")

savefig("det_extended.png")
