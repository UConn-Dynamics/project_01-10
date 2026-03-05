### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 57b99493-d454-425f-9398-e18283e27acb
using Markdown

# ╔═╡ e0b28f7c-fc64-46aa-b9fb-693cc505e1e0
using InteractiveUtils

# ╔═╡ f6eb02ee-544d-4169-8e36-6bc5140345d3
begin
	using Symbolics   # for symbolic math
	using ModelingToolkit
	using Latexify
end

# ╔═╡ 27e1db88-4970-4b61-a5dd-3b791523a89d
begin
	using DifferentialEquations   # for solving ODEs
	
	# Define parameter values
	length_L_val = 0.15
	mass_m_val = 0.1   # not needed for the final ODE, but handy if expanding later
	gravity_g_val = 9.81
	omega_val = 2.0    # example angular speed, rad/s
	
	# ODE system: u[1] = θ, u[2] = θ̇
	function pendulum_rotating!(du, u, params, time)
	    θ = u[1]
	    angular_velocity_theta = u[2]
	
	    L, g, ω = params
	
	    du[1] = angular_velocity_theta
	    du[2] = - (g / L) * sin(θ) +
	             (ω^2) * sin(θ) * cos(θ)
	end
end

# ╔═╡ e61429ae-2a7e-4022-9f49-419e9aee5f5f
begin
	using Plots   # for visualization
	
	# Initial conditions: small angle, zero initial angular velocity
	theta0 = 0.1        # radians
	theta_dot0 = 0.0
	initial_state = [theta0, theta_dot0]
	
	parameters = (length_L_val, gravity_g_val, omega_val)
	time_span = (0.0, 10.0)  # seconds
	
	problem = ODEProblem(pendulum_rotating!, initial_state, time_span, parameters)
	solution = solve(problem)
	
	# Plot θ(t)
	plot(solution.t, getindex.(solution.u, 1),
	     xlabel = "time (s)",
	     ylabel = "θ(t) (rad)",
	     title = "Pendulum in Rotating Frame: Angle vs Time",
	     legend = false)
end

# ╔═╡ f3892d09-568f-4ee2-8b7a-4c527a318bc7
using PlutoUI

# ╔═╡ 0bb123b9-041b-4e35-a6ad-2cc3598d509a
# ╠═╡ disabled = true
#=╠═╡
begin
	using Pkg
	Pkg.add.(["Symbolics", "DifferentialEquations", "Plots", "Latexify"])
end
  ╠═╡ =#

# ╔═╡ d8dd769d-dad9-4236-9ece-a34e972f3f33
#=╠═╡
begin
	@variables t
	# @variables angle_theta(time)  # θ(t)
	@syms θ(t)
	@parameters L m g ω
	
	# Define x, y, z of the mass in lab frame
	position_x = L * sin(θ(t)) * cos(ω * t)
	position_y = L * sin(θ(t)) * sin(ω * t)
	position_z = -L * cos(θ(t))
	
	# Compute velocity components
	velocity_x = Differential(t)(position_x)
	velocity_y = Differential(t)(position_y)
	velocity_z = Differential(t)(position_z)

	# Kinetic energy T = 1/2 m v^2
	velocity_squared = simplify(velocity_x^2 + velocity_y^2 + velocity_z^2)
	kinetic_T = (m / 2) * velocity_squared
	
	# Potential energy V = m g z
	potential_V = m * g * position_z
	
	# Lagrangian L = T - V
	lagrangian_L = simplify(kinetic_T - potential_V)

	@syms derivative(var1, var2)
	
	D = Differential(t)
	
	# Partial derivatives
	dL_dthetadot = expand_derivatives(derivative(lagrangian_L, D(θ(t))))
	dL_dtheta    = derivative(lagrangian_L, θ(t))
	
	# Euler–Lagrange: d/dt(∂L/∂θ̇) - ∂L/∂θ = 0
	eom = simplify(D(dL_dthetadot) - dL_dtheta)

	# length_L^2 * mass_m * D(D(angle_theta))
	# + mass_m * gravity_g * length_L * sin(angle_theta)
	# - mass_m * length_L^2 * omega^2 * sin(angle_theta) * cos(angle_theta) = 0
end
  ╠═╡ =#

# ╔═╡ 5ddcff18-a8dc-470a-8327-e612ac44457b
# begin
# 	# Kinetic energy T = 1/2 m v^2
# 	velocity_squared = simplify(velocity_x^2 + velocity_y^2 + velocity_z^2)
# 	kinetic_T = (mass_m / 2) * velocity_squared
# 	
# 	# Potential energy V = m g z
# 	potential_V = mass_m * gravity_g * position_z
# 	
# 	# Lagrangian L = T - V
# 	lagrangian_L = simplify(kinetic_T - potential_V)
# 
# 	@syms derivative(var1, var2)
# 	
# 	D = Differential(time)
# 	
# 	# Partial derivatives
# 	dL_dthetadot = expand_derivatives(derivative(lagrangian_L, D(angle_theta(time))))
# 	dL_dtheta    = derivative(lagrangian_L, angle_theta(time))
# 	
# 	# Euler–Lagrange: d/dt(∂L/∂θ̇) - ∂L/∂θ = 0
# 	eom = simplify(D(dL_dthetadot) - dL_dtheta)
# end

# ╔═╡ d4f931f6-47de-4ac7-8d1a-984f2ef4e222
# begin
# 	@syms derivative(var1, var2)
# 	
# 	D = Differential(time)
# 	
# 	# Partial derivatives
# 	dL_dthetadot = expand_derivatives(derivative(lagrangian_L, D(angle_theta(time))))
# 	dL_dtheta    = derivative(lagrangian_L, angle_theta(time))
# 	
# 	# Euler–Lagrange: d/dt(∂L/∂θ̇) - ∂L/∂θ = 0
# 	eom = simplify(D(dL_dthetadot) - dL_dtheta)
# end

# ╔═╡ 706014d3-b42c-42ac-bbec-67eebb84cdbb
# begin
# 	length_L^2 * mass_m * D(D(angle_theta))
# 	+ mass_m * gravity_g * length_L * sin(angle_theta)
# 	- mass_m * length_L^2 * omega^2 * sin(angle_theta) * cos(angle_theta) = 0
# end

# ╔═╡ 67c9f4c9-3c9e-4cd5-b3c4-8b9cea153a6e
begin
	ts = solution.t
	thetas = getindex.(solution.u, 1)
	thetadots = getindex.(solution.u, 2)
	nothing
end

# ╔═╡ 1b0e6eb1-90ac-461c-84ad-9c777b63c867
plot(thetas, thetadots,
	 xlabel="0 (rad)", ylabel="0 (rad/s)",
	 title="Phase Plot: 0 vs 0", legend=false)

# ╔═╡ ec5484c7-2275-413e-a112-42af1f7bfdf3
begin
	L = length_L_val
	ω = omega_val

	xs = L .* sin.(thetas) .* cos.(ω .* ts)
	ys = L .* sin.(thetas) .* sin.(ω .* ts)
	zs = -L .* cos.(thetas)
end

# ╔═╡ 3ae20080-f861-4222-8bfa-ecb4b554c608
plot(xs, ys, zs,
	xlabel="x", ylabel="y", zlabel="z",
	 title="Bob Trajectory (Lab Frame)",
	 legend=false, linewidth=2)

# ╔═╡ bd39aef5-5cf6-4ecd-b442-4f782228aaa7
@bind k Slider(1:length(ts), default=1, show_value=true)

# ╔═╡ 38bae75e-055d-4010-be7f-9ea2dbf8e330
begin
	xk, yk, zk = xs[k], ys[k], zs[k]
	trail = 1:k
	Lval = length_L_val

	plt = plot(xlim=(-Lval, Lval), ylim=(-Lval, Lval), zlim=(-Lval,  0.05),
		xlabel="x", ylabel="y", zlabel="z",
		title="Rotating Pendulum (interactive)",
		legend=false)

	plot!(xs[trail], ys[trail], zs[trail], linewidth=2)
	plot!([0, xk], [0, yk], [0, zk], linewidth=4)
	scatter!([xk], [yk], [zk], markersize=8)

	plt
end

# ╔═╡ Cell order:
# ╠═57b99493-d454-425f-9398-e18283e27acb
# ╠═e0b28f7c-fc64-46aa-b9fb-693cc505e1e0
# ╠═0bb123b9-041b-4e35-a6ad-2cc3598d509a
# ╠═f6eb02ee-544d-4169-8e36-6bc5140345d3
# ╠═d8dd769d-dad9-4236-9ece-a34e972f3f33
# ╟─5ddcff18-a8dc-470a-8327-e612ac44457b
# ╟─d4f931f6-47de-4ac7-8d1a-984f2ef4e222
# ╟─706014d3-b42c-42ac-bbec-67eebb84cdbb
# ╠═27e1db88-4970-4b61-a5dd-3b791523a89d
# ╠═e61429ae-2a7e-4022-9f49-419e9aee5f5f
# ╠═67c9f4c9-3c9e-4cd5-b3c4-8b9cea153a6e
# ╠═1b0e6eb1-90ac-461c-84ad-9c777b63c867
# ╠═ec5484c7-2275-413e-a112-42af1f7bfdf3
# ╠═3ae20080-f861-4222-8bfa-ecb4b554c608
# ╠═f3892d09-568f-4ee2-8b7a-4c527a318bc7
# ╠═bd39aef5-5cf6-4ecd-b442-4f782228aaa7
# ╠═38bae75e-055d-4010-be7f-9ea2dbf8e330
