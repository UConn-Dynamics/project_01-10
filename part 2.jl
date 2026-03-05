### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ aa4029bd-2d7d-4d30-8ec1-16a6f00e2f48
begin
	using Pkg
	Pkg.add.(["Symbolics", "DifferentialEquations", "Plots", "Latexify"])
end

# ╔═╡ 994e4b4b-74a2-4847-aec9-d33f0c8d97be
using Markdown

# ╔═╡ a268d2fa-e7b9-4c5e-84e1-f37b7c738e0a
using InteractiveUtils

# ╔═╡ c63781e1-529a-4357-a11b-8fe3ecc99057
begin
	using Symbolics   # for symbolic math
	using ModelingToolkit
	using Latexify
end

# ╔═╡ 3594dfba-0564-4ab2-9b08-e4d8bba8a247
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

# ╔═╡ 26f48500-2caf-4c8e-80ec-a372befccdac
begin
	using Plots   # for visualization
	
	# Initial conditions: small angle, zero initial angular velocity
	θ₀ = 0.1        # radians
	θ̇₀ = 0.0
	initial_state = [θ₀, θ̇₀]
	
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

# ╔═╡ bc78a3ce-f15d-4e08-94a0-82f68b3eb4b1
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

# ╔═╡ 0a3c8ae6-ec78-4ae7-8d61-918a4d766ff6
begin
	# Physical constants and initial conditions
	length_L = 0.15              # pendulum length (meters)
	gravity_g = 9.81             # gravity (m/s^2)
	θ₀2 = 0.15                    # initial angle (radians)
	θ̇₀2 = 0.0                     # initial angular velocity (rad/s)
	# time_span = (0.0, 10.0)      # solve for 10 seconds
	
	# Define parameter sets
	slow_omega = 0.5              # slow rotation speed
	fast_omega = 8.0              # fast rotation speed
end

# ╔═╡ f3bea490-12a4-11f1-a034-6fb93eb60805
# Equation of motion for pendulum in a rotating frame
# u[1] = θ(t), u[2] = θ̇(t)
function rotating_pendulum!(du, u, parameters, time)
    L, g, ω = parameters

    angle_theta = u[1]
    angular_velocity = u[2]

    du[1] = angular_velocity
    du[2] = - (gravity_g / length_L) * sin(angle_theta) +
             (ω^2) * sin(angle_theta) * cos(angle_theta)
end

# ╔═╡ fa9910bb-acb6-4c4d-a5c0-ee658f7a06e8
begin
	# ODE problems for slow and fast rotation
	problem_slow = ODEProblem(rotating_pendulum!, initial_state, time_span,
	                          (L, g, slow_omega))
	
	problem_fast = ODEProblem(rotating_pendulum!, initial_state, time_span,
	                          (L, g, fast_omega))
	
	solution_slow = solve(problem_slow)
	solution_fast = solve(problem_fast)
end

# ╔═╡ 2d7a4371-7dc0-48f9-90de-21e3402728f8
begin
	# Plot θ(t) for slow and fast rotation
	plot(solution_slow.t, getindex.(solution_slow.u, 1),
	     label="ω = 0.5 rad/s", xlabel="time (s)", ylabel="θ(t) (rad)",
	     title="Pendulum Motion: Slow vs Fast Rotation")
	
	plot!(solution_fast.t, getindex.(solution_fast.u, 1),
	      label="ω = 8.0 rad/s")
end

# ╔═╡ Cell order:
# ╠═994e4b4b-74a2-4847-aec9-d33f0c8d97be
# ╠═a268d2fa-e7b9-4c5e-84e1-f37b7c738e0a
# ╠═aa4029bd-2d7d-4d30-8ec1-16a6f00e2f48
# ╠═c63781e1-529a-4357-a11b-8fe3ecc99057
# ╠═bc78a3ce-f15d-4e08-94a0-82f68b3eb4b1
# ╠═3594dfba-0564-4ab2-9b08-e4d8bba8a247
# ╠═26f48500-2caf-4c8e-80ec-a372befccdac
# ╠═f3bea490-12a4-11f1-a034-6fb93eb60805
# ╠═0a3c8ae6-ec78-4ae7-8d61-918a4d766ff6
# ╠═fa9910bb-acb6-4c4d-a5c0-ee658f7a06e8
# ╠═2d7a4371-7dc0-48f9-90de-21e3402728f8