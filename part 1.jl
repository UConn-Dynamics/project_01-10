### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 0bb123b9-041b-4e35-a6ad-2cc3598d509a
begin
	using Pkg
	Pkg.add.(["Symbolics", "DifferentialEquations", "Plots", "Latexify"])
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

# ╔═╡ d8dd769d-dad9-4236-9ece-a34e972f3f33
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