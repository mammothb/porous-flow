clear;
clf;
% Model parameters (fixed)
porosity = 0.1;  % Porosity of porous medium
Re = 0.01;  % Reynolds number
Da = 1e-5;  % Darcy number
J = 1.0;  % Viscosity ratio
tau = 0.8;  % Relaxation time
dx = 1.0;  % Lattice spacing
dt = 1.0;  % Time step
cs_sqr = dx * dx / dt / dt / 3.0;
nu = (tau - 0.5) * cs_sqr * dt;  % Kinematic viscosity
nu_e = J * nu;  % Effective parameter
H = 80;  % Height of domain

% Derived parameters
K = Da * H * H;  % Permeability
F_epsilon = 1.75 / sqrt(150 * porosity^3);  % Geometric function
u0 = Re * nu / H;  % Peak velocity
r = sqrt(nu * porosity / K / nu_e);
% Error in paper?
G = u0 * nu / K / (1 - 1 / cosh(r * H / 2));

% Simulation parameters
h = 0.08;  % Space step
n = H / h + 1;  % Number of elements
hh = h * h / 3;  % Convenience variable
g0 = 0.0;  % BC at y = 0
gH = 0.0;  % BC at y = H
A = h * h * F_epsilon * porosity / nu_e / sqrt(K);  % Coefficient for u^2
B = h * h * nu * porosity / K / nu_e;  % Coefficient for u
C = h * h * porosity / nu_e;  % Source term

y = linspace(0, H, n);  % Discretized domain
u = zeros(1, n);  % Initialize u
u_new = u;  % Initialize another u
iterations = 0;  % Number of iterations
tolerance = 1e-12;  % Convergence tolerance
is_unsteady = true;
u(1) = g0;  % Apply boundary conditions
u(end) = gH;  % Apply boundary conditions

% Start looping
while is_unsteady
  iterations = iterations + 1;
  u_new(2 : end - 1) = (u(3 : end) + u(2 : end - 1) + u(1 : end - 2) - A .*...
      u(2 : end - 1) .^ 2 - B .* u(2 : end - 1) + C * G) / 3.0;
  relative_err = max(abs(u_new - u));
  % fprintf('%d\n', relative_err);
  if relative_err < tolerance, is_unsteady = false; end
  u = u_new;
end
fprintf('Iterations: %d\n', iterations);
plot (y, u ./ u0);
