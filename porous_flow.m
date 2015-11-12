clear;
clf;
% Model parameters
H = 80;  % Height of domain
porosity = 0.1;  % Porosity of porous medium
Re = 0.01;  % Reynolds number
Da = 1e-5;  % Darcy number
J = 1.0;  % Viscosity ratio
tau = 0.8;  % Relaxation time
% dx = 0.0316;  % Lattice spacing
% dt = 0.001;  % Time step
dx = 1.0;  % Lattice spacing
dt = 1.0;  % Time step
nu = (0.8 - 0.5) * dt * dx * dx / dt / dt / 3.0;  % Kinematic viscosity
nu_e = J * nu;  % Effective parameter
G = 10.0;  % External force
K = Da * H * H;  % Permeability
F_epsilon = 1.75 / sqrt(150 * porosity^3);  % Geometric function
r = sqrt(nu * porosity / K / nu_e);
u0 = G * K / nu * (1 - acosh(r * H / 2));
% u0 = Re * nu / H / dx;
% disp(u0);
disp(acosh(r * H / 2));

% Simulation parameters
h = 0.08;  % Space step
n = H / h + 1;  % Number of elements
hh = h * h / 3;  % Convenience variable
g0 = 0.0;  % BC at y = 0
gH = 0.0;  % BC at y = H
A = h * h * F_epsilon * porosity / nu_e / sqrt(K);  % Coefficient for u^2
B = h * h * nu * porosity / K / nu_e;  % Coefficient for u
C = h * h * G * porosity / nu_e;  % Source term

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
      u(2 : end - 1) .^ 2 - B .* u(2 : end - 1) + C) / 3.0;
  relative_err = max(abs(u_new - u));
  % fprintf('%d\n', relative_err);
  if relative_err < tolerance, is_unsteady = false; end
  u = u_new;
end
fprintf('Iterations: %d\n', iterations);
plot (y, u);
