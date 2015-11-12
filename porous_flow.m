clear;
clf;
L = 4;  % Length of domain
h = 0.5;  % Space step
n = L / h + 1;  % Number of elements
hh = h * h / 3;  % Convenience variable
u0 = 0.0;  % BC at x = 0
uL = 0.0;  % BC at x = L
ub = 0.5;  % Temperature of the background
ub4 = ub * ub * ub * ub;  % Convenience variable
x = linspace(0, L, n);  % Discretized domain
g = sin(pi * x / 4);  % Source term
u = zeros(1, n);  % Initialize u
u_new = u;  % Initialize another u
iterations = 0;  % Number of iterations
tolerance = 1e-6;  % Convergence tolerance
is_unsteady = true;
u(1) = u0;  % Apply boundary conditions
u(end) = uL;  % Apply boundary conditions

% Start looping
while is_unsteady
  iterations = iterations + 1;
  u_new(2 : end - 1) = (u(3 : end) + u(2 : end - 1) + u(1 : end - 2)) / 3 +...
      hh * (-u(2 : end - 1) .^ 4 + ub4 + g(2 : end - 1));
  relative_err = max(abs((u_new(2 : end - 1) - u(2 : end - 1)) ./...
      u_new(2 : end - 1)));
  fprintf('%d\n', relative_err);
  if relative_err < tolerance, is_unsteady = false; end
  u = u_new;
end
fprintf('Iterations: %d\n', iterations);
plot (x, u);
