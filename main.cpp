#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

int main()
{
  /// Model parameters (fixed)
  auto porosity = 0.1;  // Porosity of porous medium
  auto reynold = 0.01;  // Reynolds number
  auto darcy = 1e-5;  // Darcy number
  auto visco_ratio = 1.0;  // Viscosity ratio
  auto tau = 0.8;  // Relaxation time
  auto dx = 0.1;  // Space step
  auto dt = 0.01;  // Time step
  auto cs_sqr = dx * dx / dt / dt / 3.0;
  auto k_visco = (tau - 0.5) * cs_sqr * dt;  // Kinematic viscosity
  auto k_visco_e = visco_ratio * k_visco;  // Effective parameter
  auto height = 80.0;  // Height of domain

  /// Derived parameters
  auto porosity_cube = porosity * porosity * porosity;
  auto permeability = darcy * height * height;
  auto solid_particle_diameter = std::sqrt(permeability * 150.0 *
      (1 - porosity) * (1 - porosity) / porosity_cube);
  auto geometric_fun = 1.75 / std::sqrt(150.0 * porosity_cube);
  auto u_max = reynold * k_visco / height;  // Peak velocity
  // r-term to simplify equation
  auto r = std::sqrt(k_visco * porosity / permeability / k_visco_e);
  auto body_force = u_max * k_visco / permeability / (1 - 1 /
      cosh(r * height / 2.0));

  std::cout << "Permeability: " << permeability << std::endl <<
      "Solid particle diameter: " << solid_particle_diameter << std::endl <<
      "Body force: " << body_force << std::endl;

  /// Simulation parameters
  auto num_nodes = 1001u;
  auto h = height / static_cast<double>(num_nodes - 1);
  auto bc_1 = 0.0;  // Boundary condition at y = 0
  auto bc_2 = 0.0;  // Boundary condition at y = H
  // Coefficients
  auto a = h * h * geometric_fun * porosity / k_visco_e /
      std::sqrt(permeability);
  auto b = h * h * k_visco * porosity / permeability / k_visco_e;
  auto c = h * h * porosity / k_visco_e;

  /// Preparation for iterating
  std::vector<double> u(num_nodes, 0.0);
  auto u_new = u;
  auto iterations = 0;
  auto tolerance = 1e-12;
  auto is_unsteady = true;
  // Apply boundary conditions
  u[0] = bc_1;
  u[num_nodes - 1] = bc_2;

  // Start iterating
  while (is_unsteady) {
    ++iterations;
    for (auto i = 1u; i < num_nodes - 1; ++i) {
      u_new[i] = (u[i + 1] + u[i] + u[i - 1] - a * u[i] * u[i] - b * u[i] + c *
          body_force) / 3.0;
    }  // i
    auto relative_err = 0.0;
    for (auto i = 0u; i < num_nodes; ++i) {
      auto abs_err = fabs(u_new[i] - u[i]);
      if (abs_err > relative_err) relative_err = abs_err;
    }  // i
//    std::cout << relative_err << std::endl;
    if (relative_err < tolerance) is_unsteady = false;
    u = u_new;
  }  // is_unsteady

  std::ofstream data_file;
  data_file.open("porous_flow.dat");
  auto y = 0.0;
  for (auto node : u) {
    data_file << y << " " << node / u_max << std::endl;
    y += h;
  }  // node


  return 0;
}
