#!/usr/bin/env python3
"""
Clean FEniCS implementation of the 2D Posidonia oceanica model.

This script contains only the procedure used for the manuscript simulations:
    * nondimensionalisation from the selected dimensional parameter set;
    * square 2D computational domain;
    * two Gaussian vegetation patches as initial condition;
    * mixed finite-element formulation for the density n and its Laplacian;
    * time stepping with fixed time step;
    * XDMF output for post-processing.

The removed parts of the original research script include geographic polygons,
bathymetry-dependent parameters, stochastic mortality, adaptive time stepping,
reload/restart utilities, and plotting/debug blocks that were not active in the
simulation configuration used here.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path

from dolfin import *  # FEniCS/DOLFIN finite-element interface
from mshr import Polygon, generate_mesh


@dataclass
class SimulationConfig:
    """User-configurable numerical settings in dimensional units."""

    output_dir: Path = Path("data")
    mesh_resolution: int = 80
    final_time: float = 1000.0
    dt: float = 0.5
    save_every: float = 1.0
    domain_length: float = 1250.0


@dataclass
class ModelParameters:
    """Nondimensional model parameters and scaling factors."""

    omega_b: float
    omega_d0: float
    k: float
    a: float
    b: float
    d0: float
    d1: float
    alpha: float
    beta: float
    L: float
    Ly: float
    T: float
    dt: float
    save_every: float
    time_scale: float
    length_scale: float
    density_scale: float


def positive_equilibrium(a: float, omega: float, b: float) -> float:
    """Positive homogeneous equilibrium of the kinetic equation."""
    return (-a + math.sqrt(a * a + 4.0 * omega * b)) / (2.0 * b)


def nondimensional_parameters(config: SimulationConfig) -> ModelParameters:
    """
    Build the nondimensional parameter set used in the simulations.

    Dimensional starting values:
        omega_b  = 0.06  1/y
        omega_d0 = 0.072 1/y
        k        = 0.048 1/y
        a_tilde  = 100.41 cm^2
        b_tilde  = 12.5  cm^4/y
        d0_tilde = 508   cm^2/y
        d1_tilde = 6508  cm^4/y
        alpha_tilde = 2.0e4 cm^4/y
        beta_tilde  = 3.5e7 cm^6/y
    """

    # Dimensional kinetic parameters
    omega_b_dim = 0.06
    omega_d0_dim = 0.072
    k_dim = 0.048

    # Dimensional local-interaction and spatial parameters
    a_tilde = 100.41
    b_tilde = 12.5
    d0_tilde = 508.0
    d1_tilde = 6508.0
    alpha_tilde = 2.0e4
    beta_tilde = 3.5e7

    # Scaling used in the manuscript
    time_scale = a_tilde**2 / b_tilde
    length_scale = math.sqrt(d1_tilde * a_tilde / b_tilde)
    density_scale = 1.0 / a_tilde

    omega_b = omega_b_dim * time_scale
    omega_d0 = omega_d0_dim * time_scale
    k = k_dim * time_scale

    d0 = d0_tilde * a_tilde / d1_tilde
    alpha = alpha_tilde / d1_tilde
    beta = beta_tilde * b_tilde / (a_tilde * d1_tilde**2)

    # In the nondimensional form, a_tilde and b_tilde become 1.
    a = k - omega_b
    b = 1.0
    d1 = 1.0

    return ModelParameters(
        omega_b=omega_b,
        omega_d0=omega_d0,
        k=k,
        a=a,
        b=b,
        d0=d0,
        d1=d1,
        alpha=alpha,
        beta=beta,
        L=config.domain_length / length_scale,
        Ly=config.domain_length / length_scale,
        T=config.final_time / time_scale,
        dt=config.dt / time_scale,
        save_every=config.save_every / time_scale,
        time_scale=time_scale,
        length_scale=length_scale,
        density_scale=density_scale,
    )


class TwoGaussianInitialCondition(UserExpression):
    """Two localized vegetation patches centred symmetrically in the domain."""

    def __init__(self, params: ModelParameters, **kwargs):
        super().__init__(**kwargs)
        self.params = params
        omega = params.omega_b - params.omega_d0
        self.n_star = positive_equilibrium(params.a, omega, params.b)

    def eval(self, values, x):
        L = self.params.L
        sigma = L / 10.0
        c = L / 3.33

        patch_1 = math.exp(-((x[0] - c) ** 2 + (x[1] - c) ** 2) / sigma**2)
        patch_2 = math.exp(-((x[0] + c) ** 2 + (x[1] + c) ** 2) / sigma**2)

        values[0] = self.n_star * (patch_1 + patch_2)

    def value_shape(self):
        return ()


def create_square_mesh(params: ModelParameters, resolution: int):
    """Create the square domain used in the simulations."""
    half_side = 1.25 * params.L
    vertices = [
        Point(-half_side, -half_side),
        Point(half_side, -half_side),
        Point(half_side, half_side),
        Point(-half_side, half_side),
    ]
    return generate_mesh(Polygon(vertices), resolution)


def print_diagnostics(params: ModelParameters, rank: int) -> None:
    """Print the main derived quantities used to check the run."""
    if rank != 0:
        return

    omega = params.omega_b - params.omega_d0
    n_plus = positive_equilibrium(params.a, omega, params.b)
    n_minus = (-params.a - math.sqrt(params.a**2 + 4.0 * omega * params.b)) / (2.0 * params.b)

    kinetic_lambda = -params.a * n_plus - 2.0 * params.b * n_plus**2
    beta_critical = ((n_plus * params.alpha - params.d0) ** 2) / (-4.0 * n_plus * kinetic_lambda)
    q_critical = math.sqrt((-params.d0 + n_plus * params.alpha) / (2.0 * n_plus * params.beta))

    print("\n=== Derived nondimensional quantities ===")
    print(f"n_+                       = {n_plus:.8g}")
    print(f"n_-                       = {n_minus:.8g}")
    print(f"alpha_min = d0/n_+        = {params.d0 / n_plus:.8g}")
    print(f"beta_critical             = {beta_critical:.8g}")
    print(f"eta^2 = (beta_c-beta)/bc  = {(beta_critical - params.beta) / beta_critical:.8g}")
    print(f"kinetic lambda            = {kinetic_lambda:.8g}")
    print(f"q_critical                = {q_critical:.8g}")
    print("=========================================\n")


def save_parameter_file(params: ModelParameters, output_dir: Path, rank: int) -> None:
    """Save the numerical parameters used in the run."""
    if rank != 0:
        return

    omega = params.omega_b - params.omega_d0
    n_plus = positive_equilibrium(params.a, omega, params.b)
    kinetic_lambda = -params.a * n_plus - 2.0 * params.b * n_plus**2
    beta_critical = ((n_plus * params.alpha - params.d0) ** 2) / (-4.0 * n_plus * kinetic_lambda)
    q_critical = math.sqrt((-params.d0 + n_plus * params.alpha) / (2.0 * n_plus * params.beta))

    with open(output_dir / "parameters.txt", "w", encoding="utf-8") as f:
        f.write("Nondimensional Posidonia model parameters\n")
        f.write(f"omega_b       = {params.omega_b}\n")
        f.write(f"omega_d0      = {params.omega_d0}\n")
        f.write(f"k             = {params.k}\n")
        f.write(f"a             = {params.a}\n")
        f.write(f"b             = {params.b}\n")
        f.write(f"d0            = {params.d0}\n")
        f.write(f"d1            = {params.d1}\n")
        f.write(f"alpha         = {params.alpha}\n")
        f.write(f"beta          = {params.beta}\n")
        f.write(f"L             = {params.L}\n")
        f.write(f"Ly            = {params.Ly}\n")
        f.write(f"T             = {params.T}\n")
        f.write(f"dt            = {params.dt}\n")
        f.write(f"save_every    = {params.save_every}\n")
        f.write("\nDerived quantities\n")
        f.write(f"n_plus        = {n_plus}\n")
        f.write(f"beta_critical = {beta_critical}\n")
        f.write(f"q_critical    = {q_critical}\n")
        f.write(f"time_scale    = {params.time_scale}\n")
        f.write(f"length_scale  = {params.length_scale}\n")
        f.write(f"density_scale = {params.density_scale}\n")


def write_last_checkpoint(function, mesh, output_dir: Path, time_value: float, name: str) -> None:
    """Overwrite the last-step checkpoint in a clean way."""
    filename = output_dir / f"{name}_LAST.xdmf"
    with XDMFFile(mesh.mpi_comm(), str(filename)) as xdmf:
        xdmf.write_checkpoint(function, name, time_value, XDMFFile.Encoding.HDF5, True)


def run_simulation(config: SimulationConfig) -> None:
    """Run the finite-element simulation."""
    parameters["allow_extrapolation"] = True

    comm = MPI.comm_world
    rank = MPI.rank(comm)

    if rank == 0:
        config.output_dir.mkdir(parents=True, exist_ok=True)
    MPI.barrier(comm)

    params = nondimensional_parameters(config)
    print_diagnostics(params, rank)
    save_parameter_file(params, config.output_dir, rank)

    mesh = create_square_mesh(params, config.mesh_resolution)

    with XDMFFile(comm, str(config.output_dir / "domain.xdmf")) as xdmf:
        xdmf.write(mesh)

    # Mixed formulation: u = (n, laplacian(n)).
    element_n = FiniteElement("CG", triangle, 1)
    element_lap = FiniteElement("CG", triangle, 1)
    mixed_element = MixedElement([element_n, element_lap])

    V = FunctionSpace(mesh, mixed_element)
    Q = V.sub(0).collapse()

    v_n, v_lap = TestFunctions(V)

    u = Function(V)
    u_old = Function(V)

    n, lap_n = split(u)
    n_old, _ = split(u_old)

    # Initial condition for n. The Laplacian component starts from zero.
    initial_condition = TwoGaussianInitialCondition(params, degree=2)
    assign(u_old.sub(0), interpolate(initial_condition, Q))

    dt = Constant(params.dt)
    omega_b = Constant(params.omega_b)
    omega_d0 = Constant(params.omega_d0)
    a = Constant(params.a)
    b = Constant(params.b)
    d0 = Constant(params.d0)
    d1 = Constant(params.d1)
    alpha = Constant(params.alpha)
    beta = Constant(params.beta)

    # Weak formulation of the nondimensional model.
    residual = (
        ((n - n_old) / dt) * v_n * dx
        - ((omega_b - omega_d0) * n - a * n**2 - b * n**3) * v_n * dx
        + d0 * dot(grad(0.5 * (n + n_old)), grad(v_n)) * dx
        - alpha * dot(grad(n), grad(n)) * v_n * dx
        - alpha * n * dot(grad(n), grad(v_n)) * dx
        - d1 * dot(grad(n), grad(n)) * v_n * dx
        - beta * dot(grad(lap_n), grad(n)) * v_n * dx
        - beta * n * dot(grad(lap_n), grad(v_n)) * dx
        + dot(grad(n), grad(v_lap)) * dx
        + lap_n * v_lap * dx
    )

    xdmf_n = XDMFFile(comm, str(config.output_dir / "density.xdmf"))
    xdmf_lap = XDMFFile(comm, str(config.output_dir / "laplacian_density.xdmf"))
    for xdmf in (xdmf_n, xdmf_lap):
        xdmf.parameters["flush_output"] = True
        xdmf.parameters["functions_share_mesh"] = True

    time_log = None
    mortality_log = None
    if rank == 0:
        time_log = open(config.output_dir / "time_log.txt", "w", encoding="utf-8")
        mortality_log = open(config.output_dir / "mortality_log.txt", "w", encoding="utf-8")

    def save_state(time_value: float) -> None:
        density, lap_density = u_old.split()
        xdmf_n.write(density, time_value)
        xdmf_lap.write(lap_density, time_value)
        write_last_checkpoint(density, mesh, config.output_dir, time_value, "density")
        write_last_checkpoint(lap_density, mesh, config.output_dir, time_value, "laplacian_density")

        if rank == 0:
            time_log.write(f"{time_value}\n")
            mortality_log.write(f"{params.omega_d0}\n")
            time_log.flush()
            mortality_log.flush()

    save_state(0.0)

    t = 0.0
    last_save = 0.0
    step = 0
    wall_time = time.time()

    set_log_active(False)

    while t < params.T - 1.0e-14:
        solve(residual == 0, u)
        u_old.assign(u)

        t += params.dt
        step += 1

        if rank == 0:
            density_values = u_old.sub(0).compute_vertex_values(mesh)
            print(
                f"step={step:05d}  "
                f"t={t:.6g}/{params.T:.6g}  "
                f"max(n)={density_values.max():.6g}  "
                f"min(n)={density_values.min():.6g}  "
                f"elapsed={time.time() - wall_time:.2f}s"
            )
            wall_time = time.time()

        if t - last_save >= params.save_every - 1.0e-14:
            save_state(t)
            last_save = t

    if t > last_save:
        save_state(t)

    if rank == 0:
        time_log.close()
        mortality_log.close()

    xdmf_n.close()
    xdmf_lap.close()


def parse_args() -> SimulationConfig:
    parser = argparse.ArgumentParser(
        description="Run the cleaned 2D Posidonia oceanica FEniCS simulation."
    )
    parser.add_argument("--output-dir", default="data", type=Path)
    parser.add_argument("--mesh-resolution", default=80, type=int)
    parser.add_argument("--final-time", default=1000.0, type=float)
    parser.add_argument("--dt", default=0.5, type=float)
    parser.add_argument("--save-every", default=1.0, type=float)
    parser.add_argument("--domain-length", default=1250.0, type=float)

    args = parser.parse_args()

    return SimulationConfig(
        output_dir=args.output_dir,
        mesh_resolution=args.mesh_resolution,
        final_time=args.final_time,
        dt=args.dt,
        save_every=args.save_every,
        domain_length=args.domain_length,
    )


if __name__ == "__main__":
    run_simulation(parse_args())
