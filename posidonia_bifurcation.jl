#!/usr/bin/env julia

# Posidonia oceanica 1D bifurcation code
#
# Clean single-file version:
#   - model definition;
#   - finite-difference operators with Neumann boundary conditions;
#   - continuation of the homogeneous branch in omega_d0;
#   - detection of the first bifurcation point;
#   - optional continuation of the branch emerging from that point;
#   - JSON export for MATLAB plotting.

using LinearAlgebra
using SparseArrays
using Statistics
using JSON3
using Accessors
using BifurcationKit

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

"""
    nondimensional_parameters(; kwargs...)

Convert the dimensional parameter set to the nondimensional parameters
used by the 1D model.
"""
function nondimensional_parameters(;
    tilde_wb      = 0.05,
    tilde_wd0     = 0.001,
    tilde_k       = 0.041,
    tilde_a       = 100.0,
    tilde_b       = 8.0525,
    tilde_d0      = 508.0,
    tilde_d1      = 6508.0,
    tilde_alpha   = 530000.0,
    tilde_beta    = 1.6842e11,
    L_dim         = 33000.0,
    omega_d0_start = 0.3077
)
    time_scale   = tilde_a^2 / tilde_b
    length_scale = sqrt(tilde_d1 * tilde_a / tilde_b)

    omega_b = tilde_wb * time_scale
    omega_d0 = omega_d0_start
    k = tilde_k * time_scale

    d0 = tilde_a * tilde_d0 / tilde_d1
    alpha = tilde_alpha / tilde_d1
    beta = tilde_b * tilde_beta / (tilde_a * tilde_d1^2)
    L = L_dim / length_scale

    return (
        omega_b = omega_b,
        omega_d0 = omega_d0,
        k = k,
        d0 = d0,
        alpha = alpha,
        beta = beta,
        L = L,
        time_scale = time_scale,
        length_scale = length_scale,
    )
end


function homogeneous_equilibrium(p)
    delta = (p.k - p.omega_b)^2 + 4.0 * (p.omega_b - p.omega_d0)

    if delta < 0.0
        return 0.0
    end

    return ((p.omega_b - p.k) + sqrt(delta)) / 2.0
end


# ------------------------------------------------------------
# Finite-difference operators
# ------------------------------------------------------------

function build_operators(N::Int, L::Float64)
    dx = 2.0 * L / (N - 1)

    main_d2 = fill(-2.0, N)
    off_d2 = fill(1.0, N - 1)
    D2 = spdiagm(0 => main_d2, 1 => off_d2, -1 => off_d2)

    # Homogeneous Neumann boundary conditions
    D2[1, 1] = -2.0
    D2[1, 2] =  2.0
    D2[N, N] = -2.0
    D2[N, N - 1] = 2.0
    D2 ./= dx^2

    main_d1 = zeros(N)
    upper_d1 = fill(1.0, N - 1)
    lower_d1 = fill(-1.0, N - 1)
    D1 = spdiagm(0 => main_d1, 1 => upper_d1, -1 => lower_d1)

    # Zero derivative at the boundary
    D1[1, :] .= 0.0
    D1[N, :] .= 0.0
    D1 ./= 2.0 * dx

    D4 = D2 * D2

    return D1, D2, D4, dx
end


# ------------------------------------------------------------
# Model residual and Jacobian
# ------------------------------------------------------------

function F_posidonia(n, p, D1, D2, D4)
    nx = D1 * n
    nxx = D2 * n
    n4x = D4 * n

    reaction = (p.omega_b - p.omega_d0) .* n .-
               (p.k - p.omega_b) .* n.^2 .-
               n.^3

    diffusion = (p.d0 .- p.alpha .* n) .* nxx
    gradient_term = nx.^2
    biharmonic_term = -p.beta .* n .* n4x

    return reaction .+ diffusion .+ gradient_term .+ biharmonic_term
end


function J_posidonia(n, p, D1, D2, D4)
    nx = D1 * n
    nxx = D2 * n
    n4x = D4 * n

    diag_reaction = (p.omega_b - p.omega_d0) .-
                    2.0 * (p.k - p.omega_b) .* n .-
                    3.0 .* n.^2

    J_reaction = spdiagm(0 => diag_reaction)

    J_diffusion_operator = spdiagm(0 => (p.d0 .- p.alpha .* n)) * D2
    J_diffusion_local = spdiagm(0 => (-p.alpha .* nxx))

    J_gradient = 2.0 .* spdiagm(0 => nx) * D1

    J_biharmonic_operator = spdiagm(0 => (-p.beta .* n)) * D4
    J_biharmonic_local = spdiagm(0 => (-p.beta .* n4x))

    return J_reaction +
           J_diffusion_operator +
           J_diffusion_local +
           J_gradient +
           J_biharmonic_operator +
           J_biharmonic_local
end


# ------------------------------------------------------------
# JSON export
# ------------------------------------------------------------

function _branch_point_to_dict(pt)
    out = Dict{String, Any}()

    for field in (:step, :param, :l2, :maxn, :minn, :amp, :meanu, :stable)
        if hasproperty(pt, field)
            out[String(field)] = getproperty(pt, field)
        end
    end

    return out
end


function save_branch_json(br, p; filename = "data/br_primary.json")
    branch = [_branch_point_to_dict(pt) for pt in br.branch]

    special_points = Any[]
    if hasproperty(br, :specialpoint) && br.specialpoint !== nothing
        for (i, sp) in enumerate(br.specialpoint)
            step = hasproperty(sp, :step) ? Int(sp.step) : nothing

            param_value = nothing
            if step !== nothing && 1 <= step <= length(br.branch)
                param_value = br.branch[step].param
            end

            sp_type = hasproperty(sp, :type) ? string(sp.type) : string(typeof(sp))

            push!(
                special_points,
                Dict(
                    "index" => i,
                    "type" => sp_type,
                    "step" => step,
                    "param" => param_value,
                ),
            )
        end
    end

    fixed_params = Dict{String, Any}()
    for key in (:omega_b, :omega_d0, :k, :d0, :alpha, :beta, :L)
        fixed_params[String(key)] = getproperty(p, key)
    end

    out = Dict(
        "fixed_params" => fixed_params,
        "continuation_parameter" => "omega_d0",
        "branch" => branch,
        "special_points" => special_points,
    )

    mkpath(dirname(filename))
    open(filename, "w") do io
        JSON3.write(io, out)
    end

    println("Saved branch: $filename")
end


function first_bifurcation_index(br)
    if !hasproperty(br, :specialpoint) || br.specialpoint === nothing
        return nothing
    end

    for (i, sp) in enumerate(br.specialpoint)
        if hasproperty(sp, :type) && sp.type == :bp
            return i
        end
    end

    return nothing
end


# ------------------------------------------------------------
# Main run
# ------------------------------------------------------------

function main()
    mkpath("data")

    # Grid and parameters
    N = 120
    p = nondimensional_parameters()

    D1, D2, D4, dx = build_operators(N, p.L)

    println("Domain: [-L,L], L = $(p.L)")
    println("N = $N, dx = $dx")
    println("Continuation parameter: omega_d0")

    n_star = homogeneous_equilibrium(p)
    u0 = fill(n_star, N)

    println("Initial homogeneous equilibrium n* = $n_star")

    F = (u, par) -> F_posidonia(u, par, D1, D2, D4)
    J = (u, par) -> J_posidonia(u, par, D1, D2, D4)

    prob = BifurcationProblem(
        F,
        u0,
        p,
        @optic _.omega_d0;
        J = J,
        record_from_solution = (u, par; kwargs...) -> (
            l2 = norm(u) / sqrt(length(u)),
            maxn = maximum(u),
            minn = minimum(u),
            amp = (maximum(u) - minimum(u)) / 2.0,
            meanu = mean(u),
        ),
    )

    opts = ContinuationPar(
        ds = 1e-3,
        dsmin = 1e-8,
        dsmax = 1e-1,
        p_min = 0.0,
        p_max = 700.0,
        max_steps = 15000,
        detect_bifurcation = 3,
        nev = 120,
        tol_stability = 1e-8,
        n_inversion = 8,
        newton_options = NewtonPar(
            tol = 1e-8,
            max_iterations = 400,
            verbose = false,
            linesearch = true,
        ),
        save_sol_every_step = 20,
    )

    # Primary homogeneous branch
    br = continuation(prob, PALC(), opts; verbosity = 1)
    save_branch_json(br, p; filename = "data/br_primary.json")

    # First secondary branch, if a bifurcation point is detected.
    #
    # This is the minimal way to continue one branch emerging from the first
    # bifurcation point. If more branches are needed, repeat the same command
    # using a different index in br.specialpoint.
    ind_bp = first_bifurcation_index(br)

    if ind_bp === nothing
        println("No bifurcation point detected on the primary branch.")
        return
    end

    println("First bifurcation point detected at specialpoint index = $ind_bp")
    println("Continuing one branch from this point...")

    br_secondary = continuation(br, ind_bp, opts; verbosity = 1)
    save_branch_json(br_secondary, p; filename = "data/br_from_first_bp.json")
end


main()
