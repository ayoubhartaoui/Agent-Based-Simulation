# Packages to import. Check the righ version in Project.toml
using Random, StatsBase, Distributions
using DataFrames, CSV
using Printf

# ============================================================
# Configuration of the simulation
# ============================================================


struct Config
    G::Int
    T::Int
    Cw::Float64

    env::Symbol          # :low | :uncorrelated | :tradeoff | :richpoor
    trade::Symbol        # :minimal_law | :cobb_douglas
    avoidance::Symbol    # :normal | :constrained

    mu_b::Float64
    mu_l::Float64
    sigma_b::Float64
    sigma_l::Float64


    # modify this number to control the number of replciates
    # CAREFUL: the number of replicate has an high impact on computation time and output file size
    num_runs::Int

    # window time for series recoding  
    record_every::Int                   # timeseries sampling step
    eq_start::Int                       # equilibrium window start (for boxplot means)
    eq_end::Int                         # equilibrium window end   (for boxplot means)

    # window time for trait distributions (histograms/heatmaps)
    snapshot_start::Int
    snapshot_every::Int
    sample_size::Int

    outdir::String
end

mutable struct State
    b::Vector{Float64}
    λ::Vector{Float64}
    r1::Vector{Float64}
    r2::Vector{Float64}
    fitness::Vector{Float64}
    parents::Vector{Int}
end

# ============================================================
# Helpers
# ============================================================

# Minimal resources and clamping
min2(a,b) = a < b ? a : b
clamp01(x) = clamp(x, 0.0, 1.0)

function ensure_dir(path::String)
    isdir(path) || mkpath(path)
    return path
end

# ============================================================
# Environment draws
# ============================================================

struct EnvDists
    beta_rich::Beta
    beta_poor::Beta
    high_dist::Beta
    low_dist::Beta
    low_var_norm::Distribution
end

function make_env_dists()
    return EnvDists(
        Beta(15.0, 1.5),
        Beta(1.5, 15.0),
        Beta(5.0, 1.5),
        Beta(1.5, 5.0),
        truncated(Normal(0.5, 0.1), 0.0, 1.0)
    )
end

# ---- For the different environment types see the ReadMe or refer to the master thesis ----

function draw_environment!(st::State, cfg::Config, dists::EnvDists)
    G = cfg.G

    if cfg.env == :uncorrelated
        @inbounds for i in 1:G
            st.r1[i] = rand()
            st.r2[i] = rand()
        end

    elseif cfg.env == :low
        @inbounds for i in 1:G
            st.r1[i] = rand(dists.low_var_norm)
            st.r2[i] = rand(dists.low_var_norm)
        end

    elseif cfg.env == :richpoor
        n_rich = G ÷ 2
        is_rich = falses(G)
        is_rich[sample(1:G, n_rich; replace=false)] .= true
        @inbounds for i in 1:G
            if is_rich[i]
                st.r1[i] = rand(dists.beta_rich)
                st.r2[i] = rand(dists.beta_rich)
            else
                st.r1[i] = rand(dists.beta_poor)
                st.r2[i] = rand(dists.beta_poor)
            end
        end

    elseif cfg.env == :tradeoff
        @inbounds for i in 1:G
            if rand() < 0.5
                st.r1[i] = rand(dists.high_dist)
                st.r2[i] = rand(dists.low_dist)
            else
                st.r1[i] = rand(dists.low_dist)
                st.r2[i] = rand(dists.high_dist)
            end
        end

    else
        error("Unknown env = $(cfg.env)")
    end
end

# ============================================================
# Decisions & Interaction Resolution
# action codes: 0 = avoidance, 1 = Conflict, 2 = cooperative exchange
# ============================================================

# ----- Decision Rule based the limiting resource -----
@inline function choose_action(b::Float64, λ::Float64, r1::Float64, r2::Float64)
    thr = min2(r1, r2)  # Define minimal resource
    if λ >= thr
        return rand() < b ? 1 : 2 # Reminder: the belligernece parameter b is a probability to choose conflict
    else
        return 0
    end
end

# ---------------------
# Trade rules
# ---------------------

function trade_minimal_law!(st::State, j::Int, k::Int)
    R1_1, R2_1 = st.r1[j], st.r2[j]
    R1_2, R2_2 = st.r1[k], st.r2[k]

    surplus1 = max(R1_1, R2_1) - min(R1_1, R2_1)
    surplus2 = max(R1_2, R2_2) - min(R1_2, R2_2)

    receive1 = :r1
    receive2 = :r2
    t1 = 0.0
    t2 = 0.0

    if (R2_1 > R1_1) && (R1_2 > R2_2)
        receive1 = :r1
        receive2 = :r2
        t1 = min(R1_2, R2_2) / (min(R1_2, R2_2) + min(R1_1, R2_1))
        t2 = 1 - t1
    elseif (R1_1 > R2_1) && (R2_2 > R1_2)
        receive1 = :r2
        receive2 = :r1
        t1 = max(R1_2, R2_2) / (max(R1_2, R2_2) + max(R1_1, R2_1))
        t2 = 1 - t1
    end

    g1_gives = surplus1 * t1
    g2_gives = surplus2 * t2

    new_R1_1 = receive1 == :r1 ? (R1_1 + g2_gives) : (R1_1 - g1_gives)
    new_R2_1 = receive1 == :r2 ? (R2_1 + g2_gives) : (R2_1 - g1_gives)
    new_R1_2 = receive2 == :r1 ? (R1_2 + g1_gives) : (R1_2 - g2_gives)
    new_R2_2 = receive2 == :r2 ? (R2_2 + g1_gives) : (R2_2 - g2_gives) 

    
    min1_before = min2(R1_1, R2_1)
    min2_before = min2(R1_2, R2_2)
    min1_after  = min2(new_R1_1, new_R2_1)
    min2_after  = min2(new_R1_2, new_R2_2)

    if (min1_after <= min1_before) || (min2_after <= min2_before)
        return false
    end

    st.r1[j], st.r2[j] = new_R1_1, new_R2_1
    st.r1[k], st.r2[k] = new_R1_2, new_R2_2
    return true
end

function trade_cobb_douglas!(st::State, j::Int, k::Int)
    # 
    r1j, r2j = st.r1[j], st.r2[j]
    r1k, r2k = st.r1[k], st.r2[k]

    if (r1j == r1k) && (r2j == r2k)
        return false
    end

    Rj = r1j / r2j
    Rk = r1k / r2k

    eR1 = 0.0
    eR2 = 0.0

    if Rj > Rk
        α = r2k / r2j
        β = r1j / r1k
        p = (-1 + sqrt(α * β)) / (β + sqrt(α * β))
        q = (-1 + sqrt(α * β)) / (α + sqrt(α * β))
        eR1 = -p * r1j
        eR2 =  q * r2k
    else
        α = r1k / r1j
        β = r2j / r2k
        p = (-1 + sqrt(α * β)) / (β + sqrt(α * β))
        q = (-1 + sqrt(α * β)) / (α + sqrt(α * β))
        eR1 =  q * r1k
        eR2 = -p * r2j
    end

    st.r1[j] += eR1
    st.r2[j] += eR2
    st.r1[k] -= eR1
    st.r2[k] -= eR2
    return true
end

# ---------------------
# Interaction resolution between the pair of agents
# ---------------------

function resolve_pair!(st::State, cfg::Config, j::Int, k::Int, min_resource::Float64)
    aj = choose_action(st.b[j], st.λ[j], st.r1[j], st.r2[j])
    ak = choose_action(st.b[k], st.λ[k], st.r1[k], st.r2[k])

    war = 0.0
    raid = 0.0
    trade = 0.0
    avoid = 0.0

    if (aj == 1) && (ak == 1)
        # conflict-conflict (war): winner takes all
        winner, loser = rand(Bool) ? (j, k) : (k, j)
        st.r1[winner] += st.r1[loser]
        st.r2[winner] += st.r2[loser]
        st.r1[loser]  = -min_resource           # correspond to the cost of war paid on the payoff (r'1)
        st.r2[loser]  = -min_resource           # correspond to the cost of war paid on the payoff (r'2)
        war = 1.0

    elseif (aj == 1) && (ak != 1)
        can_raid = (cfg.avoidance == :constrained) || (cfg.avoidance == :normal && ak == 2)
        if can_raid
            st.r1[j] += st.r1[k];  st.r1[k] = 0.0
            st.r2[j] += st.r2[k];  st.r2[k] = 0.0
            raid = 1.0
        else
            avoid = 1.0
        end

    elseif (ak == 1) && (aj != 1)
        can_raid = (cfg.avoidance == :constrained) || (cfg.avoidance == :normal && aj == 2)
        if can_raid
            st.r1[k] += st.r1[j];  st.r1[j] = 0.0
            st.r2[k] += st.r2[j];  st.r2[j] = 0.0
            raid = 1.0
        else
            avoid = 1.0
        end


    elseif (aj == 2) && (ak == 2)
        ok = (cfg.trade == :minimal_law) ? trade_minimal_law!(st, j, k) :
             (cfg.trade == :cobb_douglas) ? trade_cobb_douglas!(st, j, k) :
             error("Unknown trade = $(cfg.trade)")
        trade = ok ? 1.0 : 0.0
        avoid = ok ? 0.0 : 1.0

    else
        avoid = 1.0
    end

    # This trick ensures that we cannot go have negative fitness value when cost are paid
    st.r1[j] += min_resource; st.r2[j] += min_resource
    st.r1[k] += min_resource; st.r2[k] += min_resource

    # fitness depends on trade production rule
    if cfg.trade == :minimal_law
        st.fitness[j] = min2(st.r1[j], st.r2[j])
        st.fitness[k] = min2(st.r1[k], st.r2[k])
    elseif cfg.trade == :cobb_douglas
        st.fitness[j] = sqrt(st.r1[j] * st.r2[j])
        st.fitness[k] = sqrt(st.r1[k] * st.r2[k])
    else
        error("Unknown trade = $(cfg.trade)")
    end


    return war, raid, trade, avoid
end

# ============================================================
# Evolution (mutation and reproduction)
# ============================================================

function reproduce_and_mutate!(st::State, cfg::Config, d_b::Binomial, d_l::Binomial)
    sample!(1:cfg.G, Weights(st.fitness), st.parents)
    st.b .= st.b[st.parents]
    st.λ .= st.λ[st.parents]

    nb_mut_b = rand(d_b)
    nb_mut_l = rand(d_l)

    if nb_mut_b > 0
        muts = sample(1:cfg.G, nb_mut_b; replace=false, ordered=true)
        for idx in muts
            st.b[idx] = clamp01(st.b[idx] + randn() * cfg.sigma_b)
        end
    end

    if nb_mut_l > 0
        muts = sample(1:cfg.G, nb_mut_l; replace=false, ordered=true)
        for idx in muts
            st.λ[idx] = clamp01(st.λ[idx] + randn() * cfg.sigma_l)
        end
    end
end

# ============================================================
# Recording
# ============================================================

mutable struct RunRecorder
    # timeseries sampled every record_every
    t::Vector{Int}
    mean_b::Vector{Float64}
    sd_b_within::Vector{Float64}
    mean_λ::Vector{Float64}
    sd_λ_within::Vector{Float64}
    mean_war::Vector{Float64}
    mean_raid::Vector{Float64}
    mean_trade::Vector{Float64}
    mean_avoid::Vector{Float64}

    # equilibrium window for boxplots (single value per run)
    eq_war_sum::Float64
    eq_raid_sum::Float64
    eq_trade_sum::Float64
    eq_avoid_sum::Float64
    eq_n::Int

    # snapshots (for hist/heatmap)
    snap_gen::Vector{Int}
    snap_b::Vector{Vector{Float64}}
    snap_λ::Vector{Vector{Float64}}
end

function RunRecorder()
    return RunRecorder(Int[], Float64[], Float64[], Float64[], Float64[],
                       Float64[], Float64[], Float64[], Float64[],
                       0.0, 0.0, 0.0, 0.0, 0,
                       Int[], Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
end

function maybe_record_timeseries!(rec::RunRecorder, st::State, cfg::Config, gen::Int,
                                 p_w::Float64, p_r::Float64, p_t::Float64, p_a::Float64)
    if gen % cfg.record_every != 0
        return
    end
    push!(rec.t, gen)
    push!(rec.mean_b, mean(st.b))
    push!(rec.sd_b_within, std(st.b))
    push!(rec.mean_λ, mean(st.λ))
    push!(rec.sd_λ_within, std(st.λ))
    push!(rec.mean_war, p_w)
    push!(rec.mean_raid, p_r)
    push!(rec.mean_trade, p_t)
    push!(rec.mean_avoid, p_a)
end

function maybe_accumulate_equilibrium!(rec::RunRecorder, cfg::Config, gen::Int,
                                      p_w::Float64, p_r::Float64, p_t::Float64, p_a::Float64)
    if (gen >= cfg.eq_start) && (gen <= cfg.eq_end)
        rec.eq_war_sum += p_w
        rec.eq_raid_sum += p_r
        rec.eq_trade_sum += p_t
        rec.eq_avoid_sum += p_a
        rec.eq_n += 1
    end
end

function maybe_snapshot!(rec::RunRecorder, st::State, cfg::Config, gen::Int)
    if gen < cfg.snapshot_start
        return
    end
    if (gen - cfg.snapshot_start) % cfg.snapshot_every != 0
        return
    end
    idx = sample(1:cfg.G, cfg.sample_size; replace=false)
    push!(rec.snap_gen, gen)
    push!(rec.snap_b, st.b[idx])
    push!(rec.snap_λ, st.λ[idx])
end

# ============================================================
# One run
# ============================================================

function run_one(cfg::Config; seed::Int=1)
    Random.seed!(seed)
    dists = make_env_dists()

    st = State(fill(0.5, cfg.G), zeros(cfg.G),
               fill(0.5, cfg.G), fill(0.5, cfg.G),
               zeros(cfg.G), zeros(Int, cfg.G))

    d_b = Binomial(cfg.G, cfg.mu_b)
    d_l = Binomial(cfg.G, cfg.mu_l)

    n_pairs = cfg.G ÷ 2
    min_resource = cfg.Cw / 2

    rec = RunRecorder()

    for gen in 1:cfg.T
        draw_environment!(st, cfg, dists)

        cnt_w = 0.0; cnt_r = 0.0; cnt_t = 0.0; cnt_a = 0.0

        @inbounds for j in 1:n_pairs
            k = j + n_pairs
            w, r, t, a = resolve_pair!(st, cfg, j, k, min_resource)
            cnt_w += w; cnt_r += r; cnt_t += t; cnt_a += a
        end

        p_w = cnt_w / n_pairs
        p_r = cnt_r / n_pairs
        p_t = cnt_t / n_pairs
        p_a = cnt_a / n_pairs

        maybe_record_timeseries!(rec, st, cfg, gen, p_w, p_r, p_t, p_a)
        maybe_accumulate_equilibrium!(rec, cfg, gen, p_w, p_r, p_t, p_a)
        maybe_snapshot!(rec, st, cfg, gen)

        reproduce_and_mutate!(st, cfg, d_b, d_l)
    end

    # equilibrium means for boxplot (single per run)
    eq_war = rec.eq_war_sum / max(rec.eq_n, 1)
    eq_raid = rec.eq_raid_sum / max(rec.eq_n, 1)
    eq_trade = rec.eq_trade_sum / max(rec.eq_n, 1)
    eq_avoid = rec.eq_avoid_sum / max(rec.eq_n, 1)

    return rec, (eq_war=eq_war, eq_raid=eq_raid, eq_trade=eq_trade, eq_avoid=eq_avoid)
end

# ============================================================
# Multiple runs + outputs
# ============================================================

function runs_to_timeseries_df(cfg::Config, recs::Vector{RunRecorder})
    # align on recorded timepoints (assumes same record_every, so same grid)
    tgrid = recs[1].t
    R = length(recs)
    M = length(tgrid)

    mean_over_runs(xgetter) = [mean([xgetter(recs[r])[i] for r in 1:R]) for i in 1:M]
    sd_over_runs(xgetter)   = [std([xgetter(recs[r])[i] for r in 1:R])  for i in 1:M]

    df = DataFrame(
        t = tgrid,
        mean_b = mean_over_runs(r -> r.mean_b),
        sd_b_runs = sd_over_runs(r -> r.mean_b),
        mean_sd_b_within = mean_over_runs(r -> r.sd_b_within),

        mean_λ = mean_over_runs(r -> r.mean_λ),
        sd_λ_runs = sd_over_runs(r -> r.mean_λ),
        mean_sd_λ_within = mean_over_runs(r -> r.sd_λ_within),

        mean_war = mean_over_runs(r -> r.mean_war),
        sd_war_runs = sd_over_runs(r -> r.mean_war),

        mean_raid = mean_over_runs(r -> r.mean_raid),
        sd_raid_runs = sd_over_runs(r -> r.mean_raid),

        mean_trade = mean_over_runs(r -> r.mean_trade),
        sd_trade_runs = sd_over_runs(r -> r.mean_trade),

        mean_avoid = mean_over_runs(r -> r.mean_avoid),
        sd_avoid_runs = sd_over_runs(r -> r.mean_avoid),
    )

    # metadata columns
    df.env .= string(cfg.env)
    df.trade .= string(cfg.trade)
    df.avoidance .= string(cfg.avoidance)
    df.Cw .= cfg.Cw
    df.G .= cfg.G
    df.T_total .= cfg.T

    return df
end

function runs_to_boxplot_df(cfg::Config, eqs::Vector{NamedTuple})
    df = DataFrame(run_id=Int[], war=Float64[], raid=Float64[], trade=Float64[], avoid=Float64[])
    for (i, e) in enumerate(eqs)
        push!(df, (run_id=i, war=e.eq_war, raid=e.eq_raid, trade=e.eq_trade, avoid=e.eq_avoid))
    end
    df.env .= string(cfg.env)
    df.trade .= string(cfg.trade)
    df.avoidance .= string(cfg.avoidance)
    df.Cw .= cfg.Cw
    return df
end

function runs_to_snapshots_df(cfg::Config, recs::Vector{RunRecorder})
    # long format: each row = one sampled individual in one snapshot in one run
    df = DataFrame(run_id=Int[], gen=Int[], b=Float64[], λ=Float64[])
    for (rid, rec) in enumerate(recs)
        for s in 1:length(rec.snap_gen)
            gen = rec.snap_gen[s]
            bs  = rec.snap_b[s]
            ls  = rec.snap_λ[s]
            for i in eachindex(bs)
                push!(df, (run_id=rid, gen=gen, b=bs[i], λ=ls[i]))
            end
        end
    end
    df.env .= string(cfg.env)
    df.trade .= string(cfg.trade)
    df.avoidance .= string(cfg.avoidance)
    df.Cw .= cfg.Cw
    return df
end

function run_experiment(cfg::Config)
    ensure_dir(cfg.outdir)

    recs = Vector{RunRecorder}(undef, cfg.num_runs)
    eqs  = Vector{NamedTuple}(undef, cfg.num_runs)

    for r in 1:cfg.num_runs
        recs[r], eqs[r] = run_one(cfg; seed=r)
        @printf("Run %d/%d done\n", r, cfg.num_runs)
    end

    df_ts = runs_to_timeseries_df(cfg, recs)
    df_box = runs_to_boxplot_df(cfg, eqs)
    df_snap = runs_to_snapshots_df(cfg, recs)

    tag = "$(cfg.env)_$(cfg.trade)_$(cfg.avoidance)_Cw$(cfg.Cw)"
    path_ts   = joinpath(cfg.outdir, "timeseries_$tag.csv")
    path_box  = joinpath(cfg.outdir, "boxplot_$tag.csv")
    path_snap = joinpath(cfg.outdir, "snapshots_$tag.csv")

    CSV.write(path_ts, df_ts)
    CSV.write(path_box, df_box)
    CSV.write(path_snap, df_snap)

    println("Wrote:")
    println(" - $path_ts")
    println(" - $path_box")
    println(" - $path_snap")

    return (timeseries=path_ts, boxplot=path_box, snapshots=path_snap)
end

# ============================================================
# Example usage
# ============================================================

cfg = Config(
    10_000, 200_000, 4.0,
    :richpoor,
    :minimal_law,
    :normal,
    0.01, 0.01, 0.02, 0.02,
    20,
    1000,
    180_001, 200_000,
    180_001, 1000, 2000,
    "out"
)

run_experiment(cfg)
