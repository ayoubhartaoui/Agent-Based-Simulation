# ============================================================
# Master Thesis – Simulation Code
#
# This script implements the evolutionary simulation model used in the Master thesis.
# 
#
# - Supports multiple environmental variability regimes
# - Supports minimal-law and Cobb–Douglas trade
# - Supports normal and constrained avoidance (for this see in the run simulation function line 193 and 200)
# - Produces CSV outputs for downstream R analysis
#
# Outputs are written to the `output/` directory.
const OUTPUT_DIR = "output"
isdir(OUTPUT_DIR) || mkdir(OUTPUT_DIR)
# ============================================================




# ==============================================
# Libraries
# ==============================================

using Random, StatsBase, Distributions
using DataFrames, CSV


# ==============================================
# Default parameters
# ==============================================

const DEFAULTS = (

    G = 10_000,            # population size (number of agents)

    T = 100_000,           # Number of generations

    Cw = 7.0,              # conflict cost (only a placeholder here)

    env_case = 1,          # environment: 1 = low variability | 2 = rich–poor | 3 = trade-off | 4 = uncorrelated
                        
    trade_case = 1,        # trade rule:  1 = minimal-law trade | 2 = Cobb–Douglas trade                 

    mu_b = 0.01,           # mutation rate for belligerence (b)

    mu_h = 0.01,           # mutation rate for dissatisfaction threshold (λ)

    sigma_b = 0.02,        # mutation step size for belligerence

    sigma_h = 0.02,        # mutation step size for dissatisfaction threshold

    num_runs = 10,         # number of replicate 

    window = 10_000,       # Size of the final time window used for: histograms | heatmaps | boxplots

    sample_size = 2000     # individuals sampled per generation in final window (histograms and heatmaps)
)



# ==============================================
# Simulation functions (main code)
# ==============================================

function run_simulation(
    G::Int64,
    T::Int64,
    Cw::Float64,
    env_case::Int64,
    trade_case::Int64,
    mu_b::Float64,
    mu_h::Float64,
    sigma_b::Float64,
    sigma_h::Float64;
    window::Int64 ,
    sample_size::Int64 ,
)
    bias    = fill(0.5, G)
    lambda  = zeros(G)
    r1      = fill(0.5, G)
    r2      = fill(0.5, G)
    fitness = zeros(G)
    parents = zeros(Int, G)

    d = Binomial(G, mu_b)
    e = Binomial(G, mu_h)

    n_pairs = Int(G ÷ 2)

    mean_bias_hist    = zeros(T)
    mean_lambda_hist  = zeros(T)
    mean_war_hist     = zeros(T)
    mean_raid_hist    = zeros(T)
    mean_trade_hist   = zeros(T)
    mean_inac_hist    = zeros(T)

    std_bias_hist       = zeros(T)
    std_lambda_hist     = zeros(T)
    std_war_within      = zeros(T)
    std_raid_within     = zeros(T)
    std_trade_within    = zeros(T)
    std_inaction_within = zeros(T)

    min_resource = Cw

    beta_rich = Beta(15.0, 1.5)
    beta_poor = Beta(1.5, 15.0)
    high_dist = Beta(5.0, 1.5)
    low_dist  = Beta(1.5, 5.0)
    n_lowvar  = truncated(Normal(0.5, 0.1), 0.0, 1.0)

    w = min(window, T)
    start_capture = T - w + 1

    sample_size_eff = min(sample_size, G)
    bias_sample_mat   = Matrix{Float32}(undef, w, sample_size_eff)
    lambda_sample_mat = Matrix{Float32}(undef, w, sample_size_eff)

    war_last   = zeros(Float32, w)
    raid_last  = zeros(Float32, w)
    trade_last = zeros(Float32, w)
    inac_last  = zeros(Float32, w)

    for t in 1:T
        if env_case == 1
            for i in 1:G
                r1[i] = rand(n_lowvar)
                r2[i] = rand(n_lowvar)
                fitness[i] = 0.0
            end
        elseif env_case == 2
            for i in 1:G
                r1[i] = rand()
                r2[i] = rand()
                fitness[i] = 0.0
            end
        elseif env_case == 3
            n_rich = G ÷ 2
            is_rich = falses(G)
            is_rich[sample(1:G, n_rich, replace=false)] .= true
            for i in 1:G
                if is_rich[i]
                    r1[i] = rand(beta_rich)
                    r2[i] = rand(beta_rich)
                else
                    r1[i] = rand(beta_poor)
                    r2[i] = rand(beta_poor)
                end
                fitness[i] = 0.0
            end
        elseif env_case == 4
            for i in 1:G
                if rand() < 0.5
                    r1[i] = rand(high_dist)
                    r2[i] = rand(low_dist)
                else
                    r1[i] = rand(low_dist)
                    r2[i] = rand(high_dist)
                end
                fitness[i] = 0.0
            end
        else
            error("env_case must be 1, 2, 3, or 4. Got: $env_case")
        end

        sum_bias   = 0.0
        sum_lambda = 0.0
        cnt_w      = 0.0
        cnt_t      = 0.0
        cnt_r      = 0.0
        cnt_inac   = 0.0

        for j in 1:n_pairs
            k = j + n_pairs

            bias_j = bias[j]
            bias_k = bias[k]

            sum_bias   += bias_j + bias_k
            sum_lambda += lambda[j] + lambda[k]

            r_to_compare_j = r1[j] < r2[j] ? r1[j] : r2[j]
            action_j = lambda[j] >= r_to_compare_j ? (rand() < bias_j ? 1 : 2) : 0

            r_to_compare_k = r1[k] < r2[k] ? r1[k] : r2[k]
            action_k = lambda[k] >= r_to_compare_k ? (rand() < bias_k ? 1 : 2) : 0

            if action_j == 1 && action_k == 1
                winner, loser = rand(Bool) ? (j, k) : (k, j)
                r1[winner] += r1[loser]
                r2[winner] += r2[loser]
                r1[loser] = -min_resource
                r2[loser] = -min_resource
                cnt_w += 1

            elseif action_j == 1 && action_k == 2       # for the constrained case just get rid of && action_k == 2
                r1[j] += r1[k]
                r2[j] += r2[k]
                r1[k] = 0.0
                r2[k] = 0.0
                cnt_r += 1

            elseif action_k == 1 && action_j == 2       # for the constrained case just get rid of && action_k == 2
                r1[k] += r1[j]
                r2[k] += r2[j]
                r1[j] = 0.0
                r2[j] = 0.0
                cnt_r += 1

            elseif action_j == 2 && action_k == 2
                if trade_case == 1
                    R1_1 = r1[j]
                    R2_1 = r2[j]
                    R1_2 = r1[k]
                    R2_2 = r2[k]

                    surplus1 = max(R1_1, R2_1) - min(R1_1, R2_1)
                    surplus2 = max(R1_2, R2_2) - min(R1_2, R2_2)

                    if R2_1 > R1_1 && R1_2 > R2_2
                        receive1 = :r_1
                        receive2 = :r_2
                        denom = min(R1_2, R2_2) + min(R1_1, R2_1)
                        t1 = denom == 0.0 ? 0.0 : min(R1_2, R2_2) / denom
                        t2 = 1 - t1
                    elseif R1_1 > R2_1 && R2_2 > R1_2
                        receive1 = :r_2
                        receive2 = :r_1
                        denom = max(R1_2, R2_2) + max(R1_1, R2_1)
                        t1 = denom == 0.0 ? 0.0 : max(R1_2, R2_2) / denom
                        t2 = 1 - t1
                    else
                        receive1 = :r_1
                        receive2 = :r_2
                        t1 = 0.0
                        t2 = 0.0
                    end

                    g1_gives = surplus1 * t1
                    g2_gives = surplus2 * t2

                    new_R1_1 = receive1 == :r_1 ? R1_1 + g2_gives : R1_1 - g1_gives
                    new_R2_1 = receive1 == :r_2 ? R2_1 + g2_gives : R2_1 - g1_gives
                    new_R1_2 = receive2 == :r_1 ? R1_2 + g1_gives : R1_2 - g2_gives
                    new_R2_2 = receive2 == :r_2 ? R2_2 + g1_gives : R2_2 - g1_gives

                    min_r1_before = min(R1_1, R2_1)
                    min_r2_before = min(R1_2, R2_2)
                    min_r1_after  = min(new_R1_1, new_R2_1)
                    min_r2_after  = min(new_R1_2, new_R2_2)

                    if !(min_r1_after <= min_r1_before || min_r2_after <= min_r2_before)
                        r1[j] = new_R1_1
                        r2[j] = new_R2_1
                        r1[k] = new_R1_2
                        r2[k] = new_R2_2
                    end

                    cnt_t += 1

                elseif trade_case == 2
                    if r1[j] == r1[k] && r2[j] == r2[k]
                        eR1 = 0.0
                        eR2 = 0.0
                    else
                        Rj = r1[j] / r2[j]
                        Rk = r1[k] / r2[k]

                        if Rj > Rk
                            alpha = r2[k] / r2[j]
                            beta  = r1[j] / r1[k]
                            p = (-1 + sqrt(alpha * beta)) / (beta + sqrt(alpha * beta))
                            q = (-1 + sqrt(alpha * beta)) / (alpha + sqrt(alpha * beta))
                            eR1 = -p * r1[j]
                            eR2 =  q * r2[k]
                        else
                            alpha = r1[k] / r1[j]
                            beta  = r2[j] / r2[k]
                            p = (-1 + sqrt(alpha * beta)) / (beta + sqrt(alpha * beta))
                            q = (-1 + sqrt(alpha * beta)) / (alpha + sqrt(alpha * beta))
                            eR1 =  q * r1[k]
                            eR2 = -p * r2[j]
                        end
                    end

                    r1[j] += eR1
                    r2[j] += eR2
                    r1[k] -= eR1
                    r2[k] -= eR2

                    cnt_t += 1
                else
                    error("trade_case must be 1 or 2. Got: $trade_case")
                end
            else
                cnt_inac += 1
            end

            r1[j] += min_resource
            r2[j] += min_resource
            r1[k] += min_resource
            r2[k] += min_resource

            if trade_case == 2
                fitness[j] = sqrt(r1[j] * r2[j])
                fitness[k] = sqrt(r1[k] * r2[k])
            else
                fitness[j] = r1[j] < r2[j] ? r1[j] : r2[j]
                fitness[k] = r1[k] < r2[k] ? r1[k] : r2[k]
            end
        end

        mean_bias_hist[t]   = sum_bias / G
        std_bias_hist[t]    = std(bias)
        mean_lambda_hist[t] = sum_lambda / G
        std_lambda_hist[t]  = std(lambda)

        p_w    = cnt_w / n_pairs
        p_raid = cnt_r / n_pairs
        p_t    = cnt_t / n_pairs
        p_inac = cnt_inac / n_pairs

        mean_war_hist[t]   = p_w
        mean_raid_hist[t]  = p_raid
        mean_trade_hist[t] = p_t
        mean_inac_hist[t]  = p_inac

        std_war_within[t]      = sqrt(p_w * (1 - p_w))
        std_raid_within[t]     = sqrt(p_raid * (1 - p_raid))
        std_trade_within[t]    = sqrt(p_t * (1 - p_t))
        std_inaction_within[t] = sqrt(p_inac * (1 - p_inac))

        if t >= start_capture
            row = t - start_capture + 1
            idxs = sample(1:G, sample_size_eff, replace=false, ordered=true)
            bias_sample_mat[row, :]   .= bias[idxs]
            lambda_sample_mat[row, :] .= lambda[idxs]
            war_last[row]   = p_w
            raid_last[row]  = p_raid
            trade_last[row] = p_t
            inac_last[row]  = p_inac
        end

        sample!(1:G, Weights(fitness), parents)
        bias   = bias[parents]
        lambda = lambda[parents]

        nb_mutants_bias = rand(d)
        nb_mutants_h    = rand(e)

        if nb_mutants_bias > 0
            mutants_bias = sample(1:G, nb_mutants_bias, replace=false, ordered=true)
            for curr_mutant in mutants_bias
                bias[curr_mutant] = clamp(bias[curr_mutant] + randn() * sigma_b, 0.0, 1.0)
            end
        end

        if nb_mutants_h > 0
            mutants_h = sample(1:G, nb_mutants_h, replace=false, ordered=true)
            for curr_mutant in mutants_h
                lambda[curr_mutant] = clamp(lambda[curr_mutant] + randn() * sigma_h, 0.0, 1.0)
            end
        end
    end

    return (
        bias               = mean_bias_hist,
        lambda             = mean_lambda_hist,
        war                = mean_war_hist,
        raid               = mean_raid_hist,
        trade              = mean_trade_hist,
        inaction           = mean_inac_hist,
        sd_bias_within     = std_bias_hist,
        sd_lambda_within   = std_lambda_hist,
        sd_war_within      = std_war_within,
        sd_raid_within     = std_raid_within,
        sd_trade_within    = std_trade_within,
        sd_inaction_within = std_inaction_within,
        bias_last_window   = bias_sample_mat,
        lambda_last_window = lambda_sample_mat,
        war_last_window    = war_last,
        raid_last_window   = raid_last,
        trade_last_window  = trade_last,
        inaction_last_window = inac_last,
        window             = w,
        sample_size        = sample_size_eff,
    )
end



# =============================================
# Multiple runs for the replicates
# =============================================

function run_multiple(
    num_runs::Int;
    G::Int,
    T::Int,
    Cw::Float64,
    env_case::Int,
    trade_case::Int,
    mu_b::Float64,
    mu_h::Float64,
    sigma_b::Float64,
    sigma_h::Float64,
    window::Int ,
    sample_size::Int ,
)
    hist = [
        run_simulation(
            G, T, Cw, env_case, trade_case,
            mu_b, mu_h, sigma_b, sigma_h;
            window=window, sample_size=sample_size
        )
        for _ in 1:num_runs
    ]

    TT = length(hist[1].bias)

    all_bias     = zeros(TT, num_runs)
    all_lambda   = zeros(TT, num_runs)
    all_war      = zeros(TT, num_runs)
    all_raid     = zeros(TT, num_runs)
    all_trade    = zeros(TT, num_runs)
    all_inaction = zeros(TT, num_runs)

    all_sd_bias_within     = zeros(TT, num_runs)
    all_sd_lambda_within   = zeros(TT, num_runs)
    all_sd_war_within      = zeros(TT, num_runs)
    all_sd_raid_within     = zeros(TT, num_runs)
    all_sd_trade_within    = zeros(TT, num_runs)
    all_sd_inaction_within = zeros(TT, num_runs)

    for (k, h) in enumerate(hist)
        all_bias[:, k]     .= h.bias
        all_lambda[:, k]   .= h.lambda
        all_war[:, k]      .= h.war
        all_raid[:, k]     .= h.raid
        all_trade[:, k]    .= h.trade
        all_inaction[:, k] .= h.inaction

        all_sd_bias_within[:, k]     .= h.sd_bias_within
        all_sd_lambda_within[:, k]   .= h.sd_lambda_within
        all_sd_war_within[:, k]      .= h.sd_war_within
        all_sd_raid_within[:, k]     .= h.sd_raid_within
        all_sd_trade_within[:, k]    .= h.sd_trade_within
        all_sd_inaction_within[:, k] .= h.sd_inaction_within
    end

    mean_bias_over_runs     = vec(mean(all_bias, dims=2))
    sd_bias_over_runs       = vec(std(all_bias, dims=2))
    mean_lambda_over_runs   = vec(mean(all_lambda, dims=2))
    sd_lambda_over_runs     = vec(std(all_lambda, dims=2))

    mean_war_over_runs      = vec(mean(all_war, dims=2))
    sd_war_over_runs        = vec(std(all_war, dims=2))
    mean_raid_over_runs     = vec(mean(all_raid, dims=2))
    sd_raid_over_runs       = vec(std(all_raid, dims=2))
    mean_trade_over_runs    = vec(mean(all_trade, dims=2))
    sd_trade_over_runs      = vec(std(all_trade, dims=2))
    mean_inaction_over_runs = vec(mean(all_inaction, dims=2))
    sd_inaction_over_runs   = vec(std(all_inaction, dims=2))

    mean_sd_bias_within_over_runs     = vec(mean(all_sd_bias_within, dims=2))
    mean_sd_lambda_within_over_runs   = vec(mean(all_sd_lambda_within, dims=2))
    mean_sd_war_within_over_runs      = vec(mean(all_sd_war_within, dims=2))
    mean_sd_raid_within_over_runs     = vec(mean(all_sd_raid_within, dims=2))
    mean_sd_trade_within_over_runs    = vec(mean(all_sd_trade_within, dims=2))
    mean_sd_inaction_within_over_runs = vec(mean(all_sd_inaction_within, dims=2))

    return hist,
           mean_bias_over_runs, sd_bias_over_runs,
           mean_lambda_over_runs, sd_lambda_over_runs,
           mean_war_over_runs, sd_war_over_runs,
           mean_raid_over_runs, sd_raid_over_runs,
           mean_trade_over_runs, sd_trade_over_runs,
           mean_inaction_over_runs, sd_inaction_over_runs,
           mean_sd_bias_within_over_runs,
           mean_sd_lambda_within_over_runs,
           mean_sd_war_within_over_runs,
           mean_sd_raid_within_over_runs,
           mean_sd_trade_within_over_runs,
           mean_sd_inaction_within_over_runs
end



# =============================================
# Write to CSV functions
# =============================================

function write_timeseries_csv(
    filepath::AbstractString;
    G::Int,
    T::Int,
    Cw::Float64,
    env_case::Int,
    trade_case::Int,
    mu_b::Float64,
    mu_h::Float64,
    sigma_b::Float64,
    sigma_h::Float64,
    num_runs::Int,
    window::Int,
    sample_size::Int,
)
    hist,
    mean_bias, sd_bias,
    mean_lambda, sd_lambda,
    mean_war, sd_war,
    mean_raid, sd_raid,
    mean_trade, sd_trade,
    mean_inaction, sd_inaction,
    sd_bias_within,
    sd_lambda_within,
    sd_war_within,
    sd_raid_within,
    sd_trade_within,
    sd_inaction_within =
        run_multiple(
            num_runs;
            G=G, T=T, Cw=Cw,
            env_case=env_case, trade_case=trade_case,
            mu_b=mu_b, mu_h=mu_h, sigma_b=sigma_b, sigma_h=sigma_h,
            window=window, sample_size=sample_size
        )

    TT = length(mean_bias)
    df = DataFrame(
        t = collect(1:TT),
        mean_bias = Float32.(mean_bias),
        sd_bias_runs = Float32.(sd_bias),
        sd_bias_within = Float32.(sd_bias_within),
        mean_lambda = Float32.(mean_lambda),
        sd_lambda_runs = Float32.(sd_lambda),
        sd_lambda_within = Float32.(sd_lambda_within),
        mean_war = Float32.(mean_war),
        sd_war_runs = Float32.(sd_war),
        sd_war_within = Float32.(sd_war_within),
        mean_raid = Float32.(mean_raid),
        sd_raid_runs = Float32.(sd_raid),
        sd_raid_within = Float32.(sd_raid_within),
        mean_trade = Float32.(mean_trade),
        sd_trade_runs = Float32.(sd_trade),
        sd_trade_within = Float32.(sd_trade_within),
        mean_inaction = Float32.(mean_inaction),
        sd_inaction_runs = Float32.(sd_inaction),
        sd_inaction_within = Float32.(sd_inaction_within),
    )

    CSV.write(filepath, df)
    return filepath
end

function write_lastwindow_samples_csv(
    filepath::AbstractString,
    hist;
)
    if isfile(filepath)
        rm(filepath)
    end

    first = true
    for (run_id, h) in enumerate(hist)
        w = h.window
        s = h.sample_size
        bias_mat = h.bias_last_window
        lam_mat  = h.lambda_last_window

        run_col  = repeat([run_id], w * s)
        gen_col  = repeat(collect(1:w), inner=s)
        sid_col  = repeat(collect(1:s), w)
        bias_col = vec(bias_mat)
        lam_col  = vec(lam_mat)

        df = DataFrame(
            run_id = run_col,
            gen_in_window = gen_col,
            sample_id = sid_col,
            bias = Float32.(bias_col),
            lambda = Float32.(lam_col),
        )

        if first
            CSV.write(filepath, df)
            first = false
        else
            CSV.write(filepath, df; append=true, writeheader=false)
        end
    end

    return filepath
end

function write_boxplot_csv(
    filepath::AbstractString;
    G::Int,
    T::Int,
    Cw::Float64,
    env_case::Int,
    trade_case::Int,
    mu_b::Float64,
    mu_h::Float64,
    sigma_b::Float64,
    sigma_h::Float64,
    num_runs::Int,
    window::Int,
    sample_size::Int,
)
    hist,
    mean_bias, sd_bias,
    mean_lambda, sd_lambda,
    mean_war, sd_war,
    mean_raid, sd_raid,
    mean_trade, sd_trade,
    mean_inaction, sd_inaction,
    sd_bias_within,
    sd_lambda_within,
    sd_war_within,
    sd_raid_within,
    sd_trade_within,
    sd_inaction_within =
        run_multiple(
            num_runs;
            G=G, T=T, Cw=Cw,
            env_case=env_case, trade_case=trade_case,
            mu_b=mu_b, mu_h=mu_h, sigma_b=sigma_b, sigma_h=sigma_h,
            window=window, sample_size=sample_size
        )

    df = DataFrame(
        Cw = Float64[],
        run_id = Int[],
        war_mean = Float64[],
        raid_mean = Float64[],
        trade_mean = Float64[],
        inaction_mean = Float64[],
        war_sd_time = Float64[],
        raid_sd_time = Float64[],
        trade_sd_time = Float64[],
        inaction_sd_time = Float64[],
        bias_within_mean = Float64[],
        lambda_within_mean = Float64[],
    )

    for (run_id, h) in enumerate(hist)
        warv   = h.war_last_window
        raidv  = h.raid_last_window
        tradev = h.trade_last_window
        inacv  = h.inaction_last_window

        push!(df, (
            Cw = Cw,
            run_id = run_id,
            war_mean = mean(warv),
            raid_mean = mean(raidv),
            trade_mean = mean(tradev),
            inaction_mean = mean(inacv),
            war_sd_time = std(warv),
            raid_sd_time = std(raidv),
            trade_sd_time = std(tradev),
            inaction_sd_time = std(inacv),
            bias_within_mean = mean(h.sd_bias_within[(end - h.window + 1):end]),
            lambda_within_mean = mean(h.sd_lambda_within[(end - h.window + 1):end]),
        ))
    end

    CSV.write(filepath, df)
    return filepath, hist
end

# =============================================
# Run the different Conflict Costs (the only function to call))
# =============================================

function run_scenarios(scenarios::AbstractVector{<:NamedTuple}, base_label::AbstractString)
    out_timeseries = joinpath(OUTPUT_DIR, "$(base_label)_timeseries_summary_ALL.csv")
    out_boxplot    = joinpath(OUTPUT_DIR, "$(base_label)_boxplot_lastwindow_ALL.csv")
    out_samples    = joinpath(OUTPUT_DIR, "$(base_label)_lastwindow_bias_lambda_samples_ALL.csv")


    if isfile(out_timeseries); rm(out_timeseries); end
    if isfile(out_boxplot); rm(out_boxplot); end
    if isfile(out_samples); rm(out_samples); end

    df_timeseries = DataFrame(
        Cw = Float32[],
        t = Int32[],
        mean_bias = Float32[],
        sd_bias_runs = Float32[],
        sd_bias_within = Float32[],
        mean_lambda = Float32[],
        sd_lambda_runs = Float32[],
        sd_lambda_within = Float32[],
        mean_war = Float32[],
        sd_war_runs = Float32[],
        sd_war_within = Float32[],
        mean_raid = Float32[],
        sd_raid_runs = Float32[],
        sd_raid_within = Float32[],
        mean_trade = Float32[],
        sd_trade_runs = Float32[],
        sd_trade_within = Float32[],
        mean_inaction = Float32[],
        sd_inaction_runs = Float32[],
        sd_inaction_within = Float32[],
    )

    df_boxplot = DataFrame(
        Cw = Float32[],
        run_id = Int32[],
        war_mean = Float32[],
        raid_mean = Float32[],
        trade_mean = Float32[],
        inaction_mean = Float32[],
        war_sd_time = Float32[],
        raid_sd_time = Float32[],
        trade_sd_time = Float32[],
        inaction_sd_time = Float32[],
        bias_within_mean = Float32[],
        lambda_within_mean = Float32[],
    )

    first_samples_write = true

    for sc in scenarios
        cw = Float32(get(sc, :Cw, DEFAULTS.Cw))

        hist,
        mean_bias, sd_bias,
        mean_lambda, sd_lambda,
        mean_war, sd_war,
        mean_raid, sd_raid,
        mean_trade, sd_trade,
        mean_inaction, sd_inaction,
        sd_bias_within,
        sd_lambda_within,
        sd_war_within,
        sd_raid_within,
        sd_trade_within,
        sd_inaction_within =
            run_multiple(
                DEFAULTS.num_runs;
                G=DEFAULTS.G,
                T=DEFAULTS.T,
                Cw=Float64(cw),
                env_case=DEFAULTS.env_case,
                trade_case=DEFAULTS.trade_case,
                mu_b=DEFAULTS.mu_b,
                mu_h=DEFAULTS.mu_h,
                sigma_b=DEFAULTS.sigma_b,
                sigma_h=DEFAULTS.sigma_h,
                window=DEFAULTS.window,
                sample_size=DEFAULTS.sample_size
            )

        TT = length(mean_bias)

        append!(df_timeseries, DataFrame(
            Cw = fill(cw, TT),
            t = Int32.(collect(1:TT)),
            mean_bias = Float32.(mean_bias),
            sd_bias_runs = Float32.(sd_bias),
            sd_bias_within = Float32.(sd_bias_within),
            mean_lambda = Float32.(mean_lambda),
            sd_lambda_runs = Float32.(sd_lambda),
            sd_lambda_within = Float32.(sd_lambda_within),
            mean_war = Float32.(mean_war),
            sd_war_runs = Float32.(sd_war),
            sd_war_within = Float32.(sd_war_within),
            mean_raid = Float32.(mean_raid),
            sd_raid_runs = Float32.(sd_raid),
            sd_raid_within = Float32.(sd_raid_within),
            mean_trade = Float32.(mean_trade),
            sd_trade_runs = Float32.(sd_trade),
            sd_trade_within = Float32.(sd_trade_within),
            mean_inaction = Float32.(mean_inaction),
            sd_inaction_runs = Float32.(sd_inaction),
            sd_inaction_within = Float32.(sd_inaction_within),
        ))

        for (run_id, h) in enumerate(hist)
            append!(df_boxplot, DataFrame(
                Cw = cw,
                run_id = Int32(run_id),
                war_mean = Float32(mean(h.war_last_window)),
                raid_mean = Float32(mean(h.raid_last_window)),
                trade_mean = Float32(mean(h.trade_last_window)),
                inaction_mean = Float32(mean(h.inaction_last_window)),
                war_sd_time = Float32(std(h.war_last_window)),
                raid_sd_time = Float32(std(h.raid_last_window)),
                trade_sd_time = Float32(std(h.trade_last_window)),
                inaction_sd_time = Float32(std(h.inaction_last_window)),
                bias_within_mean = Float32(mean(h.sd_bias_within[end - h.window + 1:end])),
                lambda_within_mean = Float32(mean(h.sd_lambda_within[end - h.window + 1:end])),
            ))
        end

        for (run_id, h) in enumerate(hist)
            w = h.window
            s = h.sample_size
            n = w * s

            df_sp = DataFrame(
                Cw = fill(cw, n),
                run_id = fill(Int32(run_id), n),
                gen_in_window = Int32.(repeat(collect(1:w), inner=s)),
                sample_id = Int32.(repeat(collect(1:s), w)),
                bias = Float32.(vec(h.bias_last_window)),
                lambda = Float32.(vec(h.lambda_last_window)),
            )

            if first_samples_write
                CSV.write(out_samples, df_sp; floatformat="%.6g")
                first_samples_write = false
            else
                CSV.write(out_samples, df_sp; append=true, writeheader=false, floatformat="%.6g")
            end
        end
    end

    CSV.write(out_timeseries, df_timeseries; floatformat="%.6g")
    CSV.write(out_boxplot, df_boxplot; floatformat="%.6g")

    return out_timeseries, out_boxplot, out_samples
end

# ------ Choice of the differnet Conflict Costs (Cw) (High costs corresponds to high violence or physiological costs) ------
SCENARIOS_INPUT = [
    (Cw = 1.25,),
    (Cw = 4.0,),
    (Cw = 7.0,),
]

trade_label = DEFAULTS.trade_case == 1 ? "Minimal_law" : "Cobb_Douglas"

env_label =
    DEFAULTS.env_case == 1 ? "low_variability" :
    DEFAULTS.env_case == 2 ? "uncorrelated" :
    DEFAULTS.env_case == 3 ? "rich_poor" :
    DEFAULTS.env_case == 4 ? "trade_off" :
    "env$(DEFAULTS.env_case)"

base_label = "$(trade_label)_$(env_label)"


# =============================================
# Example of a call
# =============================================
out_ts, out_bp, out_sp = run_scenarios(SCENARIOS_INPUT, base_label)

println("Wrote: $out_ts")
println("Wrote: $out_bp")
println("Wrote: $out_sp")
println("Simulation complete.")
