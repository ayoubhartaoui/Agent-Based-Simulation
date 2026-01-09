using Random, StatsBase, Distributions
using DataFrames
using CSV
using Printf


function run_simulation(G::Int64, T::Int64, Cw::Float64, variable_sim::Int64,
                        mu_b::Float64, mu_h::Float64, sigma_b::Float64, sigma_h::Float64)

    bias   = fill(0.5, G)
    lambda = zeros(G)
    r1 = fill(0.5, G)
    r2 = fill(0.5, G)
    fitness = zeros(G)
    parents = zeros(Int, G)

    d = Binomial(G, mu_b)
    e = Binomial(G, mu_h)

    n_pairs = Int(G ÷ 2)

    mean_bias_hist   = zeros(T)
    mean_lambda_hist = zeros(T)
    mean_war_hist    = zeros(T)
    mean_raid_hist   = zeros(T)
    mean_trade_hist  = zeros(T)
    mean_inac_hist   = zeros(T)
    v_mean_per_gen   = zeros(T)

    std_bias_hist      = zeros(T)  # within-pop SD(bias)
    std_lambda_hist    = zeros(T)  # within-pop SD(lambda)
    std_war_within     = zeros(T)  # within-pop SD(war)
    std_raid_within    = zeros(T)  # within-pop SD(raid)
    std_trade_within   = zeros(T)  # within-pop SD(trade)
    std_inaction_within = zeros(T) # within-pop SD(inaction)

    min_resource = Cw / 2

    beta_rich = Beta(15.0, 1.5)   
    beta_poor = Beta(1.5, 15.0)

    for t in 1:T
        if variable_sim == 1
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
        end


        sum_bias   = 0.0
        sum_lambda = 0.0
        cnt_w      = 0.0
        cnt_t      = 0.0
        cnt_r      = 0.0
        cnt_inac   = 0.0
        v          = 0.0

        for j in 1:n_pairs
            k = j + n_pairs

            bias_j = bias[j]
            bias_k = bias[k]

            sum_bias   += bias_j + bias_k
            sum_lambda += lambda[j] + lambda[k]

            r_to_compare_j = r1[j] < r2[j] ? r1[j] : r2[j]
            v += r_to_compare_j

            if lambda[j] >= r_to_compare_j
                if rand() < bias_j
                    action_j = 1
                else
                    action_j = 2
                end
            else
                action_j = 0
            end

            r_to_compare_k = r1[k] < r2[k] ? r1[k] : r2[k]

            if lambda[k] >= r_to_compare_k
                if rand() < bias_k
                    action_k = 1
                else
                    action_k = 2
                end
            else
                action_k = 0
            end

            if action_j == 1 && action_k == 1
                winner, loser = rand(Bool) ? (j, k) : (k, j)
                r1[winner] += r1[loser]
                r1[loser]  = -min_resource
                r2[winner] += r2[loser]
                r2[loser]  = -min_resource
                cnt_w += 1

            elseif action_j == 1
                r1[j] += r1[k]
                r1[k]  = 0.0
                r2[j] += r2[k]
                r2[k]  = 0.0
                cnt_r += 1

            elseif action_k == 1
                r1[k] += r1[j]
                r1[j]  = 0.0
                r2[k] += r2[j]
                r2[j]  = 0.0
                cnt_r += 1

            elseif action_j == 2 && action_k == 2
                R1_1 = r1[j]
                R2_1 = r2[j]
                R1_2 = r1[k]
                R2_2 = r2[k]

                surplus1 = max(R1_1, R2_1) - min(R1_1, R2_1)
                surplus2 = max(R1_2, R2_2) - min(R1_2, R2_2)

                if R2_1 > R1_1 && R1_2 > R2_2
                    receive1 = :r_1
                    receive2 = :r_2
                    t1 = min(R1_2, R2_2) / (min(R1_2, R2_2) + min(R1_1, R2_1))
                    t2 = 1 - t1
                elseif R1_1 > R2_1 && R2_2 > R1_2
                    receive1 = :r_2
                    receive2 = :r_1
                    t1 = max(R1_2, R2_2) / (max(R1_2, R2_2) + max(R1_1, R2_1))
                    t2 = 1 - t1
                else
                    t1 = 0.0
                    t2 = 0.0
                    receive1 = :r_1
                    receive2 = :r_2
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

                if min_r1_after <= min_r1_before || min_r2_after <= min_r2_before
                    r1[j] += 0.0
                    r1[k] += 0.0
                    r2[j] += 0.0
                    r2[k] += 0.0
                end

                r1[j] = new_R1_1
                r1[k] = new_R1_2
                r2[j] = new_R2_1
                r2[k] = new_R2_2

                cnt_t += 1

            else
                cnt_inac += 1
            end

            r1[j] += min_resource
            r2[j] += min_resource
            r1[k] += min_resource
            r2[k] += min_resource

            fitness[j] = r1[j] < r2[j] ? r1[j] : r2[j]
            fitness[k] = r1[k] < r2[k] ? r1[k] : r2[k]
        end

        mean_bias_hist[t]   = sum_bias / G
        std_bias_hist[t]    = std(bias)
        mean_lambda_hist[t] = sum_lambda / G
        std_lambda_hist[t]  = std(lambda)

        p_w    = cnt_w / n_pairs
        p_raid = cnt_r / n_pairs
        p_t    = cnt_t / n_pairs
        p_inac = cnt_inac / n_pairs

        mean_war_hist[t]    = p_w
        mean_raid_hist[t]   = p_raid
        mean_trade_hist[t]  = p_t
        mean_inac_hist[t]   = p_inac

        std_war_within[t]      = sqrt(p_w * (1 - p_w))
        std_raid_within[t]     = sqrt(p_raid * (1 - p_raid))
        std_trade_within[t]    = sqrt(p_t * (1 - p_t))
        std_inaction_within[t] = sqrt(p_inac * (1 - p_inac))

        v_mean_per_gen[t]   = v / G

        sample!(1:G, Weights(fitness), parents)
        bias   = bias[parents]
        lambda = lambda[parents]

        nb_mutants_bias = rand(d)
        nb_mutants_h    = rand(e)

        if nb_mutants_bias > 0
            mutants_bias = sample(1:G, nb_mutants_bias, replace=false, ordered=true)
            for curr_mut = 1:nb_mutants_bias
                curr_mutant = mutants_bias[curr_mut]
                bias[curr_mutant] = clamp(bias[curr_mutant] + randn() * sigma_b, 0.0, 1.0)
            end
        end

        if nb_mutants_h > 0
            mutants_h = sample(1:G, nb_mutants_h, replace=false, ordered=true)
            for curr_mut = 1:nb_mutants_h
                curr_mutant = mutants_h[curr_mut]
                lambda[curr_mutant] = clamp(lambda[curr_mutant] + randn() * sigma_h, 0.0, 1.0)
            end
        end
    end

    return (
        bias              = mean_bias_hist,
        lambda            = mean_lambda_hist,
        war               = mean_war_hist,
        raid              = mean_raid_hist,
        trade             = mean_trade_hist,
        inaction          = mean_inac_hist,
        sd_bias_within    = std_bias_hist,
        sd_lambda_within  = std_lambda_hist,
        sd_war_within     = std_war_within,
        sd_raid_within    = std_raid_within,
        sd_trade_within   = std_trade_within,
        sd_inaction_within = std_inaction_within
    )
end

function run_multiple(num_runs::Int; G, T, Cw, variable_sim, mu_b, mu_h, sigma_b, sigma_h)
    hist = [run_simulation(G, T, Cw, variable_sim, mu_b, mu_h, sigma_b, sigma_h)
            for _ in 1:num_runs]

    T = length(hist[1].bias)

    all_bias      = zeros(T, num_runs)
    all_lambda    = zeros(T, num_runs)
    all_war       = zeros(T, num_runs)
    all_raid      = zeros(T, num_runs)
    all_trade     = zeros(T, num_runs)
    all_inaction  = zeros(T, num_runs)

    all_sd_bias_within     = zeros(T, num_runs)
    all_sd_lambda_within   = zeros(T, num_runs)
    all_sd_war_within      = zeros(T, num_runs)
    all_sd_raid_within     = zeros(T, num_runs)
    all_sd_trade_within    = zeros(T, num_runs)
    all_sd_inaction_within = zeros(T, num_runs)

    for (k, h) in enumerate(hist)
        all_bias[:, k]      .= h.bias
        all_lambda[:, k]    .= h.lambda
        all_war[:, k]       .= h.war
        all_raid[:, k]      .= h.raid
        all_trade[:, k]     .= h.trade
        all_inaction[:, k]  .= h.inaction

        all_sd_bias_within[:, k]     .= h.sd_bias_within
        all_sd_lambda_within[:, k]   .= h.sd_lambda_within
        all_sd_war_within[:, k]      .= h.sd_war_within
        all_sd_raid_within[:, k]     .= h.sd_raid_within
        all_sd_trade_within[:, k]    .= h.sd_trade_within
        all_sd_inaction_within[:, k] .= h.sd_inaction_within
    end

    mean_bias_over_runs       = mean(all_bias, dims=2)
    sd_bias_over_runs         = std(all_bias, dims=2)
    mean_lambda_over_runs     = mean(all_lambda, dims=2)
    sd_lambda_over_runs       = std(all_lambda, dims=2)
    mean_war_over_runs        = mean(all_war, dims=2)
    sd_war_over_runs          = std(all_war, dims=2)
    mean_raid_over_runs       = mean(all_raid, dims=2)
    sd_raid_over_runs         = std(all_raid, dims=2)
    mean_trade_over_runs      = mean(all_trade, dims=2)
    sd_trade_over_runs        = std(all_trade, dims=2)
    mean_inaction_over_runs   = mean(all_inaction, dims=2)
    sd_inaction_over_runs     = std(all_inaction, dims=2)

    mean_sd_bias_within_over_runs     = mean(all_sd_bias_within, dims=2)
    mean_sd_lambda_within_over_runs   = mean(all_sd_lambda_within, dims=2)
    mean_sd_war_within_over_runs      = mean(all_sd_war_within, dims=2)
    mean_sd_raid_within_over_runs     = mean(all_sd_raid_within, dims=2)
    mean_sd_trade_within_over_runs    = mean(all_sd_trade_within, dims=2)
    mean_sd_inaction_within_over_runs = mean(all_sd_inaction_within, dims=2)

    return all_bias, vec(mean_bias_over_runs), vec(sd_bias_over_runs),
           all_lambda, vec(mean_lambda_over_runs), vec(sd_lambda_over_runs),
           all_war, vec(mean_war_over_runs), vec(sd_war_over_runs),
           all_raid, vec(mean_raid_over_runs), vec(sd_raid_over_runs),
           all_trade, vec(mean_trade_over_runs), vec(sd_trade_over_runs),
           all_inaction, vec(mean_inaction_over_runs), vec(sd_inaction_over_runs),
           vec(mean_sd_bias_within_over_runs),
           vec(mean_sd_lambda_within_over_runs),
           vec(mean_sd_war_within_over_runs),
           vec(mean_sd_raid_within_over_runs),
           vec(mean_sd_trade_within_over_runs),
           vec(mean_sd_inaction_within_over_runs)
end

const DEFAULTS = (G=10_000, T=200_000, Cw=8.0, variable_sim=1,
                  mu_b=0.01, mu_h=0.01, sigma_b=0.02, sigma_h=0.02, num_runs=20)

G = DEFAULTS.G
T = DEFAULTS.T
Cw = DEFAULTS.Cw
num_runs = DEFAULTS.num_runs
variable_sim = DEFAULTS.variable_sim
mu_b = DEFAULTS.mu_b
mu_h = DEFAULTS.mu_h
sigma_b = DEFAULTS.sigma_b
sigma_h = DEFAULTS.sigma_h

SCENARIOS_INPUT = [
    (Cw = 2.5,),
    (Cw = 8.0,),
    (Cw = 14.0,),
]

function run_scenarios(scenarios::AbstractVector{<:NamedTuple})
    results = Vector{Dict}(undef, length(scenarios))
    for (idx, sc) in enumerate(scenarios)
        cw = get(sc, :Cw, DEFAULTS.Cw)
        nr = DEFAULTS.num_runs

        all_bias, mean_over_runs, sd_over_runs,
        all_lambda, mean_lambda_over_runs, sd_lambda_over_runs,
        all_war, mean_war_over_runs, sd_war_over_runs,
        all_raid, mean_raid_over_runs, sd_raid_over_runs,
        all_trade, mean_trade_over_runs, sd_trade_over_runs,
        all_inaction, mean_inaction_over_runs, sd_inaction_over_runs,
        sd_bias_within_over_runs,
        sd_lambda_within_over_runs,
        sd_war_within_over_runs,
        sd_raid_within_over_runs,
        sd_trade_within_over_runs,
        sd_inaction_within_over_runs =
            run_multiple(nr; G=DEFAULTS.G,
                            T=DEFAULTS.T,
                            Cw=cw,
                            variable_sim=DEFAULTS.variable_sim,
                            mu_b=DEFAULTS.mu_b,
                            mu_h=DEFAULTS.mu_h,
                            sigma_b=DEFAULTS.sigma_b,
                            sigma_h=DEFAULTS.sigma_h)

        results[idx] = Dict(
            :G => DEFAULTS.G,
            :T => DEFAULTS.T,
            :Cw => cw,
            :num_runs => nr,
            :all_bias => all_bias,
            :mean => mean_over_runs,
            :std => sd_over_runs,
            :all_lambda => all_lambda,
            :mean_lambda => mean_lambda_over_runs,
            :std_lambda => sd_lambda_over_runs,
            :all_war => all_war,
            :mean_war => mean_war_over_runs,
            :std_war => sd_war_over_runs,
            :all_raid => all_raid,
            :mean_raid => mean_raid_over_runs,
            :std_raid => sd_raid_over_runs,
            :all_trade => all_trade,
            :mean_trade => mean_trade_over_runs,
            :std_trade => sd_trade_over_runs,
            :all_inaction => all_inaction,
            :mean_inaction => mean_inaction_over_runs,
            :std_inaction => sd_inaction_over_runs,
            :sd_bias_within => sd_bias_within_over_runs,
            :sd_lambda_within => sd_lambda_within_over_runs,
            :sd_war_within => sd_war_within_over_runs,
            :sd_raid_within => sd_raid_within_over_runs,
            :sd_trade_within => sd_trade_within_over_runs,
            :sd_inaction_within => sd_inaction_within_over_runs
        )
    end
    return results
end

function scenario_to_long_df(r::Dict)
    T = length(r[:mean])
    return DataFrame(
        G                = fill(r[:G], T),
        T_total          = fill(r[:T], T),
        Cw               = fill(r[:Cw], T),
        num_runs         = fill(r[:num_runs], T),
        t                = collect(1:T),
        mean_bias        = Float32.(r[:mean]),
        sd_bias_runs     = Float32.(r[:std]),
        sd_bias_within   = Float32.(r[:sd_bias_within]),
        mean_lambda      = Float32.(r[:mean_lambda]),
        sd_lambda_runs   = Float32.(r[:std_lambda]),
        sd_lambda_within = Float32.(r[:sd_lambda_within]),
        mean_war         = Float32.(r[:mean_war]),
        sd_war_runs      = Float32.(r[:std_war]),
        sd_war_within    = Float32.(r[:sd_war_within]),
        mean_raid        = Float32.(r[:mean_raid]),
        sd_raid_runs     = Float32.(r[:std_raid]),
        sd_raid_within   = Float32.(r[:sd_raid_within]),
        mean_trade       = Float32.(r[:mean_trade]),
        sd_trade_runs    = Float32.(r[:std_trade]),
        sd_trade_within  = Float32.(r[:sd_trade_within]),
        mean_inaction    = Float32.(r[:mean_inaction]),
        sd_inaction_runs = Float32.(r[:std_inaction]),
        sd_inaction_within = Float32.(r[:sd_inaction_within]),
    )
end

function write_results_summary_csv(results::Vector{Dict}, filepath::AbstractString)
    if isfile(filepath)
        rm(filepath)
    end
    first = true
    for r in results
        df = scenario_to_long_df(r)
        if first
            CSV.write(filepath, df)
            first = false
        else
            CSV.write(filepath, df; append = true, writeheader = false)
        end
    end
    return filepath
end

scenario_results = run_scenarios(SCENARIOS_INPUT)
csv_path = "Minimal_law_rich_poor_nonavoid.csv"
write_results_summary_csv(scenario_results, csv_path)
println("Wrote summarized results to: $csv_path")
println("Simulation complete.")
