# Load relevant packages

using Random
using DataFrames
using CSV
using Plots

# Import the data

msc_data = DataFrame(CSV.File("path\\to\\HC.SitesDiffs.txt"));


# Log prior

function lnprior(tau::Float64,
                 theta::Float64,
                 time::Vector{Float64},
                 mu_tau::Float64,
                 mu_theta::Float64)

    prior_tau_theta = -(tau / mu_tau) - (theta / mu_theta)

    # .- and .* are element-wise operators
    prior_ti = sum( log(2 / theta) .- (2 / theta .* time) )
    
    return prior_tau_theta + prior_ti          

end 



# Log likelihood

function lnlikelihood(tau::Float64, time::Vector{Float64}, data::DataFrame) 

    p = 3/4 .- 3/4 * exp.(-8/3 * (tau .+ time)) # time is a vector

    lnp = sum( data[:, :xi] .* log.(p) .+ (data[:, :ni] .- data[:, :xi]) .* log.(1 .- p) )

    return lnp
    
end

# Lets break down the first line of the likelihood function

#= tau_test .+ time_vec: This performs element-wise addition between tau_test and each element of time_vec, resulting in a new array of the same size as time_vec.

-8/3 * (tau_test .+ time_vec): This multiplies each element of the new array from the previous step by -8/3.

exp.(-8/3 * (tau_test .+ time_vec)): This calculates the exponential of each element of the array from the previous step. We use exp. and not exp as we are performing this operation element-wise. 

3/4 * exp.(-8/3 * (tau_test .+ time_vec)): This multiplies each element of the array from the previous step by 3/4.

3/4 .- 3/4 * exp.(-8/3 * (tau_test .+ time_vec)): Finally, this subtracts each element of the array from the previous step from 3/4, resulting in the final array. =#



# Acceptance ratio

function ln_accept_ratio_time_j(j::Int64,
                                tau::Float64,
                                theta::Float64,
                                time_j::Float64,
                                new_time_j::Float64,
                                data::DataFrame)
  
    p_current = 3/4 - 3/4 * exp(-8/3 * (tau + time_j)) # current prob of data
    
    p_update = 3/4 - 3/4 * exp(-8/3 * (tau + new_time_j)) # prob of data with new t_j
    
    accept_ratio = -2 / theta * (new_time_j - time_j) + 
                            data[j, :xi] * log(p_update / p_current) + 
                            (data[j, :ni] - data[j, :xi]) * log((1 - p_update) / (1 - p_current))
    
    return accept_ratio

end



# MCMC algorithm 

function mcmc_MSC(;num_samples::Int64, 
                  tau::Float64, 
                  theta::Float64, 
                  mu_tau::Float64, 
                  mu_theta::Float64, 
                  time_0::Float64, 
                  w_tau::Float64, 
                  w_theta::Float64, 
                  w_time::Float64, 
                  nloci::Int64,
                  data::DataFrame)

        # Subset the data frame to have nrows = nloci

        data = data[1:nloci, :]

        coalescent_times = fill(time_0, nloci) # Starting value for t_i's
        sample_tau = [tau; zeros(num_samples)] # Array of length N+1 with 'tau' as the first value
        sample_theta = [theta; zeros(num_samples)]

        accept_tau = accept_theta = accept_time = 0

        tau_current = tau
        theta_current = theta

        # Initialise log prior and log likelihood

        lnp = lnprior(tau_current, theta_current, coalescent_times, mu_tau, mu_theta)  
        L = lnlikelihood(tau_current, coalescent_times, data)

        for i in 1:num_samples

            # Change tau

            tau_proposed = abs( tau_current + w_tau * (rand() - 0.5) )
            
            lnp_new = lnprior(tau_current, theta_current, coalescent_times, mu_tau, mu_theta)
            L_new = lnlikelihood(tau_proposed, coalescent_times, data)
            lnaccept = lnp_new + L_new - lnp - L 
        
            if lnaccept >= 0.0 || rand() < exp(lnaccept)

                tau_current = tau_proposed
                lnp = lnp_new
                L = L_new
                accept_tau += 1

            end

            # Change theta

            theta_proposed = abs( theta_current + w_theta * (rand() - 0.5) )

            lnp_new = lnprior(tau_current, theta_proposed, coalescent_times, mu_tau, mu_theta) # Only update the log prior because theta does not appear in the likelihood computation. 
            lnaccept = lnp_new - lnp

            if lnaccept >= 0.0 || rand() < exp(lnaccept)
                theta_current = theta_proposed
                lnp = lnp_new
                accept_theta += 1

            end
        
            # Update the j coalescent times

            for j in 1:nloci

                time_j_proposed = abs( coalescent_times[j] + w_time * (0.5 - rand()) )

                lnaccept = ln_accept_ratio_time_j(j, tau_current, theta_current, coalescent_times[j], time_j_proposed, data)

                if lnaccept >= 0.0 || rand() < exp(lnaccept)

                    coalescent_times[j] = time_j_proposed
                    accept_time += 1

                end

            end

        # Record the current parameter values

        sample_tau[i+1] = tau_current
        sample_theta[i+1] = theta_current

        # Update log prior and log likelihood

        lnp = lnprior(tau_current, theta_current, coalescent_times, mu_tau, mu_theta)   
        L = lnlikelihood(tau_current, coalescent_times, data)

        # Progress

        if i % (num_samples ÷ 20) == 0 

            println("$(round(Int64, i / num_samples * 100))% completed. Acceptance ratios:  tau $(accept_tau / i),  theta $(accept_theta / i),  time $(accept_time / (i * nloci))")

        end


    end

    # Calculate the proportion of accepted proposals

    accept_tau /= num_samples 
    accept_theta /= num_samples
    accept_time /= num_samples * nloci
        
    # Return results as a dictionary

    return Dict(
        "Tau samples" => sample_tau, 
        "Theta samples" => sample_theta, 
        "Acceptance proportion (Tau)" => accept_tau, 
        "Acceptance proportion (Theta)" => accept_theta, 
        "Acceptance proportion (Coalescent times)" => accept_time
    )

end



# Store results in a dictionary

result_dict = mcmc_MSC(num_samples=10000, tau=0.01, theta=0.001, mu_tau=0.005, mu_theta=0.001, time_0=0.001, w_tau=0.00075, w_theta=0.0012, w_time=0.02, nloci=1000, data = msc_data)

# Extract the samples

theta_samples = result_dict["Theta samples"]
tau_samples = result_dict["Tau samples"]

# Create the plots

p1 = plot(theta_samples, title = "Theta Samples", xlabel = "Iteration", ylabel = "Theta", legend = false, ylim = (0.001, 0.006))
p2 = plot(tau_samples, title = "Tau Samples", xlabel = "Iteration", ylabel = "Tau", legend = false, ylim = (0.003, 0.01))

# Place the plots in one window

plot(p1, p2, layout = (2, 1))
