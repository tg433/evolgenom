using DelimitedFiles
using SpecialFunctions
using Distributions
using Statistics
using Optim



function main()

	start_time = time_ns()
	
	n_trees, mutationfile, readdepthfile, optim_type, errorrate_const, targetdir, treedir, n_clonal_mutations_const, contamination, min_clonal_peak, max_clonal_peak, min_errorrate, max_errorrate, ploidy = read_args()
	
	branchlengths = Array{Dict{Rational, Float64}, 1}()
	read_trees(branchlengths, n_trees, treedir)

	elapsed_time = time_ns() - start_time
	println(stderr, "trees read\n$(Int(elapsed_time)/10^9) sec elapsed")

	start_time = time_ns()
	println(stderr, "processing $(mutationfile)...")

	mutations, log_fact = Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}(), Float64[]
	read_mutations(mutations, mutationfile, targetdir)
	
	ratio_of_total_reads, wildtypeonly_total_reads = Dict{Int64, Float64}(), Dict{Int64, Int64}()
	n_sites = read_total_readcount_distr(ratio_of_total_reads, wildtypeonly_total_reads, log_fact, mutations, readdepthfile, targetdir)

	loglikelihoods, mus = Array{Float64, 1}(undef, n_trees), Array{Float64, 1}(undef, n_trees)
	if optim_type == "error"
		n_clonal_mutations = n_clonal_mutations_const
		
		result = optimize(x -> optimize_error(x, mutations, branchlengths, n_trees, loglikelihoods, log_fact, ploidy, contamination, n_clonal_mutations, n_sites, ratio_of_total_reads, wildtypeonly_total_reads), min_errorrate, max_errorrate)
	
		errorrate = Optim.minimizer(result)
		
	elseif optim_type == "clonalpeak"
		errorrate = errorrate_const
		
		result = optimize(x -> optimize_clonalpeak(x, mutations, branchlengths, n_trees, loglikelihoods, errorrate, log_fact, ploidy, contamination, n_sites, ratio_of_total_reads, wildtypeonly_total_reads), min_clonal_peak, max_clonal_peak)
		
		n_clonal_mutations = round(Int64, Optim.minimizer(result))
		
	elseif optim_type == "off"
		n_clonal_mutations, errorrate = n_clonal_mutations_const, errorrate_const
		
	else
		error("possible values of optim_type are:   error | clonalpeak | off")
		
	end
	
	avg_loglikelihood, avg_mu, std_mu = evaluate_average_likelihood(mutations, branchlengths, n_trees, mus, loglikelihoods, errorrate, log_fact, ploidy, contamination, n_clonal_mutations, n_sites, ratio_of_total_reads, wildtypeonly_total_reads)

	print_results(targetdir*"estimations_"*mutationfile, targetdir*"treewise_estimations_"*mutationfile, avg_mu, std_mu, errorrate, n_clonal_mutations, avg_loglikelihood, loglikelihoods, mus, start_time, contamination, treedir, optim_type, min_errorrate, max_errorrate, min_clonal_peak, max_clonal_peak, n_trees, ploidy, mutationfile, readdepthfile)
end


function read_args()

	help_message = "\nusage: julia estimate_parameters.jl [ARGUMENT LIST]\n\narguments:\n\n\t-n (int):\tnumber of trees used. MANDATORY PARAMETER. No default value.\n\n\t-im (string):\tfile containing the mutant site readcounts. MANDATORY PARAMETER. No default value.\n\n\t-ir (string):\tfile containing the read depths of all sites. MANDATORY PARAMETER. No default value.\n\n\t-o (string):\tfor what parameter should the estimation optimize? Possible values are error, clonalpeak, off. Default: off\n\n\t-e (float):\terror rate. Should be between 0.0 and 1.0. Default: 0.0\n\n\t-dtarget (string):\ttarget directory, location of the input and output files. Default: current directory\n\n\t-dtree (string):\tdirectory of the trees. Default: current directory\n\n\t-cp (int):\tnumber of clonal mutations. Default: 0\n\n\t-cont (float):\tcontamination of the sample by normal cells. Should be between 0.0 and 1.0. Default: 0.0\n\n\t-minc (int):\tminimum number of clonal mutations. Only used if the number of clonal mutations is optimized. Default: 0\n\n\t-maxc (int):\tmaximum number of clonal mutations. Only used if the number of clonal mutations is optimized. Default: 1000\n\n\t-mine (float):\tminimum value of the error rate. Only used if the error rate is optimized. Default: 10.0^(-10)\n\n\t-maxe (float):\tmaximum value of the error rate. Only used if the error rate is optimized. Default: 0.01\n\n\t-p (int):\tploidy. Default: 2\n\n\t-h, --help, -help, help:\tdisplay this help message\n"

	n_trees = 0
	mutationfile = ""
	readdepthfile = ""
	optim_type = "off"
	errorrate_const = 0.0
	targetdir = "./"
	treedir = "./"
	n_clonal_mutations_const = 0
	contamination = 0.0
	min_clonal_peak = 0
	max_clonal_peak = 1000
	min_errorrate = 10.0^(-10) 
	max_errorrate = 0.01
	ploidy = 2
	
	println(stderr, "")

	if length(ARGS) == 0  ||  ARGS[1] == "-h"  ||  ARGS[1] == "-help"  ||  ARGS[1] == "help"  ||  ARGS[1] == "--help"
		println(stdout, help_message)
		exit()
	else
		for arg_id in 1:2:length(ARGS)
			flag, value = ARGS[arg_id], ARGS[arg_id + 1]
			
			if flag == "-n"
				n_trees = parse(Int64, value)
			elseif flag == "-im"
				mutationfile = value
			elseif flag == "-ir"
				readdepthfile = value
			elseif flag == "-o"
				optim_type = value
			elseif flag == "-e"
				errorrate_const = parse(Float64, value)
			elseif flag == "-dtarget"
				targetdir = value
			elseif flag == "-dtree"
				treedir = value
			elseif flag == "-cp"
				n_clonal_mutations_const = parse(Int64, value)
			elseif flag == "-cont"
				contamination = parse(Float64, value)
			elseif flag == "-minc"
				min_clonal_peak = parse(Int64, value)
			elseif flag == "-maxc"
				max_clonal_peak = parse(Int64, value)
			elseif flag == "-mine"
				min_errorrate = parse(Float64, value)
			elseif flag == "-maxe"
				max_errorrate = parse(Float64, value)
			elseif flag == "-p"
				ploidy = parse(Int64, value)
			end
		end
	end

	if n_trees == 0
		println(stderr, "n_trees must be defined!")
		exit()
	end
	if mutationfile == ""
		println(stderr, "the filename containing the mutant site readcounts must be defined!")
		exit()
	end
	if readdepthfile == ""
		println(stderr, "the filename containing the readdepths must be defined!")
		exit()
	end
	
	if targetdir[end] != '/'
		targetdir *= "/"
	end
	if treedir[end] != '/'
		treedir *= "/"
	end
	
	return n_trees, mutationfile, readdepthfile, optim_type, errorrate_const, targetdir, treedir, n_clonal_mutations_const, contamination, min_clonal_peak, max_clonal_peak, min_errorrate, max_errorrate, ploidy
end


function read_trees(branchlengths::Array{Dict{Rational, Float64}, 1}, n_trees::Int64, treedir::String)

	for treeid in 1:n_trees
		push!(branchlengths, Dict{Rational, Float64}())
		
		treefile = treedir*"peaks_tree_n10000_sampled1_0.99_$(treeid).txt"
# 		treefile = treedir*"tree_$(treeid).txt"

		fin = open(treefile, "r")
		lines = readlines(fin)
		close(fin)

		for line in lines
			a, b, c = split(line, '\t')
			numerator, denominator, reallength = parse(Int64, a), parse(Int64, b), parse(Float64, c)
			branch = numerator // denominator
			if haskey(branchlengths[treeid], branch)
				branchlengths[treeid][branch] += reallength
			else
				branchlengths[treeid][branch]  = reallength
			end
		end
	end
end


function read_mutations(mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, mutationfile::String, targetdir::String)

	empty!(mutations)

	readcountdata = readdlm(targetdir*mutationfile, '\t', Int64)
	for i in 1:size(readcountdata, 1)
		total_reads, n_mutant_sites = readcountdata[i,1], readcountdata[i,5]
		mutant_reads1, mutant_reads2, mutant_reads3 = sort(readcountdata[i,2:4])
		if n_mutant_sites > 0  &&  mutant_reads3 > 0
			if haskey(mutations, total_reads)
				if haskey(mutations[total_reads], (mutant_reads1, mutant_reads2, mutant_reads3))
					mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3] += n_mutant_sites
				else
					mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3]  = n_mutant_sites
				end
			else
				mutations[total_reads] = Dict{Tuple{Int64, Int64, Int64}, Int64}((mutant_reads1, mutant_reads2, mutant_reads3) => n_mutant_sites)
			end
		end
	end
end


function read_total_readcount_distr(ratio_of_total_reads::Dict{Int64, Float64}, wildtypeonly_total_reads::Dict{Int64, Int64}, log_fact::Array{Float64, 1}, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, readdepthfile::String, targetdir::String)

	empty!(ratio_of_total_reads)
	empty!(log_fact)
	
	total_readcount_data = readdlm(targetdir*readdepthfile, '\t', Int64)

	sum_sites_in_total_reads = sum(total_readcount_data[:,2])
	for i in 1:size(total_readcount_data, 1)
		ratio_of_total_reads[total_readcount_data[i,1]] = total_readcount_data[i,2] / sum_sites_in_total_reads
		
		if !haskey(mutations, total_readcount_data[i,1])
			wildtypeonly_total_reads[total_readcount_data[i,1]] = total_readcount_data[i,2]
		end
	end
	
	for i in 1:maximum(keys(ratio_of_total_reads))
		push!(log_fact, lfactorial(i))
	end
	
	doesitsumtoone = sum(values(ratio_of_total_reads))
	if doesitsumtoone <= 0.999999999999  ||  doesitsumtoone >= 1.000000000001
		error("ratio_of_total_reads does not sum to 1! ($doesitsumtoone)")
	end
	
	return sum_sites_in_total_reads
end


function optimize_error(errorrate::Float64, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, branchlengths::Array{Dict{Rational, Float64}, 1}, n_trees::Int64, loglikelihoods::Array{Float64, 1}, log_fact::Array{Float64, 1}, ploidy::Int64, contamination::Float64, n_clonal_mutations::Int64, n_sites::Int64, ratio_of_total_reads::Dict{Int64, Float64}, wildtypeonly_total_reads::Dict{Int64, Int64})

	if 0 <= n_clonal_mutations <= n_sites
		log_p_mutant, log_p_wildtype, log_branchrelsize = Float64[], Float64[], Float64[]
		freq_threshold_error = Dict{Int64, Float64}()
		clonalpeakposition = 1.0 / ploidy * (1.0 - contamination)

		n_nonclonal_mutations = count_n_nonclonal_mutations_and_load_freq_threshold_error(freq_threshold_error, ratio_of_total_reads, n_sites, errorrate, mutations, n_clonal_mutations)	
		
		if n_nonclonal_mutations > 0
			for treeid in 1:n_trees
				mu, L = load_tree(log_p_mutant, log_p_wildtype, log_branchrelsize, branchlengths[treeid], errorrate, clonalpeakposition, ploidy, contamination, n_clonal_mutations, n_nonclonal_mutations, freq_threshold_error, ratio_of_total_reads)
				
				loglikelihoods[treeid] = calculate_likelihood(mutations, log_p_mutant, log_p_wildtype, log_branchrelsize, errorrate, log_fact, mu, L, n_sites, ratio_of_total_reads, wildtypeonly_total_reads)
			end

			return -Statistics.mean(loglikelihoods)			# minus is because the optimization algorithm looks for the minimum of the loglikelihood, and we need the maximum
		else
			return Inf
		end
	else
		return Inf
	end
end


function optimize_clonalpeak(n_clonal_mutations_float::Float64, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, branchlengths::Array{Dict{Rational, Float64}, 1}, n_trees::Int64, loglikelihoods::Array{Float64, 1}, errorrate::Float64, log_fact::Array{Float64, 1}, ploidy::Int64, contamination::Float64, n_sites::Int64, ratio_of_total_reads::Dict{Int64, Float64}, wildtypeonly_total_reads::Dict{Int64, Int64})

	n_clonal_mutations = round(Int64, n_clonal_mutations_float)

	if 0 <= n_clonal_mutations <= n_sites
		log_p_mutant, log_p_wildtype, log_branchrelsize = Float64[], Float64[], Float64[]
		freq_threshold_error = Dict{Int64, Float64}()
		clonalpeakposition = 1.0 / ploidy * (1.0 - contamination)

		n_nonclonal_mutations = count_n_nonclonal_mutations_and_load_freq_threshold_error(freq_threshold_error, ratio_of_total_reads, n_sites, errorrate, mutations, n_clonal_mutations)	
		
		if n_nonclonal_mutations > 0
			for treeid in 1:n_trees
				mu, L = load_tree(log_p_mutant, log_p_wildtype, log_branchrelsize, branchlengths[treeid], errorrate, clonalpeakposition, ploidy, contamination, n_clonal_mutations, n_nonclonal_mutations, freq_threshold_error, ratio_of_total_reads)
				
				loglikelihoods[treeid] = calculate_likelihood(mutations, log_p_mutant, log_p_wildtype, log_branchrelsize, errorrate, log_fact, mu, L, n_sites, ratio_of_total_reads, wildtypeonly_total_reads)
			end

			return -Statistics.mean(loglikelihoods)			# minus is because the optimization algorithm looks for the minimum of the loglikelihood, and we need the maximum
		else
			return Inf
		end
	else
		return Inf
	end
end


function evaluate_average_likelihood(mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, branchlengths::Array{Dict{Rational, Float64}, 1}, n_trees::Int64, mus::Array{Float64, 1}, loglikelihoods::Array{Float64, 1}, errorrate::Float64, log_fact::Array{Float64, 1}, ploidy::Int64, contamination::Float64, n_clonal_mutations::Int64, n_sites::Int64, ratio_of_total_reads::Dict{Int64, Float64}, wildtypeonly_total_reads::Dict{Int64, Int64})

	log_p_mutant, log_p_wildtype, log_branchrelsize = Float64[], Float64[], Float64[]
	freq_threshold_error = Dict{Int64, Float64}()
	clonalpeakposition = 1.0 / ploidy * (1.0 - contamination)

	n_nonclonal_mutations = count_n_nonclonal_mutations_and_load_freq_threshold_error(freq_threshold_error, ratio_of_total_reads, n_sites, errorrate, mutations, n_clonal_mutations)	
	
	if n_nonclonal_mutations > 0
		for treeid in 1:n_trees
			mus[treeid], L = load_tree(log_p_mutant, log_p_wildtype, log_branchrelsize, branchlengths[treeid], errorrate, clonalpeakposition, ploidy, contamination, n_clonal_mutations, n_nonclonal_mutations, freq_threshold_error, ratio_of_total_reads)
			
			loglikelihoods[treeid] = calculate_likelihood(mutations, log_p_mutant, log_p_wildtype, log_branchrelsize, errorrate, log_fact, mus[treeid], L, n_sites, ratio_of_total_reads, wildtypeonly_total_reads)
		end

		return Statistics.mean(loglikelihoods), Statistics.mean(mus), Statistics.std(mus)
	else
		return Inf, 0.0, 0.0
	end
end


function count_n_nonclonal_mutations_and_load_freq_threshold_error(freq_threshold_error::Dict{Int64, Float64}, ratio_of_total_reads::Dict{Int64, Float64}, n_sites::Int64, errorrate::Float64, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, n_clonal_mutations::Int64)

		# count the nonclonal mutations which are not affected by errors
	
	n_nonclonal_mutations = 0
	
	for total_reads in keys(ratio_of_total_reads)
		for false_mutant_reads in 1:total_reads
			expected_false_mutant_sites = n_sites * pdf(Binomial(total_reads, errorrate), false_mutant_reads)
			if expected_false_mutant_sites < 0.999999999999  &&  false_mutant_reads >= total_reads * errorrate		# second term is to get after the maximum of the expected_false_mutant_sites(false_mutant_reads) function
				freq_threshold_error[total_reads] = false_mutant_reads / total_reads
				break
			end
		end
		if !haskey(freq_threshold_error, total_reads)
			freq_threshold_error[total_reads] = 2.0
		end

		if haskey(mutations, total_reads)
			for (mutant_reads1, mutant_reads2, mutant_reads3) in keys(mutations[total_reads])
				if mutant_reads1 / total_reads >= freq_threshold_error[total_reads] - 0.000000000001
					n_nonclonal_mutations += mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3]
				end
				if mutant_reads2 / total_reads >= freq_threshold_error[total_reads] - 0.000000000001
					n_nonclonal_mutations += mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3]
				end
				if mutant_reads3 / total_reads >= freq_threshold_error[total_reads] - 0.000000000001
					n_nonclonal_mutations += mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3]
				end
			end
		end
	end
	n_nonclonal_mutations -= n_clonal_mutations
	
	return n_nonclonal_mutations
end


function load_tree(log_p_mutant::Array{Float64, 1}, log_p_wildtype::Array{Float64, 1}, log_branchrelsize::Array{Float64, 1}, branchlengths_current_tree::Dict{Rational, Float64}, errorrate::Float64, clonalpeakposition::Float64, ploidy::Int64, contamination::Float64, n_clonal_mutations::Int64, n_nonclonal_mutations::Int64, freq_threshold_error::Dict{Int64, Float64}, ratio_of_total_reads::Dict{Int64, Float64})

	empty!(log_p_mutant)
	empty!(log_p_wildtype)
	empty!(log_branchrelsize)
	
	errorrate_third, errorrate_fourthird = errorrate / 3.0, errorrate * 4.0 / 3.0

	push!(log_p_mutant, log(clonalpeakposition + errorrate_third - clonalpeakposition * errorrate_fourthird))
	push!(log_p_wildtype, log(1.0 - (clonalpeakposition + errorrate - clonalpeakposition * errorrate_fourthird)))
	
	branchlist = sort(collect(keys(branchlengths_current_tree)))
	reduced_branchlengthsum = 0.0

	for branch in branchlist
		branch_corrected = branch / ploidy * (1.0 - contamination)
		push!(log_p_mutant, log(branch_corrected + errorrate_third - branch_corrected * errorrate_fourthird))
		push!(log_p_wildtype, log(1.0 - (branch_corrected + errorrate - branch_corrected * errorrate_fourthird)))
		
		for total_reads in keys(ratio_of_total_reads)
			max_limit = max(1, round(Int64, freq_threshold_error[total_reads] * total_reads)) - 1
			reduced_branchlengthsum += branchlengths_current_tree[branch]  *  ccdf(Binomial(total_reads, branch_corrected), max_limit)  *  ratio_of_total_reads[total_reads]
		end
	end
	
	mu = n_nonclonal_mutations / reduced_branchlengthsum			# this is the mutation rate of ALL SITES COMBINED


		# 	clonalpeak_branchlength has the dimension of number_of_mutations / mutation_rate. The clonal branchlength which generates n_clonal_mutations mutations with mutation rate mu is n_clonal_mutations / mu; and there is the additional requirement that it is the number of mutations having mutant_reads > 0 that is observed.
	
	clonalpeak_branchlength = 0.0
	for total_reads in keys(ratio_of_total_reads)
		clonalpeak_branchlength += n_clonal_mutations / mu / ccdf(Binomial(total_reads, clonalpeakposition), 0) * ratio_of_total_reads[total_reads]			# ccdf is due to the fact that n_clonal_mutations is the number of visible mutations, using a certain sequencing depth
	end
	
	sum_branchlengths = sum(values(branchlengths_current_tree)) + clonalpeak_branchlength
	
	push!(log_branchrelsize, log(clonalpeak_branchlength / sum_branchlengths))
	for branch in branchlist										# it is important to go over branches in exactly the same order as when filling up log_p_mutant and log_p_wildtype
		push!(log_branchrelsize, log(branchlengths_current_tree[branch] / sum_branchlengths))
	end
	
	return mu, sum_branchlengths
end


function calculate_likelihood(mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, log_p_mutant::Array{Float64, 1}, log_p_wildtype::Array{Float64, 1}, log_branchrelsize::Array{Float64, 1}, errorrate::Float64, log_fact::Array{Float64, 1}, mu::Float64, L::Float64, n_sites::Int64, ratio_of_total_reads::Dict{Int64, Float64}, wildtypeonly_total_reads::Dict{Int64, Int64})

	loglikelihood = 0.0
	log_errorrate_third, log_1merrorrate, log_onethird = log(errorrate / 3.0), log(1.0 - errorrate), log(1.0 / 3.0)
	minus_mu_dot_branchlengthsum_per_site = -mu * L / n_sites
	one_minus_exp_minus_mu_dot_branchlengthsum_per_site = 1.0 - exp(minus_mu_dot_branchlengthsum_per_site)
	
	for total_reads in keys(mutations)
		for (mutant_reads1, mutant_reads2, mutant_reads3) in keys(mutations[total_reads])
			p, wildtype_reads = 0.0, total_reads - mutant_reads1 - mutant_reads2 - mutant_reads3
			combinatorial_term = log_fact[total_reads] - (log_fact_safe(mutant_reads1, log_fact) + log_fact_safe(mutant_reads2, log_fact) + log_fact[mutant_reads3] + log_fact[wildtype_reads])
			
			for k in 1:length(log_p_mutant)
				p += exp(log_onethird + log_branchrelsize[k] + combinatorial_term + mutant_reads1 * log_p_mutant[k] + (mutant_reads2 + mutant_reads3) * log_errorrate_third + wildtype_reads * log_p_wildtype[k]) + exp(log_onethird + log_branchrelsize[k] + combinatorial_term + mutant_reads2 * log_p_mutant[k] + (mutant_reads1 + mutant_reads3) * log_errorrate_third + wildtype_reads * log_p_wildtype[k])
				if mutant_reads1 + mutant_reads2 == 0				# the (mutant_reads1 + mutant_reads2) * log_errorrate_third might be 0 * -Inf, in which case we expect it to be 0
                    p += exp(log_onethird + log_branchrelsize[k] + combinatorial_term + mutant_reads3 * log_p_mutant[k] + wildtype_reads * log_p_wildtype[k])
				else
                    p += exp(log_onethird + log_branchrelsize[k] + combinatorial_term + mutant_reads3 * log_p_mutant[k] + (mutant_reads1 + mutant_reads2) * log_errorrate_third + wildtype_reads * log_p_wildtype[k])
                end
			end
			
			p *= one_minus_exp_minus_mu_dot_branchlengthsum_per_site
			
			p += exp(minus_mu_dot_branchlengthsum_per_site + combinatorial_term + (mutant_reads1 + mutant_reads2 + mutant_reads3) * log_errorrate_third + wildtype_reads * log_1merrorrate)
			
			loglikelihood += mutations[total_reads][mutant_reads1, mutant_reads2, mutant_reads3] * log(p)
		end
			# mutant_reads = 0 sites:
			
		p = 0.0
		for k in 1:length(log_p_wildtype)
			p += exp(log_branchrelsize[k] + total_reads * log_p_wildtype[k])
		end
		p *= one_minus_exp_minus_mu_dot_branchlengthsum_per_site
		p += exp(minus_mu_dot_branchlengthsum_per_site + total_reads * log_1merrorrate)
		
		loglikelihood += (ratio_of_total_reads[total_reads] * n_sites - sum(values(mutations[total_reads]))) * log(p)
	end

			# mutant_reads = 0 sites (continued):

	for total_reads in keys(wildtypeonly_total_reads)
		p = 0.0
		for k in 1:length(log_p_wildtype)
			p += exp(log_branchrelsize[k] + total_reads * log_p_wildtype[k])
		end
		p *= one_minus_exp_minus_mu_dot_branchlengthsum_per_site
		p += exp(minus_mu_dot_branchlengthsum_per_site + total_reads * log_1merrorrate)
		
		loglikelihood += wildtypeonly_total_reads[total_reads] * log(p)
	end
	
	return loglikelihood
end


function log_fact_safe(x::Int64, log_fact::Array{Float64, 1})

	if x == 0
		return 0.0
	else
		return log_fact[x]
	end
end


function print_results(outputfile_main::String, outputfile_treewise::String, mu::Float64, mu_std::Float64, errorrate::Float64, n_clonal_mutations::Int64, loglikelihood::Float64, loglikelihoods::Array{Float64, 1}, mus::Array{Float64, 1}, start_time::UInt64, contamination::Float64, treedir::String, optim_type::String, min_errorrate::Float64, max_errorrate::Float64, min_clonal_peak::Int64, max_clonal_peak::Int64, n_trees::Int64, ploidy::Int64, mutationfile::String, readdepthfile::String)

	fout = open(outputfile_treewise, "w")	
	for treeid in 1:n_trees
		println(fout, "$treeid\t$(loglikelihoods[treeid])\t$(mus[treeid])")
	end
	close(fout)

	fout = open(outputfile_main, "w")
	println(stderr, "")
	println(stderr, "mutationfile $mutationfile")
	println(fout,   "mutationfile $mutationfile")
	println(stderr, "readdepthfile $readdepthfile")
	println(fout,   "readdepthfile $readdepthfile")
	println(stderr, "tree directory $treedir")
	println(fout,   "tree directory $treedir")
	println(stderr, "mu $mu")
	println(fout,   "mu $mu")
	println(stderr, "mu std dev $mu_std")
	println(fout,   "mu std dev $mu_std")
	println(stderr, "errorrate $errorrate")
	println(fout,   "errorrate $errorrate")
	println(stderr, "clonal peak size $n_clonal_mutations")
	println(fout,   "clonal peak size $n_clonal_mutations")
	println(stderr, "contamination $contamination")
	println(fout,   "contamination $contamination")
	println(stderr, "ploidy $ploidy")
	println(fout,   "ploidy $ploidy")
	println(stderr, "optimization type $optim_type")
	println(fout,   "optimization type $optim_type")
	if optim_type == "error"
		println(stderr, "min error $min_errorrate, max error $max_errorrate")
		println(fout,   "min error $min_errorrate, max error $max_errorrate")
	elseif optim_type == "clonalpeak"
		println(stderr, "min clonal peak $min_clonal_peak, max clonal peak $max_clonal_peak")
		println(fout,   "min clonal peak $min_clonal_peak, max clonal peak $max_clonal_peak")
	end
	println(stderr, "$loglikelihood loglikelihood\n")
	println(fout,   "$loglikelihood loglikelihood\n")

	elapsed_time = time_ns() - start_time
	println(stderr, "$(Int(elapsed_time)/10^9) sec elapsed\n")
	println(fout,   "$(Int(elapsed_time)/10^9) sec elapsed\n")

	close(fout)

end


main()
