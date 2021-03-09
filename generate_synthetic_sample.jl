using DelimitedFiles
using Random
using Distributions
using EmpiricalCDFs



function main()
	mu, population_size, seed_normal, seed_error, errorrate, n_samples, contamination, clonal_peak, seq_depth, n_sites, ploidy, outputdir, treefile, cdffile, seqdepthdistr_type, logerr = read_args()


	Random.seed!(seed_normal)
	rng_error = MersenneTwister(seed_error)

	seqdepthdistr = EmpiricalCDF()
	if seqdepthdistr_type == "uniform"
		read_empiricalcdf(seq_depth, seqdepthdistr)
	elseif seqdepthdistr_type == "empirical"
		read_empiricalcdf(cdffile, seqdepthdistr)
	end

	branchlengths = Dict{Rational, Float64}()
	
	read_tree(treefile, branchlengths, clonal_peak, mu)

	for nrun in 1:n_samples
		mutations, total_reads_hist = Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}(), Dict{Int64, Int64}()
	
		n_observed_all_mutations, n_observed_true_mutations = simulate_data(n_sites, branchlengths, mu, seqdepthdistr, errorrate, mutations, rng_error, total_reads_hist, logerr, population_size, ploidy, contamination)

		print_output(outputdir, mutations, treefile, nrun, n_observed_all_mutations, n_observed_true_mutations, total_reads_hist, mu, population_size, contamination, clonal_peak, n_sites, errorrate, cdffile, seed_normal, seed_error, seq_depth, seqdepthdistr_type, ploidy)
	end
end


function read_args()

	help_message = "\nusage:\njulia generate_synthetic_sample.jl [ARGUMENT LIST]\n\narguments:\n\n\t-m (float):\tmutation rate, per whole genome. MANDATORY PARAMETER. No default value.\n\n\t-t (string):\ttreefile. Should contain the access path. MANDATORY PARAMETER. No default value.\n\n\t-n (int):\tpopulation size. MANDATORY PARAMETER. No default value.\n\n\t-sqt (string):\ttype of sequencing depth distribution. MANDATORY PARAMETER. Possible values: uniform | empirical. No default value.\n\n\t-sn (int):\tseed of the RNG. Default: 1111\n\n\t-sr (int):\tseed of the RNG, dedicated for sequencing errors. Default 1112\n\n\t-e (float):\terror rate. Default 0.0\n\n\t-s (int):\tnumber of samples to be generated. Default 1\n\n\t-cont (float):\tcontamination of the sample by normal cells. Should be between 0.0 and 1.0. Default: 0.0\n\n\t-cp (int):\tnumber of clonal mutations. Default: 0\n\n\t-sd (int):\tsequencing depth, uniform for all sites. Default 100\n\n\t-l (int):\tlength of the genome, number of sites. Default 3 088 286 401\n\n\t-p (int):\tploidy. Default 2\n\n\t-odir (string): output directory. Default is the current directory.\n\n\t-cdf (string):\tfile containing the empirical sequencing depth distribution. Should contain the access path. No default value.\n\n\t-h, --help, -help, help:\tdisplay this help message\n"
	
	mu = -1.0
	population_size = 0
	seed_normal = 1111
	seed_error = 1112
	errorrate = 0.0
	n_samples = 1
	contamination = 0.0
	clonal_peak = 0
	seq_depth = 100
	n_sites = 3088286401
	ploidy = 2
	seqdepthdistr_type = ""
	outputdir = "./" # "gen2_err_e$(logerr)_treesize$(population_size)_beta$(beta)_treesDominik"
	cdffile = ""     # "/mnt/data/tibelyg/o2n_data/varscan/readcounts_corrected_tumor1-normal-v3.snp"
	treefile = ""    # "/mnt/data/tibelyg/geneclocks/$(population_size)_$(beta)_1_treesim/peaks_tree_n$(population_size)_sampled1_$(beta)_$(treeid).txt"

	if errorrate < 0.0000000000001
		logerr = 0
	else
		logerr = ceil(Int64, -log10(errorrate))
	end

	println(stderr, "")

	if length(ARGS) == 0  ||  ARGS[1] == "-h"  ||  ARGS[1] == "-help"  ||  ARGS[1] == "help"  ||  ARGS[1] == "--help"
		println(stdout, help_message)
		exit()
	else
		for arg_id in 1:2:length(ARGS)
			flag, value = ARGS[arg_id], ARGS[arg_id + 1]
			
			if flag == "-m"
				mu = parse(Float64, value)
			elseif flag == "-n"
				population_size = parse(Int64, value)
			elseif flag == "-sn"
				seed_normal = parse(Int64, value)
			elseif flag == "-sr"
				seed_error = parse(Int64, value)
			elseif flag == "-e"
				errorrate = parse(Float64, value)
			elseif flag == "-s"
				n_samples = parse(Int64, value)
			elseif flag == "-cont"
				contamination = parse(Float64, value)
			elseif flag == "-cp"
				clonal_peak = parse(Int64, value)
			elseif flag == "-sd"
				seq_depth = parse(Int64, value)
			elseif flag == "-l"
				n_sites = parse(Int64, value)
			elseif flag == "-p"
				ploidy = parse(Int64, value)
			elseif flag == "-odir"
				outputdir = value
			elseif flag == "-t"
				treefile = value
			elseif flag == "-cdf"
				cdffile = value
			elseif flag == "-sqt"
				seqdepthdistr_type = value
			end
		end
	end
	
	if mu < -0.1
		println(stderr, "mutation rate must be defined!")
		exit()
	end
	if treefile == ""
		println(stderr, "treefile must be defined!")
		exit()
	end
	if population_size == 0
		println(stderr, "population size must be defined!")
		exit()
	end
	if seqdepthdistr_type == ""
		println(stderr, "sequencing depth distribution type must be defined!")
		exit()
	end
	
	return mu, population_size, seed_normal, seed_error, errorrate, n_samples, contamination, clonal_peak, seq_depth, n_sites, ploidy, outputdir, treefile, cdffile, seqdepthdistr_type, logerr
end


function read_empiricalcdf(cdffile::String, seqdepthdistr::EmpiricalCDF{Float64})

	cdfdata = readdlm(cdffile, '\t', Int64)
	for i in 1:size(cdfdata, 1)
		for n in 1:cdfdata[i,5]
			push!(seqdepthdistr, cdfdata[i,1])
		end
	end
	sort!(seqdepthdistr)
end


function read_empiricalcdf(seq_depth::Int64, seqdepthdistr::EmpiricalCDF{Float64})

	push!(seqdepthdistr, seq_depth)
	sort!(seqdepthdistr)
end


function read_tree(treefile::String, branchlengths::Dict{Rational, Float64}, clonal_peak::Int64, mu::Float64)

	empty!(branchlengths)

	treedata = readdlm(treefile, '\t')
	for i in 1:size(treedata, 1)
		numerator, denominator = parse(Int64, treedata[i,1]), parse(Int64, treedata[i,2])
		branch = numerator // denominator
		if haskey(branchlengths, branch)
			branchlengths[branch] += treedata[i,3]
		else
			branchlengths[branch] = treedata[i,3]
		end
	end
	
	if clonal_peak > 0
		branchlengths[1 // 1] = clonal_peak / mu
	end
end


function simulate_data(n_sites::Int64, branchlengths::Dict{Rational, Float64}, mu::Float64, seqdepthdistr::EmpiricalCDF{Float64}, errorrate::Float64, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, rng_error::MersenneTwister, total_reads_hist::Dict{Int64, Int64}, logerr::Int64, population_size::Int64, ploidy::Int64, contamination::Float64)

	empty!(mutations)

	mu_sitewise = mu / n_sites
    n_observed_all_mutations, n_observed_true_mutations = 0, 0
    errorrate_third = errorrate / 3.0

	for site in 1:n_sites
		total_reads = convert(Int64, rand(seqdepthdistr))
		if haskey(total_reads_hist, total_reads)
			total_reads_hist[total_reads] += 1
		else
			total_reads_hist[total_reads]  = 1
		end
		
		site_is_mutant_flag = false
		selected_branch = 0//1
		for branch in keys(branchlengths)
			if rand() > exp(-mu_sitewise * branchlengths[branch])
				site_is_mutant_flag = true
				selected_branch = branch
			end
		end
		
		if site_is_mutant_flag
			successes = round(Int64, population_size * selected_branch)		# there are population_size tumour cells in the sequenced sample
			failures = round(Int64, population_size * ploidy / (1.0 - contamination)) - successes
			true_mutant_reads = rand(Hypergeometric(successes, failures, total_reads))
		else
			true_mutant_reads = 0
		end

		reads_overall = Int64[true_mutant_reads, 0, 0, total_reads - true_mutant_reads]
		observed_mutant_reads = true_mutant_reads
		if logerr > 0
			multi = Multinomial(total_reads - true_mutant_reads, [errorrate_third, errorrate_third, errorrate_third, 1.0 - errorrate])
			reads_fromwildtype = rand(rng_error, multi)		# [real_mutant_allele, false_mutant_different_allele, false_mutant_different_allele, wild_type_allele]

			multi = Multinomial(true_mutant_reads, [1.0 - errorrate, errorrate_third, errorrate_third, errorrate_third])
			reads_frommutant = rand(rng_error, multi)		# [real_mutant_allele, false_mutant_different_allele, false_mutant_different_allele, wild_type_allele]

			reads_overall = reads_fromwildtype + reads_frommutant

			observed_mutant_reads = reads_overall[1] + reads_overall[2] + reads_overall[3]
		end
		
		
		if observed_mutant_reads > 0
			mut_reads1, mut_reads2, mut_reads3 = sort(reads_overall[1:3])

			if haskey(mutations, total_reads)
				if haskey(mutations[total_reads], (mut_reads1, mut_reads2, mut_reads3))
					mutations[total_reads][mut_reads1, mut_reads2, mut_reads3] += 1
				else
					mutations[total_reads][mut_reads1, mut_reads2, mut_reads3] = 1
				end
			else
				mutations[total_reads] = Dict{Tuple{Int64, Int64, Int64}, Int64}()
				mutations[total_reads][mut_reads1, mut_reads2, mut_reads3] = 1
			end
			
			n_observed_all_mutations += 1
			if site_is_mutant_flag
				n_observed_true_mutations += 1
			end
		end
	end
	
	return n_observed_all_mutations, n_observed_true_mutations
end


function print_output(outputdir::String, mutations::Dict{Int64, Dict{Tuple{Int64, Int64, Int64}, Int64}}, treefile::String, nrun::Int64, n_observed_all_mutations::Int64, n_observed_true_mutations::Int64, total_reads_hist::Dict{Int64, Int64}, mu::Float64, population_size::Int64, contamination::Float64, clonal_peak::Int64, n_sites::Int64, errorrate::Float64, cdffile::String, seed_normal::Int64, seed_error::Int64, seq_depth::Int64, seqdepthdistr_type::String, ploidy::Int64)

	if !isdir(outputdir)
		mkdir(outputdir)
	end

	println(stderr, "$n_observed_all_mutations observed mutations\n$n_observed_true_mutations observed true mutations")

	outfile = "$outputdir/syntheticdata_run$nrun.txt"
	fout = open(outfile, "w")
	for total_reads in sort(collect(keys(mutations)))
		for (observed_mutant_reads1, observed_mutant_reads2, observed_mutant_reads3) in sort(collect(keys(mutations[total_reads])))
			println(fout, "$total_reads\t$observed_mutant_reads1\t$observed_mutant_reads2\t$observed_mutant_reads3\t$(mutations[total_reads][observed_mutant_reads1, observed_mutant_reads2, observed_mutant_reads3])")
		end
	end
	close(fout)
	
	outfile = "$outputdir/total_reads_hist_run$nrun.txt"
	fout = open(outfile, "w")
	for total_reads in sort(collect(keys(total_reads_hist)))
		println(fout, "$total_reads\t$(total_reads_hist[total_reads])")
	end
	close(fout)

	outfile = "$outputdir/parameters_run$nrun.txt"
	fout = open(outfile, "w")
	println(fout, "mu $mu\npopulation size $population_size\ncontamination $contamination\nclonal peak $clonal_peak\ntreefile $treefile\ndna length $n_sites\nsequencing error rate $errorrate\nseqdepth distribution type $seqdepthdistr_type")
	if seqdepthdistr_type == "uniform"
		println(fout, "uniform seq depth $seq_depth")
	else
		println(fout, "seq depth distribution file $cdffile")
	end
	println(fout, "ploidy $ploidy\nseed_normal $seed_normal\nseed_error $seed_error\n$n_observed_all_mutations observed mutations\n$n_observed_true_mutations observed true mutations")
	close(fout)
end


main()


