using Plots, Random, LaTeXStrings, Statistics, Distributions, GLM, DataFrames, DelimitedFiles, DataStructures, LsqFit, StatsBase


function calculate_ci(data, confidence_level)
    n = length(data)
    mean_value = mean(data)
    std_error = std(data) / sqrt(n)

    z_value = quantile(Normal(), 1 - (1 - confidence_level) / 2)

    lower_bound = mean_value - z_value * std_error
    upper_bound = mean_value + z_value * std_error

    return lower_bound, upper_bound
end


#####################################################################################################

function creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

    for i in 1:Int(p_A*N)
        population_matrix[i,1]=1
        population_matrix[i,2]=1
    end

    for j in 1:Int(trunc((1-p_A)*N))
        population_matrix[j+Int(p_A*N),1]=2
        population_matrix[j+Int(p_A*N),2]=2
    end

    

    population_array_1=shuffle!(population_matrix[:,1])    # Shuffling of the columns of the matrix in order to create heterozygotes too
    population_array_2=shuffle!(population_matrix[:,2])

    population_matrix[:,1]=population_array_1
    population_matrix[:,2]=population_array_2

    return population_matrix
    
end

function next_generation(population_matrix, N)

    new_population=zeros(N,2) #rows, columns

    for j in 1:N

        for k in 1:2  # select one allele of each parent 

            selected_parent=rand(1:N) # randomized selection of aparent
            selected_allele=rand(1:2) # randomized selection of the allele 

            new_population[j,k]=population_matrix[selected_parent,selected_allele]
            
        end
        
    end 

    return new_population
    
end

function distributions_of_the_unfixed_allele_frequency(N, #population 
    p_A, #frequency for the allele A
    t_array,   #time of evolution of the population
    n_tries
)

    population_matrix=zeros(N,2) #array with one of the alleles of every individual
    


    plot(linewidth=3, xlabel="allele frequency", ylabel="probability density")

    for t_multiplier in t_array

        t=t_multiplier*N

        heterallelic_density_array=[]

        percentage_counter=0

        for _ in 1:n_tries
            
            creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

            # Loop over the time ###############################################################################

            time_array=[]
            counter_A=0

            for i in 1:t

                push!(time_array, i)

                new_population=next_generation(population_matrix, N)

                population_matrix=new_population

                counter_A=count(x -> x == 1, new_population)/(2*N)
            end


            if counter_A!=0.0 && counter_A!=1.0                 # not fixed populations
                percentage_counter+=2                           # allele frequency, not individual frequency
                push!(heterallelic_density_array, counter_A)
            end

            

        end

        percentage_counter=percentage_counter*100/n_tries
        
        count_dict = Dict()

        # Number of times each value appears
        for key in heterallelic_density_array
            count_dict[key] = get(count_dict, key, 0) + 1
        end

        sorted_keys = sort(collect(keys(count_dict)))

        ordered_dict = OrderedDict()
        for key in sorted_keys
            push!(ordered_dict, key => count_dict[key])
        end

        allele_frequency_array = collect(keys(ordered_dict))
        number_of_population_array = collect(values(ordered_dict)) .*percentage_counter ./ sum(values(ordered_dict))



        #histogram(heterallelic_density_array, title="t=$t_multiplier N", normalized=false, grid=false, label="", xlabel="allele frequency")
        plot!(allele_frequency_array, number_of_population_array, grid=false, label="t=$t_multiplier N")
        xlims!(0,1)

        savefig("Prob_distribution_N=$N-p=$p_A.png")

    end
    
    
    

end

# N/10, N/5, N/2, N, 2N, 3N


t_array_1=[1/10, 1/5, 1/2, 1, 2, 3]
t_array_2=[1/10, 1/5, 1/2, 1, 2, 4]

#distributions_of_the_unfixed_allele_frequency(100, 0.5, t_array_1, 10^7)
#distributions_of_the_unfixed_allele_frequency(100, 0.1, t_array_2, 10^7)


function evolution_plotting(N, #population 
    p_A, #frequency for the allele A
    t,   #time of evolution of the population
    n_tries)

    population_matrix=zeros(N,2) #array with one of the alleles of every individual

    time_array=[i for i in 1:1:t]

    first_momenT_matrix_fixed=zeros(t,0)
    second_momenT_matrix_fixed=zeros(t,0)

    first_moment_array=[]
    second_moment_array=[]
    variance_array=[]

    T_fixed=-4*N*(p_A*log(p_A)+(1-p_A)*log(1-p_A))
    

    plot(xlabel="t (Generations)", ylabel="x", grid=false, title="N=$N")

    for _ in 1:n_tries
        

        creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

        # Loop over the time ###############################################################################

        time_array=[]
        population_A=[]

        for i in 1:t

            push!(time_array, i)
            
            new_population=next_generation(population_matrix, N)

            population_matrix=new_population

            counter_A=count(x -> x == 1, new_population)/(2*N)

            push!(population_A, counter_A)
            
        end

        first_momenT_matrix_fixed=hcat(first_momenT_matrix_fixed, population_A)
        second_momenT_matrix_fixed=hcat(second_momenT_matrix_fixed, population_A.^2)

        plot!(time_array, population_A, label="")
        ylims!(0,1)
    
    end

    vline!([T_fixed], label=L"\bar{t}", color=:black)
    savefig("population_over_time_$N-$t-$n_tries.png")

    for i in 1:length(time_array)
        first=mean(first_momenT_matrix_fixed[i,:])
        second=mean(second_momenT_matrix_fixed[i,:])
        variance=second-first^2
        push!(first_moment_array, first)
        push!(second_moment_array, second)
        push!(variance_array, variance )
    end


    plot(time_array, first_moment_array, xlabel="t (Generations)", label=L"\left\langle x\left(t\right)\right\rangle", grid=false, title="N=$N")
    plot!(time_array, second_moment_array, label=L"\left\langle x^2\left(t\right)\right\rangle")
    plot!(time_array, variance_array, label=L"\sigma_{x}^{2}(t)")
    vline!([T_fixed], label=L"\bar{t}", color=:black)
    ylims!(-0.1, 1.1)
    savefig("moments_fig_$N-$t-$n_tries.png")

end

#evolution_plotting(50, 0.5, 500, 6)
#evolution_plotting(500, 0.5, 500, 6)
#evolution_plotting(500, 0.5, 5000, 6)

function moments_trajectories(N, p_A, t, n_tries)

    population_matrix=zeros(N,2) #array with one of the alleles of every individual

    time_array=[i for i in 1:1:t]

    T_fixed=-4*N*(p_A*log(p_A)+(1-p_A)*log(1-p_A))

    first_momenT_matrix_fixed=zeros(t,0)
    second_momenT_matrix_fixed=zeros(t,0)

    first_moment_array=[]
    second_moment_array=[]
    variance_array=[]

    for _ in 1:n_tries
        
        creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

        # Loop over the time ###############################################################################

        population_A=[]

        for _ in 1:t

            
            new_population=next_generation(population_matrix, N)

            population_matrix=new_population

            counter_A=count(x -> x == 1, new_population)/(2*N)

            push!(population_A, counter_A)
            
        end

        first_momenT_matrix_fixed=hcat(first_momenT_matrix_fixed, population_A)
        second_momenT_matrix_fixed=hcat(second_momenT_matrix_fixed, population_A.^2)
    
    end

    for i in 1:length(time_array)
        first=mean(first_momenT_matrix_fixed[i,:])
        second=mean(second_momenT_matrix_fixed[i,:])
        push!(first_moment_array, first)
        push!(second_moment_array, second)
        push!(variance_array, second-first^2 )
    end


    plot(time_array, first_moment_array, xlabel="generations", label=L"\left\langle x\left(t\right)\right\rangle", grid=false, title="x=$p_A")
    plot!(time_array, second_moment_array, label=L"\left\langle x^2\left(t\right)\right\rangle")
    plot!(time_array, variance_array, label=L"\sigma_{x}^{2}(t)")
    vline!([T_fixed], label=L"\bar{t}", color=:black)
    ylims!(-0.1, 1.1)
    savefig("2-moments_fig_$N-$t-$n_tries-$p_A.png")
    
end

#moments_trajectories(500, 0.7, 5000, 10^4)
#moments_trajectories(500, 0.3, 5000, 10^4)
#moments_trajectories(500, 0.5, 5000, 10^4)



function t_fixed_and_t_lost(N_array, #population 
    p_A, #frequency for the allele A
    t,   #time of evolution of the population
    n_tries
)

    p_A_teor=[i for i in 0.01:0.01:1]

    T_matrix_fixed=zeros(length(p_A), 0)
    T_matrix_lost=zeros(length(p_A), 0)
    fixation=true
    
    for i in 1:2
    
        if fixation==true
            plot(grid=false, xlabel="x", ylabel=L"\bar{t_1}(x)", title="Time until fixation")
        else
            plot(grid=false, xlabel="x", ylabel=L"\bar{t_0}(x)", title="Time until allele lost")
        end

        for N in N_array

            population_matrix=zeros(N,2) #array with one of the alleles of every individual
            t_fixed_teor=-4*N .*(1 .-p_A_teor) .*log.(1 .-p_A_teor)./ p_A_teor
            t_lost_teor=-4*N .*(p_A_teor).*log.(p_A_teor)./(1 .-p_A_teor)
            t_array=[]
            t_teor_fixed=[]
            t_teor_lost=[]
            t_lost_array=[]
            t_rel_fixed=[]
            t_rel_lost=[]
            conf_fix_array=[]
            conf_lost_array=[]

            for p in p_A

                T_fixed=-4*N*(1-p)*log(1-p)/p
                T_lost=-4*N*p*log(p)/(1-p)

                push!(t_teor_fixed, T_fixed)
                push!(t_teor_lost, T_lost)


                t_fixed_mean_array=[]
                t_lost_mean_array=[]

                for _ in 1:n_tries
                    
                    fixed=false
                    counter_A=0.0

                    creation_of_the_alleles_of_the_population(population_matrix, p, N)

                    # Loop over the time ###############################################################################
                    

                    for i in 1:t

                        
                        new_population=next_generation(population_matrix, N)

                        population_matrix=new_population
                        
                        counter_A=count(x -> x == 1, new_population)/(2*N)

                        if counter_A== 1.0 && fixed==false

                            push!(t_fixed_mean_array, i)
                            fixed=true
                        
                        elseif counter_A== 0.0 && fixed==false

                            push!(t_lost_mean_array, i)
                            fixed=true

                        end

                    end

                end


                mean_fixed=mean(t_fixed_mean_array)
                mean_lost=mean(t_lost_mean_array)
                rel_err_fix=abs(mean_fixed-T_fixed)*100/T_fixed
                rel_err_lost=abs(mean_lost-T_lost)*100/T_lost
                confidence_fix=calculate_ci(t_fixed_mean_array, 0.95)
                confidence_lost=calculate_ci(t_lost_mean_array, 0.95)

                push!(t_array, mean_fixed)
                push!(t_lost_array, mean_lost)
                push!(t_rel_fixed, rel_err_fix)
                push!(t_rel_lost, rel_err_lost)
                push!(conf_fix_array, confidence_fix)
                push!(conf_lost_array, confidence_lost)

                

            end
            
            if fixation==true
                scatter!(p_A, t_array, label="N=$N")
                plot!(p_A_teor, t_fixed_teor, label="Theoretical result for N=$N")
            else 
                scatter!(p_A, t_lost_array, label="N=$N")
                plot!(p_A_teor, t_lost_teor, label="Theoretical result for N=$N")
            end

            T_matrix_fixed=hcat(T_matrix_fixed, t_teor_fixed)
            T_matrix_fixed=hcat(T_matrix_fixed, t_array)
            T_matrix_fixed=hcat(T_matrix_fixed, t_rel_fixed)
            T_matrix_fixed=hcat(T_matrix_fixed, conf_fix_array)


            T_matrix_lost=hcat(T_matrix_lost, t_teor_lost)
            T_matrix_lost=hcat(T_matrix_lost, t_lost_array)
            T_matrix_lost=hcat(T_matrix_lost, t_rel_lost)
            T_matrix_lost=hcat(T_matrix_lost, conf_lost_array)


        end

        if fixation==true
            savefig("time_until_fixation.png")
            writedlm("T_fixed_-p_A.txt", T_matrix_fixed)
        else
            savefig("time_until_loss.png")
            writedlm("T_lost_-p_A.txt", T_matrix_lost)
        end
        fixation=false
    end
    
end

p_A_array=[i for i in 0.1:0.1:0.9]
N_array=[50, 100]

#t_fixed_and_t_lost(N_array, p_A_array, 2000, 1000)


function heterozygosis_evolution(N, #population 
    p_A, #frequency for the allele A
    t,   #time of evolution of the population
    n_tries)
    

    time_array=[i for i in 1:1:t]

    heterozygosity_matrix=zeros(length(time_array), 0)

    for i in 1:n_tries

        population_matrix=zeros(N,2) #array with one of the alleles of every individual
        creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

        # Loop over the time ###############################################################################

        heterozygosity_array_time=[]

        #counter_A=0
        

        for _ in 1:t

            heterozygosity=0 

            new_population=next_generation(population_matrix, N)

            population_matrix=new_population

            #counter_A=count(x -> x == 1, new_population)/(2*N)

            for n in 1:N

                if population_matrix[n,1]!=population_matrix[n,2]
                    heterozygosity +=1
                end
            end

            push!(heterozygosity_array_time, heterozygosity/N)

        end

        heterozygosity_matrix=hcat(heterozygosity_matrix, heterozygosity_array_time)
    end

    mean_heterozygosis_array=[]
    std_heterozygosis_array=[]

    for i in 1:length(time_array)
        push!(mean_heterozygosis_array, mean(heterozygosity_matrix[i,:]))
        push!(std_heterozygosis_array, std(heterozygosity_matrix[i,:]))
    end

    model = lm(@formula(y ~ x), DataFrame(x=time_array , y=log.(mean_heterozygosis_array)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    y_teor=slope.*time_array.+intercept

    scatter(time_array, log.(mean_heterozygosis_array), linewidth=2, 
    label="", xlabel="Generations", ylabel=L"\ln(H_t)", grid=false)
    plot!(time_array, y_teor,linewidth=2, label="Linear fit")
    
    savefig("Evolution_of_heterozygosis_$p_A.png")

    println("")
    println("slope= ", slope," ± ", slope_error)
    println("")
    println("r_squared= ", r_squared)
    println("")
    print("Sigmas:")
    println((-slope-1/(2*N))/slope_error)

end

heterozygosis_evolution(100, 0.5,100,100) 
heterozygosis_evolution(100, 0.1,100,100)

function variance(N, #population 
    p_A, #frequency for the allele A
    t,   #time of evolution of the population
    n_tries)

    population_matrix=zeros(N,2) #array with one of the alleles of every individual
    counter_matrix=zeros(t, 0)


    for _ in 1:n_tries
        

        creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

        # Loop over the time ###############################################################################

        time_array=[]
        population_A=[]

        for i in 1:t

            push!(time_array, i)
            
            new_population=next_generation(population_matrix, N)

            population_matrix=new_population

            counter_A=count(x -> x == 1, new_population)/(2*N)

            push!(population_A, counter_A)
            
        end

        counter_matrix=hcat(counter_matrix, population_A)
    
    end

    counter_matrix_formated = [round(i, digits=3, RoundNearestTiesAway) for i in counter_matrix]

    list_counter=[]

    for i in 1:size(counter_matrix_formated,1)
        for j in 1:size(counter_matrix_formated,2)
            if counter_matrix_formated[i,j]∉ list_counter
                push!(list_counter, counter_matrix_formated[i,j])
            end
        end
    end

    list_counter_float = Float64.(list_counter)

    counter_dict=Dict{Float64, Float64}()

    for i in list_counter_float
        var_array=[]
        for j in 1:size(counter_matrix_formated,1)
            for k in 1:size(counter_matrix_formated,2)
                if counter_matrix_formated[j,k]==i && j!=size(counter_matrix_formated,1)
                    push!(var_array, counter_matrix_formated[j+1,k])
                end
            end
        end
        variance=var(var_array)
        if isnan(variance)
            counter_dict[i]=0
        else
            counter_dict[i]=variance
        end
    end

    frequencies_array = collect(keys(counter_dict))
    variance_array = collect(values(counter_dict))


    model(x, p) = x .* (1 .- x) / (2 .*p[1])

    initial_guess = [1.0] 

    fit = curve_fit(model, frequencies_array, variance_array, initial_guess)

    # Get the fitted parameters
    params = coef(fit)
    param_errors = fit.stderror


    freq_sorted=sort(frequencies_array, rev=true)
    var_teor= freq_sorted .* (1 .-freq_sorted)/(2 *params[1])

    println("Fitted N: ", params[1])
    println("Relative error: ", abs(params[1]-N)*100/N)
    println("standard deviation", param_errors)
    scatter(frequencies_array, variance_array, xlabel="x", ylabel=L"Var(x_{t+1}|x_t=x)", grid=false, label="")
    plot!(freq_sorted, var_teor, linewidth=3, label=L"Var_{teor}")
    savefig("Var_frequencies_$N.png")

end

#variance(1000, 0.5, 5000,1000)
#variance(100, 0.5, 500,1000)

population_array=[20, 40, 60, 80, 100, 200, 500]

function homozygosis_evolution(N_array, #population 
    p_A, #frequency for the allele A
    t,   #time of evolution of the population
    n_tries)
    

    time_array=[i for i in 1:1:t]

    plot(xlabel="t (Generations)", grid=false)
    for N in N_array

        homozygosity_matrix=zeros(length(time_array), 0)

        y_teor=1 .-exp.(-1/(2*N).*time_array)

        for i in 1:n_tries

            population_matrix=zeros(N,2) #array with one of the alleles of every individual
            creation_of_the_alleles_of_the_population(population_matrix,p_A, N)

            # Loop over the time ###############################################################################

            homozygosity_array_time=[]

            #counter_A=0
            

            for _ in 1:t

                homozygosity=0 

                new_population=next_generation(population_matrix, N)

                population_matrix=new_population

                #counter_A=count(x -> x == 1, new_population)/(2*N)

                for n in 1:N

                    if population_matrix[n,1]==population_matrix[n,2]
                        homozygosity +=1
                    end
                end

                push!(homozygosity_array_time, homozygosity/N)

            end

            homozygosity_matrix=hcat(homozygosity_matrix, homozygosity_array_time)
        end

        mean_homozygosis_array=[]

        for i in 1:length(time_array)
            push!(mean_homozygosis_array, mean(homozygosity_matrix[i,:]))
        end
        
        scatter!(time_array, mean_homozygosis_array, ms=3, label="N=$N")
        plot!(time_array, y_teor, label="Theoretical result for N=$N", linewidth = 2)

    end

    savefig("Inbreeding_figure.png")

end

#homozygosis_evolution(population_array, 0.5, 300, 100)
