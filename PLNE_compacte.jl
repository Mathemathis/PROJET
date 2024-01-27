using JuMP
using CPLEX

function plne_compacte(n::Int64, s::Int64, t::Int64, S::Int64, p::Vector{Int64}, d::Dict{Any, Any}, name_instance, is_perturbated, timelimit)
    """Crée une modèle et résoute le problème compact (sans incertitudes)"""
    deltap=Dict()
    deltam=Dict()
    for i in 1:n
        deltap[i]=[]
        deltam[i]=[]
    end
    for (i,j) in keys(d)
        push!(deltap[i],j)
        push!(deltam[j],i)
    end

    if is_perturbated == "Yes"
        for key in keys(d)
            d[key]+=rand(1:1000)/100000
        end
    end

    m = Model(CPLEX.Optimizer)

    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    # @variable(m, a[1:n]>=0)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    @constraint(m,  sum(a[i]*p[i] for i in 1:n)<=S)

    
    # set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    if is_perturbated == "Yes"
        logfile_name = "txtFiles/plne_compacte/perturbation/$name_instance.txt"
    else
        logfile_name = "txtFiles/plne_compacte/no_perturbation/$name_instance.txt"
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    @objective(m, Min, sum(x[i,j]*d[i,j] for i in 1:n for j in deltap[i]))

    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Value a : ", JuMP.value.(a))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        close(logfile)
        return JuMP.value.(x), JuMP.value.(a)
    end
    close(logfile)
end