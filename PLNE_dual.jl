using JuMP
using CPLEX
include("PLNE_compacte.jl")

function plne_dual(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, name_instance, is_perturbated, with_objective_problem, initial_value, timelimit)
    """Résout le PLNE par la méthode duale (question 4)"""

    # initialisation des voisinages entrants et sortants pour chaque sommet
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

    #Dualisation objectif

    if with_objective_problem == "Yes"

        @variable(m, alpha>=0)
        @variable(m, beta[i in 1:n, j in deltap[i]]>=0)

        @constraint(m,  [i in 1:n, j in deltap[i]], alpha+beta[i,j]>=x[i,j]*d[i,j])
    end

    #Dualisation Contrainte poids

    @variable(m, gamma>=0)
    @variable(m, eta[i in 1:n]>=0)

    @constraint(m,  [i in 1:n], gamma+eta[i]>=a[i]*ph[i])
    @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma<=S)

    if with_objective_problem == "Yes"
        @objective(m, Min, sum(x[i,j]*d[i,j] for i in 1:n for j in deltap[i])+alpha*d1+sum(beta[i,j]*D[i,j] for i in 1:n for j in deltap[i]))
    else
        @objective(m, Min, sum(x[i,j]*d[i,j] for i in 1:n for j in deltap[i]) )
    end

    if initial_value=="compacte"
        initial_values=plne_compacte(n, s, t, S, p, d, name_instance, "No", 60)
        JuMP.set_start_value.(x, initial_values[1])
        JuMP.set_start_value.(a, initial_values[2])
    end
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    if is_perturbated == "Yes"
        if with_objective_problem == "Yes"
            if initial_value == "Yes"
                logfile_name = "txtFiles/plne_dual/Fonction_obj_robuste/perturbation/with_initial_value/$name_instance.txt"
            else
                logfile_name = "txtFiles/plne_dual/Fonction_obj_robuste/perturbation/without_initial_value/$name_instance.txt"
            end
        else
            if initial_value == "Yes"
                logfile_name = "txtFiles/plne_dual/Fonction_obj_non_robuste/perturbation/with_initial_value/$name_instance.txt"
            else
                logfile_name = "txtFiles/plne_dual/Fonction_obj_non_robuste/perturbation/without_initial_value/$name_instance.txt"
            end
        end
    else
        if with_objective_problem == "Yes"
            if initial_value == "Yes"
                logfile_name = "txtFiles/plne_dual/Fonction_obj_robuste/no_perturbation/with_initial_value/$name_instance.txt"
            else
                logfile_name = "txtFiles/plne_dual/Fonction_obj_robuste/no_perturbation/without_initial_value/$name_instance.txt"
            end
        else
            if initial_value == "Yes"
                logfile_name = "txtFiles/plne_dual/Fonction_obj_non_robuste/no_perturbation/with_initial_value/$name_instance.txt"
            else
                logfile_name = "txtFiles/plne_dual/Fonction_obj_non_robuste/no_perturbation/without_initial_value/$name_instance.txt"
            end
        end
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        close(logfile)
        return JuMP.value.(x), JuMP.value.(a), JuMP.value.(eta), JuMP.value.(gamma)
    end
    close(logfile)
end

function check_dual_poids(n::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, a::Vector{Float64}, timelimit)

    m = Model(CPLEX.Optimizer)

    @variable(m, delta_2[1:n]>=0)    

    @constraint(m,  [i in 1:n], delta_2[i]<=2)

    @constraint(m,  sum(delta_2[i] for i in 1:n)<=d2)
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    logfile_name = "txtFiles/check_dual_poids.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    @objective(m, Max, sum(a[i]*(p[i]+delta_2[i]*ph[i]) for i in 1:n))

    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        return ("Valeur de l’objectif : ", JuMP.objective_value(m))
    end
end

function check_dual_objective(n::Int64, d1::Int64, d::Dict{Any, Any}, D::Dict{Any, Any}, x, timelimit)

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

    m = Model(CPLEX.Optimizer)

    @variable(m, delta_1[i in 1:n, j in deltap[i]]>=0)    

    @constraint(m,  [i in 1:n, j in deltap[i]], delta_1[i,j]<=D[i,j])

    @constraint(m,  sum(delta_1[i,j] for i in 1:n for j in deltap[i])<=d1)
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    logfile_name = "txtFiles/check_dual_objective.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    @objective(m, Max, sum(x[i,j]*d[i,j]*(1+delta_1[i,j]) for i in 1:n for j in deltap[i]))

    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        return ("Valeur de l’objectif : ", JuMP.objective_value(m))
    end
end