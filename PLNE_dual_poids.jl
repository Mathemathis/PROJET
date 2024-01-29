using JuMP
using CPLEX
include("PLNE_compacte.jl")

function plne_dual_poids(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d, D)
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

    m = Model(CPLEX.Optimizer)
    set_silent(m)

    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    # @variable(m, a[1:n]>=0)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)


    #Dualisation Contrainte poids

    @variable(m, gamma>=0)
    @variable(m, eta[i in 1:n]>=0)

    @constraint(m,  [i in 1:n], gamma+eta[i]>=a[i]*ph[i])
    @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma<=S)

    @objective(m, Min, sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma-S )

    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)
    
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_Solutions", 1)


    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        return JuMP.value.(x), JuMP.value.(a), JuMP.value.(eta), JuMP.value.(gamma)
    end
end