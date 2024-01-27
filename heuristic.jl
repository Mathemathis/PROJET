using JuMP
using CPLEX
include("PLNE_compacte.jl")
include("parsing.jl")

function plne_compacte_no_print(n::Int64, s::Int64, t::Int64, S::Int64, p::Vector{Int64}, d::Dict{Any, Any}, name_instance, is_perturbated, timelimit)
    """Crée un modèle et résout le problème compact (sans incertitudes)"""
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
    set_silent(m)

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
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance)
    """Crée un modèle et trouve un chemin minimisant le poids robuste sur les sommets
    retourne les valeurs de x et a"""
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
    @variable(m, a[1:n], Bin)
    @variable(m, gamma>=0)
    @variable(m, eta[i in 1:n]>=0)

    # chemin réalisable
    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    # dual problème des sommets
    @constraint(m,  [i in 1:n], gamma+eta[i]>=a[i]*ph[i])
    @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma<=S)

    @objective(m, Min, 0)

    # solution de départ
    initial_values=plne_compacte_no_print(n, s, t, S, p, d, name_instance, "No", 60)
    JuMP.set_start_value.(x, initial_values[1])
    JuMP.set_start_value.(a, initial_values[2])
    
    # set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    # Résolution du modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Value a : ", JuMP.value.(a))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function transformSol(a, n::Int64, s::Int64, t::Int64, ph::Vector{Int64}, d::Dict{Any, Any})
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
    chemin = [s]
    aretes = []
    current_node=s
    while current_node!=t
        println("current_node : ", current_node, " ")
        for j in collect(deltap[current_node])
            if a[j]>= 1 - 1e-5
                push!(chemin, j)
                push!(aretes, (current_node, j))
                current_node = j
                break
            end
        end
    end
    println("chemin : ", chemin)
    ind_poids_croissants =sort(chemin, lt = (x, y) -> ph[x] <= ph[y])
    ind_aretes_croissantes = sort(aretes, lt = (x, y) -> d[x] <= d[y])
    return(chemin, ind_poids_croissants, ind_aretes_croissantes)
end

function getInfoSommets(ind_poids_croissants, p, ph, d2)
    res=sum([p[i] for i in collect(ind_poids_croissants)]) # somme deterministe
    println("res init ", res)

    capa=0 # budget pour augmenter les delta^2
    l = length(ind_poids_croissants)
    ind_res = l +  1# dernier indice à être rempli totalement (indice dans ind_poids_croissants)
     
    for i in reverse(ind_poids_croissants)
        if capa+2<=d2
            capa+=2 # augmentation du budget
            res+=2*ph[i] # augmentation du poids total des sommets
            ind_res -= 1 # on rempli totatlement ce sommet
            println("res =", res)
        else 
            res+=(d2-capa)*ph[i]
            break
        end
    end
    println("res final = ", res)
    return res, ind_res # total des poids max, indice du dernier sommet rempli totalement
end

function getInfoArcs(ind_aretes_croissantes, d, D, d1)
    res=sum([d[i] for i in collect(ind_aretes_croissantes)]) # somme deterministe

    capa=0 # budget pour augmenter les delta^2
    l = length(ind_aretes_croissantes)
    ind_res = l + 1  # dernier indice à être rempli totalement (indice dans ind_aretes_croissantes)
    
    for (i,j) in reverse(ind_aretes_croissantes)
        if capa+D[i,j]<=d1
            capa+=D[i,j]
            res+=d[i,j]*D[i,j]
            ind_res -= 1
        else
            res+=d[i,j]*(d1-capa)
            break
        end
    end
    return res, ind_res
end


function main()
    name_instance="20_USA-road-d.NY.gr"
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    print("d2 = ", d2)
    d2 = 1
    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance)
    chemin, ind_poids_croissants, ind_aretes_croissantes = transformSol(a, n, s, t, ph, d)
    println("sommets poids_croissants : ", ind_poids_croissants)
    println("p poids_croissants : ", [p[x] for x in collect(ind_poids_croissants)])
    println("ph poids_croissants : ", [ph[x] for x in collect(ind_poids_croissants)])
    println("aretes croissantes : ", ind_aretes_croissantes)
    println("aretes poids_croissants : ", [d[x] for x in collect(ind_aretes_croissantes)])
    res_sommets, ind_sommets = getInfoSommets(ind_poids_croissants, p, ph, d2)
    res_arcs, ind_arcs = getInfoArcs(ind_aretes_croissantes, d, D, d1)
    println("resultat sommets :", res_sommets)
    println("indice sommets :", ind_sommets)
end



main()
