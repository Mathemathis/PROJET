using JuMP
using CPLEX

include("utils/parsing.jl")
include("utils/utils_heuristic.jl")

function plneCompacte(n::Int64, s::Int64, t::Int64, S::Int64, p::Vector{Int64}, d::Dict{Any, Any}, name_instance, is_perturbated, timelimit, deltap, deltam)
    """Crée un modèle et résout le problème compact (sans incertitudes)"""

    if is_perturbated == "Yes"
        for key in keys(d)
            d[key]+=rand(1:1000)/100000
        end
    end

    m = Model(CPLEX.Optimizer)
    #set_silent(m)

    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
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

    if feasibleSolutionFound
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance, deltap, deltam)
    """Donnne une solution initiale qui respecte la contrainte robuste pour les sommets
    retourne les valeurs de x et a"""

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
    initial_values=plneCompacte(n, s, t, S, p, d, name_instance, "No", 60, deltap, deltam)
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
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function transformSol(a, n::Int64, s::Int64, t::Int64, ph::Vector{Int64}, d::Dict{Any, Any}, p, deltap, deltam)
    """Prend une solution réalisable  avec la valeur de a
    renvoie 
    - le chemin 
    - les indices des sommets empruntés par ordre décroissant de ph i_ph_dec[3] = nom du 3 eme sommet par ordre decroissant des ph
    - un dictionnaire i_to_i_ph_dec[sommet 14] = ordre de 14 dans i_ph_dec"""
    chemin = [s] # construction du chemin et des aretes
    aretes = []
    current_node=s
    a[s] = 0
    while current_node!=t
        for j in collect(deltap[current_node])
            if a[j]>= 1 - 1e-5
                push!(chemin, j)
                push!(aretes, (current_node, j))
                current_node = j
                a[j] = 0
                break
            end
        end
    end
    i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] >= ph[y])
    i_to_i_ph_dec=Dict() # passer du numero du sommet a sa position dans i_ph_dec
    for i in collect(chemin)
        i_to_i_ph_dec[i]= findfirst(x -> x == i, i_ph_dec)
    end
    return(chemin, i_ph_dec, i_to_i_ph_dec)
end

function nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, nv_noeud, old_noeud, d2, ph, p, longueur)
    """test si le nouveau voisin est admissible pour les poids des sommets"""
    i_old_noeud = i_to_i_ph_dec[old_noeud]
    nv_poids = sum_poids + p[nv_noeud] - p[old_noeud]
    if i_old_noeud <= i_lim # nv_noeud compte dans le poids
        nv_poids -= (2*ph[old_noeud]) # on enleve ce poids
        if ((i_lim == longueur) || (ph[nv_noeud] >= ph[i_ph_dec[(i_lim + 1)]])) # attention peut être placé juste à droite
            nv_poids += (2*ph[nv_noeud]) # ajout du noeud
            # a droite de i_lim reste pareil
        else # on ajoute pas le nv noeud et decalage à droite
            nv_poids +=  2* ph[i_ph_dec[(i_lim + 1)]]
            nv_poids -= ph[i_ph_dec[(i_lim + 1)]]* (d2 - 2 * i_lim)  # on ajoute la variation du sommet à droite
            if i_lim + 2 <= longueur
                if  ph[nv_noeud] >= ph[i_ph_dec[(i_lim + 2)]]
                    nv_poids += (ph[nv_noeud]* (d2 - 2 * i_lim))
                else
                    nv_poids += (ph[i_ph_dec[(i_lim + 2)]]* (d2 - 2 * i_lim))
                end
            end
        end 
    else  # ancien noeud pas enleve
        if ((i_lim > 0) && (ph[nv_noeud] > ph[i_ph_dec[i_lim]])) # sinon on ne change rien
            
            nv_poids += (2*ph[nv_noeud]) # on sature le noeud
            nv_poids -= (2*ph[i_ph_dec[i_lim]]) # on enleve le noeud limite
            nv_poids +=(ph[i_ph_dec[i_lim]] * (d2 - 2 * i_lim)) # saturation partielle du noeud limite
            nv_poids -= (ph[i_ph_dec[(i_lim+1)]]* (d2 - 2 * i_lim)) # le sommet à droite de la limite n'est plus saturé
        else 
            if ((i_old_noeud == i_lim+1) || (ph[nv_noeud] >  ph[i_ph_dec[(i_lim+1)]])) # il est placé juste à droite
                nv_poids += (ph[nv_noeud] * (d2 - 2 * i_lim))
                nv_poids -= (ph[i_ph_dec[(i_lim+1)]]* (d2 - 2 * i_lim))  # saturation partielle du noeud limite
            end
        end
    end
    return(nv_poids)
end

function VoisAmeliorant(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs)
    i = 2
    while (i <= longueur-1) 
        sommets_admissibles = intersect(deltap[chemin[i-1]], deltam[chemin[i+1]]) # il existe un chemin
        sommets_admissibles = filter(x -> x != chemin[i], sommets_admissibles) # enlever le chemin actuel
        println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
        for nv_noeud in collect(sommets_admissibles)
            nv_poids = nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, nv_noeud, chemin[i], d2, ph, p, longueur) # se fait en temps constant (youpi !)
            if nv_poids <= S
                if nvDist(chemin, nv_noeud, chemin[i], d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isAdmissible(chemin, nv_noeud, chemin[i], d2, p, ph, deltap, S))
                    println("nv distances ? ", nvDist(chemin, nv_noeud, chemin[i], d1, d, D))
                    return(chemin[i], nv_noeud)
                end
            end
        end
        i += 1
    end
    return(-1, -1)
end

function deplacement(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs, longueur, i_lim, d1, d2, p, ph, d, D, old_noeud, nv_noeud)
<<<<<<< HEAD
    """on enleve old noeud du chemin et on met nv_noeud, retourne les informations"""
    # i_lim, longueur restent égaux
    i_old_noeud = i_to_i_ph_dec[old_noeud]
    sum_poids = nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, nv_noeud, old_noeud, d2, ph, p, longueur) # utiliser le temps constant

    chemin = nvChemin(chemin, old_noeud, nv_noeud)

    i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] >= ph[y])
    i_to_i_ph_dec=Dict() # passer du numero du sommet a sa position dans i_ph_dec
    for i in collect(chemin)
        i_to_i_ph_dec[i]= findfirst(x -> x == i, i_ph_dec)
    end
    sum_arcs =  Dist(chemin, d1, d, D)
    return(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs)
end

function voisinages(name_instance)
    # preparation solution
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance, deltap, deltam)
    chemin, i_ph_dec, i_to_i_ph_dec  = transformSol(a, n, s, t, ph, d, p, deltap, deltam)
    longueur = length(chemin)
    sum_poids, i_lim = getInfoSommets(chemin, p, ph, d2)
    sum_arcs = Dist(chemin, d1, d, D)
    println("Solution admissible")
    println("sum_dist = ", sum_arcs)
    println("sum_poids = ", sum_poids, ", S = ", S,  "\n")

    # debut des voisinages
    condition = true
    while condition 
        old_noeud, nv_noeud = VoisAmeliorant(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs)
        println("old_noeud = ", old_noeud, " nv_noeud = ", nv_noeud, "\n")
        if old_noeud == -1
            condition = false
        else
            chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs = deplacement(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs, longueur, i_lim, d1, d2, p, ph, d, D, old_noeud, nv_noeud)
        end
    end
=======
    
>>>>>>> 295799a7bf43f537708437639b84bcc1b784bcd6
end

function voisAdmissibles(chemin, deltap)
    println("debut recherche voisinage admissibles")
    for i in 2:(length(chemin)-1)
        if chemin[i+1] in collect(deltap[chemin[i-1]])
            println(chemin[i-1], chemin[i], chemin[i+1])
        end
    end
end
function main()
<<<<<<< HEAD
    #name_instance="100_USA-road-d.BAY.gr"
    name_instance="900_USA-road-d.NY.gr"
    #voisinages(name_instance)

    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance, deltap, deltam)
    chemin, i_ph_dec, i_to_i_ph_dec  = transformSol(a, n, s, t, ph, d, p, deltap, deltam)
    voisAdmissibles(chemin, deltap)
=======
    name_instance="100_USA-road-d.BAY.gr"
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")

    # init delta p delta m
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

    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance, deltap, deltam)
    println("solution est construite")
    #=
    ph[14] = 7
    p[15] = 18
    ph[15] = 6
    d2 = 5 =#

    chemin, i_ph_dec, i_to_i_ph_dec  = transformSol(a, n, s, t, ph, d, p, deltap, deltam)
    println("chemin = ", chemin)
    println("i_ph_dec = ", i_ph_dec)
    println("i_to_i_ph_dec = ", i_to_i_ph_dec)

    longueur = length(chemin)

    sum_poids, i_lim = getInfoSommets(chemin, p, ph, d2)
    #=
    nv_noeud = 15
    old_noeud = 14

    nv_poids = nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, nv_noeud, old_noeud, d2, ph, p, longueur)
    println("nv_poids = ", nv_poids)

    nv_chemin = nvChemin(chemin, old_noeud, nv_noeud)
    nv_poids_test, _ = getInfoSommets(nv_chemin, p, ph, d2)
    println("nv_poids_test = ", nv_poids_test)
    =#


   
    sum_arcs = Dist(chemin, d1, d, D)
    println("sum_arcs  = ", sum_arcs)
    println("chemin = ", chemin)

   old_noeud, nv_noeud = VoisAmeliorant(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, i_lim, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs)

>>>>>>> 295799a7bf43f537708437639b84bcc1b784bcd6
end