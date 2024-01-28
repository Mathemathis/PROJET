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

function transformSol(a, n::Int64, s::Int64, t::Int64, ph::Vector{Int64}, d::Dict{Any, Any}, p)
    """Prend une solution réalisable 
    renvoie le chemin, les indices des sommets empruntés par ordre décroissant de ph et des arêtes empruntées par ordre décroissant de d"""
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
    i_poids_d =sort(chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    i_aretes_d = sort(aretes, lt = (x, y) -> d[x] <= d[y], rev = true)
    poids_h_d = [ph[i] for i in i_poids_d]
    poids_d = [p[i] for i in i_poids_d]

    aretes_d = [d[i] for i in i_aretes_d]
    return(chemin, i_poids_d, i_aretes_d, poids_d, poids_h_d, aretes_d)
end

function getInfoSommets(i_poids_d, p, ph, d2)
    """Renvoie le poids des sommets dans le cas robuste et tous les sommets de 1 à i_res sont chargés au maximum
    i_poids_d : liste des indices des sommets triées par poids decroissant de ph"""
    res=sum([p[i] for i in collect(i_poids_d)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    i_res = 0 # dernier indice à être rempli totalement (indice dans i_poids_croissants)
     
    for i in i_poids_d
        if capa+2<=d2
            capa+=2 # augmentation du budget
            res+=2*ph[i] # augmentation du poids total des sommets
            i_res += 1 # on rempli totatlement ce sommet
        else 
            res+=(d2-capa)*ph[i]
            capa = 0
            break
        end
    end
    return res, i_res, capa # total des poids max, indice du dernier sommet rempli totalement
end

function getInfoArcs(i_aretes_d, d, D, d1)
    """Renvoie la distance du chemin dans le cas robuste et tous les arcs de 1 à i_res sont chargés au maximum
    i_aretes_d : liste des indices aretes triées par poids decroissant de d"""

    res=sum([d[i] for i in collect(i_aretes_d)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    i_res = 0  # dernier indice à être rempli totalement (indice dans i_aretes_croissantes)
    
    for (i,j) in i_aretes_d
        if capa+D[i,j]<=d1
            capa+=D[i,j]
            res+=d[i,j]*D[i,j]
            i_res -= 1
        else
            res+=d[i,j]*(d1-capa)
            capa = d1
            break
        end
    end
    return res, i_res, capa
end

function voisAdmissible(i_poids_d, poids_d, poids_h_d, sum_poids, i_lim, nv_noeud, i_old_noeud, d2, ph, p, longueur)
    """test si le nouveau voisin est admissible pour les poids des sommets"""
    
    nv_poids = sum_poids + p[nv_noeud] - poids_d[i_old_noeud]

    
    if i_old_noeud <= i_lim # nv_noeud compte dans le poids
        nv_poids -= (2*poids_h_d[i_old_noeud]) # on enleve ce poids
        if ((i_lim == longueur) || (ph[nv_noeud] >= poids_h_d[i_lim + 1])) # attention peut être placé juste à droite
            nv_poids += (2*ph[nv_noeud]) # ajout du noeud
            # a droite de i_lim reste pareil
        else # on ajoute pas le nv noeud et decalage à droite
            nv_poids +=  + 2*poids_h_d[(i_lim + 1)] 
            nv_poids -= poids_h_d[(i_lim +1)]* (d2 - 2 * i_lim)  # on ajoute la variation du sommet à droite
            if i_lim + 2 <= longueur
                if  ph[nv_noeud] >= poids_h_d[(i_lim + 2)]
                    nv_poids += (ph[nv_noeud]* (d2 - 2 * i_lim))
                else
                    nv_poids += (poids_h_d[i_lim+ 2]* (d2 - 2 * i_lim))
                end
            end
        end 
    else  # ancien noeud pas enleve
        if ((i_lim > 0) && (ph[nv_noeud] > poids_h_d[i_lim])) # sinon on ne change rien
            nv_poids += (2*ph[nv_noeud]) # on sature le noeud
            nv_poids -= (2*poids_h_d[i_lim]) # on enleve le noeud limite
            nv_poids +=(poids_h_d[i_lim] * (d2 - 2 * i_lim)) # saturation partielle du noeud limite
            nv_poids -= (poids_h_d[i_lim +1]* (d2 - 2 * i_lim)) # le sommet à droite de la limite n'est plus saturé
        else 
            if ph[nv_noeud] >  poids_h_d[i_lim+1] # il est placé juste à droite
                nv_poids += ph[nv_noeud] * (d2 - 2 * i_lim)
                nv_poids -= poids_h_d[(i_lim +1)]* (d2 - 2 * i_lim)  # saturation partielle du noeud limite
            end
        end
    end
    return(nv_poids)
end

function checkChemin(chemin, nv_noeud, old_noeud, d2, p, ph)
    current_node = chemin[1]
    for i in 1:length(chemin) # calcul du nouveau chemin
        if chemin[i] == old_noeud
            chemin[i] = nv_noeud
        end
    end
    i_poids_d =sort(chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    return(getInfoSommets(i_poids_d, p, ph, d2))
end

function nvDist(chemin, nv_noeud, old_noeud, d1, d, D)
    aretes = []
    current_node = chemin[1]
    for i in 2:length(chemin) # calcul du nouveau chemin
        if chemin[i] == old_noeud
            chemin[i] = nv_noeud
        end
        push!(aretes, (current_node, chemin[i]))
        current_node = chemin[i]
    end
    println("aretes = ", aretes)
    i_aretes_d =sort(aretes, lt = (x, y) -> d[x] <= d[y], rev = true)
    return(getInfoArcs(i_aretes_d, d, D, d1))
end


function main()

    name_instance="20_USA-road-d.NY.gr"

    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
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
    print("d2 = ", d2)
    d2 = 1
    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance)
    chemin, i_poids_d, i_aretes_d, poids_d, poids_h_d, aretes_d = transformSol(a, n, s, t, ph, d, p)
    longueur = length(chemin)
    println("sommets poids_croissants : ", i_poids_d)
    println("poids decroissants : ", poids_d)

    sum_poids, i_lim, capa_sommets = getInfoSommets(i_poids_d, p, ph, d2)
    sum_arcs, arc_limite, capa_arcs = getInfoArcs(i_aretes_d, d, D, d1)
    println("resultat sommets :", sum_poids)
    println("indice sommets :", i_lim)
    println("capa = ", capa_sommets)
    i_old_noeud = 2
    old_noeud = chemin[i_old_noeud]
    nv_noeud = 8
    println("old noeud : p[$old_noeud] = ", p[old_noeud])
    println("old noeud : ph[$old_noeud] = ", ph[old_noeud])
    p[nv_noeud] = 24
    ph[nv_noeud] = 3
    D[2,8] = 1
    println("nv noeud : p[$nv_noeud] = ", p[nv_noeud])
    println("nv noeud : ph[$nv_noeud] = ", ph[nv_noeud])

    nv_poids = voisAdmissible(i_poids_d, poids_d, poids_h_d, sum_poids, i_lim, nv_noeud, i_old_noeud, d2, ph, p, longueur)
    poids_ref, _, _ = checkChemin(chemin, nv_noeud, old_noeud, d2, p, ph)
    nv_dist, _, _ = nvDist(chemin, nv_noeud, old_noeud, d1, d, D)
    println("methode rapide = ", nv_poids)
    println("methode lente = ", poids_ref)   
    println("ancienne distance = ", sum_arcs) 
    println("nv_dist = ", nv_dist)

    println("D[2, 14] = ", D[2,14])
    println("d[2, 14] = ", d[2,14])
    println("D[2, 8] = ", D[2,8])
    println("D[2, 8] = ", d[2,8])
    println("D[14, 17] = ", D[14,17])
    println("d[14, 17] = ", d[14,17])
    println("D[8, 17] = ", D[8,17])
    println("d[8, 17] = ", d[8,17])

end




