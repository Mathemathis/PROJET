using JuMP
using CPLEX
include("PLNE_compacte.jl")
include("parsing.jl")

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
    renvoie le chemin, les indices des sommets empruntés par ordre décroissant de ph et des arêtes empruntées par ordre décroissant de d"""
    chemin = [s] # construction du chemin et des aretes
    aretes = []
    current_node=s
    a[s] = 0
    while current_node!=t
        println("chemin = ", chemin)
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
    i_p_dec =sort(chemin, lt = (x, y) -> ph[x] >= ph[y])
    i_aretes_d = sort(aretes, lt = (x, y) -> d[x] >= d[y])
    ph_dec = [ph[i] for i in i_p_dec]
    p_dec = [p[i] for i in i_p_dec]

    #aretes_d = [d[i] for i in i_aretes_d]
    return(chemin, i_p_dec, p_dec, ph_dec, i_aretes_d)
end

function getInfoSommets(i_p_dec, p, ph, d2)
    """Renvoie le poids des sommets dans le cas robuste et tous les sommets de 1 à i_res sont chargés au maximum
    i_p_dec : liste des indices des sommets triées par poids decroissant de ph"""
    res=sum([p[i] for i in collect(i_p_dec)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    i_res = 0 # dernier indice à être rempli totalement (indice dans i_poids_croissants)
     
    for i in i_p_dec
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
    return res, i_res # total des poids max, indice du dernier sommet rempli totalement
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
    return res, i_res
end

function nvResPoids(i_p_dec, p_dec, ph_dec, sum_poids, i_lim, nv_noeud, i_old_noeud, d2, ph, p, longueur)
    """test si le nouveau voisin est admissible pour les poids des sommets"""
    nv_poids = sum_poids + p[nv_noeud] - p_dec[i_old_noeud]
    if i_old_noeud <= i_lim # nv_noeud compte dans le poids
        nv_poids -= (2*ph_dec[i_old_noeud]) # on enleve ce poids
        if ((i_lim == longueur) || (ph[nv_noeud] >= ph_dec[i_lim + 1])) # attention peut être placé juste à droite
            nv_poids += (2*ph[nv_noeud]) # ajout du noeud
            # a droite de i_lim reste pareil
        else # on ajoute pas le nv noeud et decalage à droite
            nv_poids +=  + 2*ph_dec[(i_lim + 1)] 
            nv_poids -= ph_dec[(i_lim +1)]* (d2 - 2 * i_lim)  # on ajoute la variation du sommet à droite
            if i_lim + 2 <= longueur
                if  ph[nv_noeud] >= ph_dec[(i_lim + 2)]
                    nv_poids += (ph[nv_noeud]* (d2 - 2 * i_lim))
                else
                    nv_poids += (ph_dec[i_lim+ 2]* (d2 - 2 * i_lim))
                end
            end
        end 
    else  # ancien noeud pas enleve
        if ((i_lim > 0) && (ph[nv_noeud] > ph_dec[i_lim])) # sinon on ne change rien
            nv_poids += (2*ph[nv_noeud]) # on sature le noeud
            nv_poids -= (2*ph_dec[i_lim]) # on enleve le noeud limite
            nv_poids +=(ph_dec[i_lim] * (d2 - 2 * i_lim)) # saturation partielle du noeud limite
            nv_poids -= (ph_dec[(i_lim +1)]* (d2 - 2 * i_lim)) # le sommet à droite de la limite n'est plus saturé
        else 
            if ((i_old_noeud == i_lim+1) || (ph[nv_noeud] >  ph_dec[(i_lim+1)])) # il est placé juste à droite
                nv_poids += (ph[nv_noeud] * (d2 - 2 * i_lim))
                nv_poids -= (ph_dec[(i_lim +1)]* (d2 - 2 * i_lim))  # saturation partielle du noeud limite
            end
        end
    end
    return(nv_poids)
end

function checkChemin(chemin, nv_noeud, old_noeud, d2, p, ph)
    """Maniere longue de calculer le poids d'un chemin voisin pour verifier les calculs"""
    nv_chemin = copy(chemin)
    for i in 1:length(chemin) # calcul du nouveau chemin
        if nv_chemin[i] == old_noeud
            nv_chemin[i] = nv_noeud
        end
    end
    i_p_dec =sort(nv_chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    return(getInfoSommets(i_p_dec, p, ph, d2))
end

function isAdmissible(chemin, nv_noeud, old_noeud, d2, p, ph, deltap, S)
    nv_chemin = copy(chemin)
    for i in 1:length(chemin) # calcul du nouveau chemin
        if nv_chemin[i] == old_noeud
            nv_chemin[i] = nv_noeud
        end
    end
    for i in 1:(length(nv_chemin)-1)
        if !(in(nv_chemin[i+1], deltap[nv_chemin[i]]))
            return(false) # on ne peut pas tracer le chemin
        end
    end 
    i_p_dec =sort(nv_chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    res_poids, _ = getInfoSommets(i_p_dec, p, ph, d2)
    return(res_poids <= S) # poids 
end

function nvDist(chemin, nv_noeud, old_noeud, d1, d, D)
    """Calcul de la nouvelle distance des aretes d'un chemin"""
    aretes = []
    nv_chemin = copy(chemin)
    current_node = nv_chemin[1]
    for i in 2:length(nv_chemin) # calcul du nouveau chemin
        if nv_chemin[i] == old_noeud
            nv_chemin[i] = nv_noeud
        end
        push!(aretes, (current_node, nv_chemin[i]))
        current_node = nv_chemin[i]
    end
    i_aretes_d =sort(aretes, lt = (x, y) -> d[x] <= d[y], rev = true)
    sum_arcs, _  = getInfoArcs(i_aretes_d, d, D, d1)
    return(sum_arcs)
end

function VoisAmeliorant(chemin, i_p_dec, p_dec, ph_dec, sum_poids, i_lim, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, current_sol)
    found = false
    i = 2
    while (i <= longueur-1) && (!found)
        droite = chemin[i-1]
        gauche = chemin[i+1]
        i_old_noeud = findfirst(x -> x == chemin[i], i_p_dec)
        println("old_noeud = ", chemin[i])
        sommets_admissibles = intersect(deltap[droite], deltam[gauche])
        sommets_admissibles = filter(x -> x != chemin[i], sommets_admissibles)
        println("sommets admissibles = ", sommets_admissibles)
        for nv_noeud in collect(sommets_admissibles)
            nv_poids = nvResPoids(i_p_dec, p_dec, ph_dec, sum_poids, i_lim, nv_noeud, i_old_noeud, d2, ph, p, longueur)
            if nv_poids <= S
                println("on a trouve un sommet qui respecte la contrainte des poids")
                if nvDist(chemin, nv_noeud, chemin[i], d1, d, D) < current_sol
                    println("on a trouve une solution ameliorante")
                    found = true
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

function main()
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
    chemin, i_p_dec, p_dec, ph_dec, i_aretes_d  = transformSol(a, n, s, t, ph, d, p, deltap, deltam)

    longueur = length(chemin)
    sum_poids, i_lim = getInfoSommets(i_p_dec, p, ph, d2)
    current_sol, arc_limite = getInfoArcs(i_aretes_d, d, D, d1)
    println("current_sol  = ", current_sol)
    println("chemin = ", chemin)

    old_noeud, nv_noeud = VoisAmeliorant(chemin, i_p_dec, p_dec, ph_dec, sum_poids, i_lim, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, current_sol) 
end




