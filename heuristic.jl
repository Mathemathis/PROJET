using JuMP
using CPLEX

include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")

function initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    bool, chemin = dijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    if bool == false
        @warn("On a pas trouve de chemin de poids <= S")
    else
        i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] >= ph[y])
        i_to_i_ph_dec=Dict() # passer du numero du sommet a sa position dans i_ph_dec
        for i in collect(chemin)
            i_to_i_ph_dec[i]= findfirst(x -> x == i, i_ph_dec)
        end
        return(chemin, i_ph_dec, i_to_i_ph_dec)
    end

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
    chemin, i_ph_dec, i_to_i_ph_dec  = initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    deltap, deltam = initDelta(d, n)

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
    #name_instance="100_USA-road-d.BAY.gr"
    name_instance="100_USA-road-d.BAY.gr"
    voisinages(name_instance)

    #=n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
    timelimit = 30
    x, a =constrSol(n, s, t, S, p, d, ph, d2, timelimit, name_instance, deltap, deltam)
    chemin, i_ph_dec, i_to_i_ph_dec  = transformSol(a, n, s, t, ph, d, p, deltap, deltam)
    voisAdmissibles(chemin, deltap)=#
end