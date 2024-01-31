include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")

"""Garde un nombre de noeud constants
chemin[i-1] -> chemin[i] = old_noeud -> chemin[i+1]
chemin[i-1] -> nv_noeud -> chemin[i+1]"""

function RechLocEch(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)

    for (j, noeud) in enumerate(chemin[2:(end-1)])
        i = j + 1
        sommets_admissibles = intersect(deltap[chemin[i-1]], deltam[chemin[i+1]]) # il existe un chemin
        sommets_admissibles = filter(x -> x != noeud, sommets_admissibles) # enlever le chemin actuel
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin = nvCheminEch(chemin, noeud, nv_noeud)
            nv_poids_lent = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids_lent <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs - 1e-5 # temps lineaire en le nombre d'aretes
                    println("____________________________")
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isChemin(chemin, deltap, s, t))
                    println("nv distances ? ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")
                    return(true, nv_chemin)
                end
            end
        end
        i += 1
    end
    return(false, -1)
end


function nvCheminEch(chemin, old_noeud, nv_noeud)
    """Nouveau chemin ou on remplace old_noeud par un nouveau noeud"""
    nv_chemin = copy(chemin) # calcul du nouveau chemin
    for i in 1:length(chemin) 
        if nv_chemin[i] == old_noeud
            nv_chemin[i] = nv_noeud
        end
    end
    return(nv_chemin)
end