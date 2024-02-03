include("parsing.jl")
include("utils_heuristic.jl")

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
    end
    return(false, -1)
end


"""Saute mouton / au lieu de chemin[i] -> chemin[i+1] -> chemin[i+2] -> chemin[i+3]
    on fait chemin[i] -> nv_noeud -> chemin[i+3]"""

function RechLocInf(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-3)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+3]]) # il existe un chemin
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin =  nvCheminInf(chemin, chemin[i], nv_noeud)
            if !isChemin(chemin, deltap, s, t)
                @warn("ceci n'est pas un chemin")
            end
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")
                    return(true, nv_chemin)
                end
            end
        end
    end
    return(false, -1)
end



"""au lieu de chemin[i] -> chemin[i+1] 
    on fait chemin[i] -> nv_noeud -> chemin[i+1]"""


function RechLocSup(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-1)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+1]]) # il existe un chemin
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin =  nvCheminSup(chemin, chemin[i], nv_noeud)
            if !isChemin(chemin, deltap, s, t)
                @warn("ceci n'est pas un chemin")
            end
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")
                    return(true, nv_chemin)
                end
            end
        end
    end
    return(false, -1)
end


