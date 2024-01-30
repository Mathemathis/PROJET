include("./utils/parsing.jl")
include("./utils/utils_heuristic.jl")
include("constrSol.jl")
include("heuristic.jl")


function RechLocSup(chemin, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs, s, t)
    i = 1
    while (i <= longueur-3) 
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+3]]) # il existe un chemin
        println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin =  nvCheminSup(chemin, chemin[i], nv_noeud)
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")
                    return(true, chemin[i], nv_noeud)
                end
            end
        end
        i += 1
    end
    return(false, -1, -1)
end