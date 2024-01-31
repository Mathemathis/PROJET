include("./utils/parsing.jl")
include("./utils/utils_heuristic.jl")
include("constrSol.jl")

function RechLocSup(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-1)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+1]]) # il existe un chemin
        println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
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
        i += 1
    end
    return(false, -1)
end


function nvCheminSup(chemin, noeud, nv_noeud)
    """Saute mouton / au lieu de chemin[i] -> chemin[i+1] -> chemin[i+2] -> chemin[i+3]
    on fait chemin[i] -> nv_noeud -> chemin[i+3]"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin)
        push!(nv_chemin, i)
        if chemin[i] == noeud
            push!(nv_chemin, nv_noeud)
        end
        i += 1
    end
    return(nv_chemin)
end