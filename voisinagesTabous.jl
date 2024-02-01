include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")


function getBestVoisTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous)
    """Retourne le meilleur mouvement qui n'est pas dans les tabous"""
    res = RechLocEchTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous)
    res = vcat(res, RechLocInfTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous))
    res = vcat(res, RechLocSupTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous))
    if res != []
        min_chemin = argmin(x -> x[2], res)
        min_index = findfirst(x -> x == min_chemin, res)
        min_chemin = res[min_index][1]
        return(true, min_chemin)
    else
        return(false, [])
    end
end

function IsAdmissibleTabous(chemin, tabous)
    aretes = cheminToAretes(chemin)
    for t in tabous  # on teste toutes les aretes
        if t in aretes
            return false
        end
    end
    return true
end


function RechLocEchTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous)
    res = []
    for (j, noeud) in enumerate(chemin[2:(end-1)])
        i = j + 1
        sommets_admissibles = intersect(deltap[chemin[i-1]], deltam[chemin[i+1]]) # il existe un chemin
        sommets_admissibles = filter(x -> x != noeud, sommets_admissibles) # enlever le chemin actuel
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin = nvCheminEch(chemin, noeud, nv_noeud)
            if IsAdmissibleTabous(nv_chemin, tabous)
                if !isChemin(chemin, deltap, s, t)
                    @warn("ceci n'est pas un chemin")
                end
                nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
                if nv_poids <= S
                    nv_dist = Dist(nv_chemin, d1, d, D) 
                    push!(res, [nv_chemin, nv_dist])
                end
            end
        end
    end
    return(res)
end


function RechLocInfTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous)
    res = []
    for (i, noeud) in enumerate(chemin[1:(end-3)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+3]]) # il existe un chemin
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin =  nvCheminInf(chemin, chemin[i], nv_noeud)
            if IsAdmissibleTabous(nv_chemin, tabous)
                if !isChemin(chemin, deltap, s, t)
                    @warn("ceci n'est pas un chemin")
                end
                nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
                if nv_poids <= S
                    nv_dist = Dist(nv_chemin, d1, d, D) 
                    push!(res, [nv_chemin, nv_dist])
                end
            end
        end
    end
    return(res)
end


function RechLocSupTabous(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, s, t, tabous)
    res = []
    for (i, noeud) in enumerate(chemin[1:(end-1)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+1]]) # il existe un chemin
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin =  nvCheminSup(chemin, chemin[i], nv_noeud)
            if IsAdmissibleTabous(nv_chemin, tabous)
                if !isChemin(chemin, deltap, s, t)
                    @warn("ceci n'est pas un chemin")
                end
                nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
                if nv_poids <= S
                    nv_dist = Dist(nv_chemin, d1, d, D) 
                    push!(res, [nv_chemin, nv_dist])
                end
            end
        end
    end
    return(res)
end