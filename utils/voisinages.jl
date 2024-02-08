include("parsing.jl")
include("utils_heuristic.jl")

"""Garde un nombre de noeud constants
chemin[i-1] -> chemin[i] = old_noeud -> chemin[i+1]
chemin[i-1] -> nv_noeud -> chemin[i+1]"""

function RechLocEch(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-2)])
        sommets_admissibles = intersect(deltap[chemin[i]], deltam[chemin[i+2]]) # il existe un chemin
        sommets_admissibles = filter(x -> x != chemin[i+1], sommets_admissibles) # enlever le chemin actuel
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_noeud in collect(sommets_admissibles)
            nv_chemin = nvCheminEch(chemin, chemin[i+1], nv_noeud)
            if !isChemin(nv_chemin, deltap, s, t)
                println("nv_chemin echange = ", nv_chemin)
                println("old nd =", chemin[i+1], " nv nd  = ", nv_noeud)
                println("deltam[$(chemin[i+2])] = ", deltam[chemin[i+2]])
                @warn("RechLocEch - ceci n'est pas un chemin")
            end
            nv_poids_lent = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids_lent <= S
                #println("Solution admissible")
                
                if Dist(nv_chemin, d1, d, D) < sum_arcs - 1e-5 # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isChemin(chemin, deltap, s, t))
                    println("nv distances ? ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")=#
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
            if !isChemin(nv_chemin, deltap, s, t)
                @warn("RechLocInf - ceci n'est pas un chemin")
            end
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")=#
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
            if !isChemin(nv_chemin, deltap, s, t)
                @warn("RechLocSup - ceci n'est pas un chemin")
            end
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________") =#
                    return(true, nv_chemin)
                end
            end
        end
    end
    return(false, -1)
end

function RechLocInf2(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-4)])
        sommets_admissibles= []
        sommets_k = deltap[chemin[i]]
        sommets_j = deltam[chemin[i+4]]
        for k in sommets_k
            for l in sommets_j
                if l in deltap[k]
                    push!(sommets_admissibles,[k,l])
                end
            end
        end
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_nds in collect(sommets_admissibles)
            k = nv_nds[1]
            l = nv_nds[2]
            nv_chemin =  nvCheminInf2(chemin, chemin[i], k, l)
            if !isChemin(nv_chemin, deltap, s, t)
                @warn("RechLocInf - ceci n'est pas un chemin")
            end
            nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("Chemin dmissible ? ",isChemin(nv_chemin, deltap, s, t))
                    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")=#
                    return(true, nv_chemin)
                end
            end
        end
    end
    return(false, -1)
end


function RechLocEch2(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-3)])
        sommets_admissibles= []
        sommets_k = deltap[chemin[i]]
        sommets_j = deltam[chemin[i+3]]
        for k in sommets_k
            for l in sommets_j
                if l in deltap[k]
                    push!(sommets_admissibles,[k,l])
                end
            end
        end
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_nds in collect(sommets_admissibles)
            k = nv_nds[1]
            l = nv_nds[2]
            nv_chemin = nvCheminEch2(chemin, noeud, k, l)
            if !isChemin(nv_chemin, deltap, s, t)
                println("chemin = ", chemin)
                println("noeud = ", noeud)
                println("k, l =", k, ", ", l)
                @warn("RechLocEch2 - ceci n'est pas un chemin")
            end
            nv_poids_lent = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids_lent <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs - 1e-5 # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isChemin(chemin, deltap, s, t))
                    println("nv distances ? ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")=#
                    return(true, nv_chemin) 
                end
            end
        end
    end
    return(false, -1)
end


function RechLocSup2(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    for (i, noeud) in enumerate(chemin[1:(end-2)])
        sommets_admissibles= []
        sommets_k = deltap[chemin[i]]
        sommets_j = deltam[chemin[i+2]]
        for k in sommets_k
            for l in sommets_j
                if l in deltap[k]
                    push!(sommets_admissibles,[k,l])
                end
            end
        end
        #println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", noeud)
        for nv_nds in collect(sommets_admissibles)
            k = nv_nds[1]
            l = nv_nds[2]
            nv_chemin = nvCheminSup2(chemin, noeud, k, l)
            if !isChemin(nv_chemin, deltap, s, t)
                @warn("RechLocEch - ceci n'est pas un chemin")
            end
            nv_poids_lent = getInfoSommets(nv_chemin, p, ph, d2)
            if nv_poids_lent <= S
                #println("Solution admissible")
                if Dist(nv_chemin, d1, d, D) < sum_arcs - 1e-5 # temps lineaire en le nombre d'aretes
                    #=println("____________________________")
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isChemin(chemin, deltap, s, t))
                    println("nv distances ? ", Dist(nv_chemin, d1, d, D))
                    println("____________________________")=#
                    return(true, nv_chemin) 
                end
            end
        end
    end
    return(false, -1)
end