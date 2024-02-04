include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")
include("utils/voisinages.jl")



function initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    bool, chemin = dijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    if bool == false
        @warn("On a pas trouve de chemin de poids <= S")
    else
        return(chemin)
    end
end

function Vk(k, chemin, d2, ph, p, d1, d, D, deltap, deltam, S, dist, s, t)
    if k == 1
        return RechLocEch(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, dist, s, t)
    else 
        if k == 2
            return RechLocInf(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, dist, s, t)
        else
            if k == 3
                return RechLocSup(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, dist, s, t)
            else
                return(false, [])
            end
        end
    end
end

function rechLoc(name_instance)
    """méthode de la descente de voisinages 
    retourne la borne superieur trouve (distance des aretes) et le chemin associé pour vérifications"""
    timelimit = time() + 60 # limite de temps fixée à 60 secondes

    # preparation solution
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    chemin  = initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    deltap, deltam = initDelta(d, n)

    dist = Dist(chemin, d1, d, D) # meilleur distance atteinte
    iter = 0

    # debut des voisinages
    k = 1
    while k <= 3
        if time() > timelimit
            #println("on sort à cause du temps")
            return(dist, chemin)
        end
        trouve, nv_chemin =  Vk(k, chemin, d2, ph, p, d1, d, D, deltap, deltam, S, dist, s, t)
        if !trouve
            k += 1
            #println("ON CHERCHE DANS LE VOISINAGE ", k)
        else
            chemin = nv_chemin
            dist = Dist(chemin, d1, d, D)
            k = 1
            iter += 1
        end
    end
    return(dist, chemin, iter)
end

function main()
    name_instance="1500_USA-road-d.COL.gr"
    
    dist, nv_chemin, _ = @time rechLoc(name_instance)

    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println("Est bien un chemin ? ", isChemin(nv_chemin, deltap, s, t))
    println("Poids du chemin = ", nv_poids)
    println("valeur de S = ", S)
    println("Distance du chemin = ", dist)
end