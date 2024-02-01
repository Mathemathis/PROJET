include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")
include("voisinages.jl")



function initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    bool, chemin = dijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    if bool == false
        @warn("On a pas trouve de chemin de poids <= S")
    else
        return(chemin)
    end
end

function Vk(k, chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    if k == 1
        return RechLocEch(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
    else 
        if k == 2
            return RechLocInf(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
        else
            if k == 3
                return RechLocSup(chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
            else
                return(false, [])
            end
        end
    end
end
function voisinages(name_instance)
    # preparation solution
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    chemin  = initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    deltap, deltam = initDelta(d, n)

    sum_arcs = Dist(chemin, d1, d, D)
    println("Distance initiale = ", sum_arcs)

    # debut des voisinages
    k = 1
    while k <= 3
        trouve, nv_chemin =  Vk(k, chemin, d2, ph, p, d1, d, D, deltap, deltam, S, sum_arcs, s, t)
        if !trouve
            k += 1
            println("ON CHERCHE DANS LE VOISINAGE ", k)
        else
            chemin = nv_chemin
            sum_arcs = Dist(chemin, d1, d, D)
            k = 1
        end
    end
    return(chemin)
end


function main()
    name_instance="100_USA-road-d.NY.gr"
    #name_instance="2500_USA-road-d.COL.gr"
    nv_chemin = voisinages(name_instance)

    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println("Est bien un chemin ? ", isChemin(nv_chemin, deltap, s, t))
    println("Poids du chemin = ", nv_poids)
    println("valeur de S = ", S)
    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
end