include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")
include("vois1.jl")
include("vois2.jl")

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









function voisinages(name_instance)
    # preparation solution
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    chemin, i_ph_dec, i_to_i_ph_dec  = initSolDijkstra(n, s, t, S, d1, d2, p, ph, d, D)
    deltap, deltam = initDelta(d, n)

    longueur = length(chemin)
    sum_poids = getInfoSommets(chemin, p, ph, d2)
    sum_arcs = Dist(chemin, d1, d, D)

    # debut des voisinages
    condition = true
    while condition 
        trouve, old_noeud, nv_noeud = RechLocEchange(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs, s, t)
        if !trouve
            condition = false
        else
            chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs = deplacementEch(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs, longueur, d1, d2, p, ph, d, D, old_noeud, nv_noeud)
        end
    end
    return(chemin)
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
    name_instance="500_USA-road-d.BAY.gr"
    nv_chemin = voisinages(name_instance)

    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_poids = getInfoSommets(nv_chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println(" Est bien un chemin ?", isChemin(nv_chemin, deltap, s, t))
    println("Poids du chemin = ", nv_poids)
    println("valeur de S = ", S)
    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))
end