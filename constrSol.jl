include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
using DataStructures

function dijkstra( n, s, t, S, d1, d2, p, ph, d, D)
    #println("s = ", s, ", t = ", t)
    deltap, deltam = initDelta(d, n)

    nds_visites = []
    distance = PriorityQueue()
    distance[s] = p[s] + min(2, d2) * ph[s]
    chemin_emprunte = Dict()
    chemin_emprunte[s] = [s]
    for _ in 1:n
        nd_courant = dequeue!(distance)
        push!(nds_visites, nd_courant)

        voisins_non_visites = setdiff(deltap[nd_courant], nds_visites)
        for i in collect(voisins_non_visites)
            chemin = vcat(chemin_emprunte[nd_courant], [i]) # mise a jour du chemin
            nv_dist  = getInfoSommets(chemin, p, ph, d2)
            if haskey(distance, i)
                if nv_dist < distance[i]
                    distance[i] = nv_dist
                    chemin_emprunte[i] = chemin 
                end
            else
                distance[i] = nv_dist
                chemin_emprunte[i] = chemin
            end
            if (i == t) 
                if  (distance[t] <= S)
                    return(true, chemin_emprunte[t]) # on sort dès que l'on trouve un chemin pour t admissible
                end
            end
        end
    end
    @warn("La methode dirjska n'a pas trouve de chemin robuste pour les poids")
    return(false, [])
end

function main()
    name_instance="2500_USA-road-d.BAY.gr"
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    bool, chemin = dijkstra( n, s, t, S, d1, d2, p, ph, d, D)


    # test
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_dist = getInfoSommets(chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    # println(" Est bien un chemin ?", isChemin(chemin, deltap, s, t))
    # println("distance du chemin = ", nv_dist)
    # println("valeur de S = ", S)
end

function get_init_sol(name_instance)

    # test
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    bool, chemin = dijkstra( n, s, t, S, d1, d2, p, ph, d, D)
    nv_dist = getInfoSommets(chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    # println(" Est bien un chemin ?", isChemin(chemin, deltap, s, t))
    # println("distance du chemin = ", nv_dist)
    # println("valeur de S = ", S)
    a=[0 for i in 1:n]
    x=Dict(key => 0.0 for key in keys(d))
    for i in 1:(length(chemin)-1)
        x[chemin[i], chemin[i+1]]=1.0
    end
    for node in chemin 
        a[node]=1
    end
    return x, a
end