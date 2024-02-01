include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
using DataStructures

function dijkstra(name_instance)
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    print("s = ", s, ", t = ", t)
    deltap, deltam = initDelta(d, n)

    nds_visites = []
    distance = PriorityQueue()
    distance[s] = 0
    chemin_emprunte = Dict()
    chemin_emprunte[s] = [s]
    for _ in 1:n
        nd_courant = dequeue!(distance)
        push!(nds_visites, nd_courant)

        voisins_non_visites = setdiff(deltap[nd_courant], nds_visites)
        for i in collect(voisins_non_visites)
            chemin = vcat(chemin_emprunte[nd_courant], [i]) # mise a jour du chemin
            nv_dist, _ = getInfoSommets(chemin, p, ph, d2)
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
                    println("chemin pour t :", chemin_emprunte[t])
                    return(true, chemin_emprunte[t])
                else
                    println("on a croise t mais son poids est : ", distance[t] )
                end
            end
        end
    end
    return(false, [])
end

function main()
    name_instance="2500_USA-road-d.BAY.gr"
    bool, chemin = dijkstra(name_instance)


    # test
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_dist, _ = getInfoSommets(chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println(" Est bien un chemin ?", isChemin(chemin, deltap))
    println("distance du chemin = ", nv_dist)
    println("valeur de S = ", S)
end

function get_init_sol(name_instance)
    bool, chemin = dijkstra(name_instance)


    # test
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    nv_dist, _ = getInfoSommets(chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println(" Est bien un chemin ?", isChemin(chemin, deltap))
    println("distance du chemin = ", nv_dist)
    println("valeur de S = ", S)
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