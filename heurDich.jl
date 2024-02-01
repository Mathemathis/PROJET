include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
using DataStructures

function dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, lambda)
    nds_visites = []
    obj = PriorityQueue()
    obj[s] = (1 - lambda) * (p[s] + min(2, d2) * ph[s])
    chemin_emprunte = Dict()
    chemin_emprunte[s] = [s]
    while !isempty(obj)
        push!(nds_visites, nd_courant)
        voisins_non_visites = setdiff(deltap[nd_courant], nds_visites)

        for i in collect(voisins_non_visites)
            chemin = vcat(chemin_emprunte[nd_courant], [i]) # mise a jour du chemin
            nv_poids  = getInfoSommets(chemin, p, ph, d2)
            nv_dist = Dist(chemin, d1, d, D)
            nv_obj = (1 - lambda) * nv_poids + lambda * nv_dist
            if haskey(obj, i)
                if nv_obj < obj[i] 
                    obj[i] = nv_obj
                    chemin_emprunte[i] = chemin 
                end
            else
                obj[i] = nv_obj
                chemin_emprunte[i] = chemin
            end
        end
    end
    return(chemin_emprunte[t])
end

function dijkstraDist( n, s, t, S, d1, d2, p, ph, d, D, deltap)
    nds_visites = []
    chemin_emprunte = Dict()
    chemin_emprunte[s] = [s]
    distance = PriorityQueue()
    distance[s] = 0
    for _ in 1:n
        nd_courant = dequeue!(distance)
        push!(nds_visites, nd_courant)
        voisins_non_visites = setdiff(deltap[nd_courant], nds_visites)

        for i in collect(voisins_non_visites)
            chemin = vcat(chemin_emprunte[nd_courant], [i]) # mise a jour du chemin
            nv_dist = Dist(chemin, d1, d, D)
            if haskey(distance, i)
                if nv_dist < distance[i]
                    chemin_emprunte[i] = chemin 
                    distance[i] = nv_dist
                end
            else
                chemin_emprunte[i] = chemin
                distance[i] = nv_dist
            end
        end
    end
    return(chemin_emprunte[t])
end

function RechDich( n, s, t, S, d1, d2, p, ph, d, D, deltap)
    chemin_inf = dijkstraDist( n, s, t, S, d1, d2, p, ph, d, D, deltap)
    b_sup, best_chemin = dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, Inf)

    

    a =  Dist(chemin_inf, d1, d, D)
    s1 = getInfoSommets(chemin_inf, p, ph, d2)

    if s1 <= S 1e-5 
        return(chemin_inf, 1e-5)
    end
   
    b =  Dist(best_chemin, d1, d, D)
    m = (a+b)/2
    iter = 0

    println("Premier solution réalisable de valeur ", b)

    # b contient toujours la valeur de la plus petite solution realisable
    # a contient toujours la valeur de la plus grande solution realisable
    while abs(a-b) > 1e-3 && iter < 100
        println(" a = ", a, ", b = ", b)
        m = (a+b)/2
        bool, chemin = dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, m)
        if bool # on a trouve un chemin
            b = m
            best_chemin = chemin
        else
            a = m
        end
        iter += 1
    end
   
    return(best_chemin, abs(a-b))


end
function main()
    name_instance="800_USA-road-d.COL.gr"
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
    chemin, gap = @time RechDich( n, s, t, S, d1, d2, p, ph, d, D, deltap)

    poids = getInfoSommets(chemin, p, ph, d2)
    valeur = Dist(chemin, d1, d, D)

    println("chemin  = ", chemin)
    println("est ce bien un chemin ?", isChemin(chemin, deltap, s, t))
    println("S = ", S)
    println("poids de notre chemin = ", poids, "\n")
    println("valeur = ", valeur)
    println("gap = ", gap)
end
