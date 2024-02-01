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
        nd_courant = dequeue!(obj)
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



function RechDich( n, s, t, S, d1, d2, p, ph, d, D, deltap)
    chemin_inf = dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, 1) # prise en compte uniquement des aretes
    chemin_sup = dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, 0) # prise en compte uniquement des poids

    borne_inf =  Dist(chemin_inf, d1, d, D)
    borne_sup = Dist(chemin_sup, d1, d, D)

    lambda = 1/2
    iter = 0
    a = 0
    b = 1
    while b - a > 1e-5 && iter < 100
        #println("b inf = ", borne_inf, ", b sup = ", borne_sup)
        #println("lambda = ", lambda)
        #println("[ ", a, ", ", b, "]")
        lambda = (a + b)/2
        chemin = dijkstraPoids( n, s, t, S, d1, d2, p, ph, d, D, deltap, lambda)
        dist =  Dist(chemin, d1, d, D)
        poids = getInfoSommets(chemin, p, ph, d2)
        if poids <= S + 1e-7 
            #println("solution realisable")
            a = lambda
            if dist < borne_sup
                borne_sup = dist
                chemin_sup = chemin
            end
        else
            b = lambda
        end
        borne_inf = max(dist + (1-lambda)* (poids - S) / lambda, borne_inf)
    end

    return(borne_sup, borne_inf, chemin_sup)

end
function main()
    name_instance="2000_USA-road-d.COL.gr"
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
    bsup, binf, chemin = @time RechDich( n, s, t, S, d1, d2, p, ph, d, D, deltap)

    poids = getInfoSommets(chemin, p, ph, d2)
    valeur = Dist(chemin, d1, d, D)

    println("chemin  = ", chemin)
    println("est ce bien un chemin ?", isChemin(chemin, deltap, s, t))
    println("S = ", S)
    println("poids de notre chemin = ", poids, "\n")
    println("valeur obj = ", valeur)
    println(" [", binf, ", ", bsup, "]")
end
