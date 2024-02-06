include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
using DataStructures

function dijkstraLambda(s, t, d1, d2, p, ph, d, D, deltap, lambda)
    timelimit = 60
    nds_visites = []
    obj = PriorityQueue()
    obj[s] = (1 - lambda) * (p[s] + min(2, d2) * ph[s])
    chemin_emprunte = Dict()
    chemin_emprunte[s] = [s]
    while !isempty(obj)
        nd_courant = dequeue!(obj)
        if nd_courant == t
            return(chemin_emprunte[t])
        end
        push!(nds_visites, nd_courant)
        voisins_non_visites = setdiff(deltap[nd_courant], nds_visites)
        for i in collect(voisins_non_visites)
            chemin = vcat(chemin_emprunte[nd_courant], [i]) # mise a jour du chemin
            nv_poids  = getInfoSommets(chemin, p, ph, d2)

            if lambda > 0
                nv_dist = Dist(chemin, d1, d, D)
                nv_obj = (1 - lambda) * nv_poids + lambda * nv_dist
            else 
                nv_obj = nv_poids
            end

            if haskey(obj, i)
                if nv_obj <= obj[i] - 1e-5
                    obj[i] = nv_obj
                    chemin_emprunte[i] = chemin 
                end
            else
                obj[i] = nv_obj
                chemin_emprunte[i] = chemin
            end
        end
    end 
end

function RechDich(name_instance)
    """ méthode de la recherche dichotomique en faisant varier lambda
    retourne la borne superieure trouvée (distance sur les aretes) et le chemin qui y est associé (pour verifications)"""
    start_time = time() + timelimit # limite de temps fixée à 60 secondes
    # preparation solution
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)
   

    chemin_sup = dijkstraLambda( s, t, d1, d2, p, ph, d, D, deltap, 0) # prise en compte uniquement des poids
    borne_sup = Dist(chemin_sup, d1, d, D)

    lambda = 1/2
    iter = 0
    a = 0
    b = 1

    iter = 0 # nb de fois ou on a une borne sup améliorante

    while b - a > 1e-5 && iter < 100
        if time() > start_time
            #println("on sort à cause du temps")
            return(borne_sup, chemin_sup, iter)
        end
        lambda = (a + b)/2
        chemin = dijkstraLambda( s, t, d1, d2, p, ph, d, D, deltap, lambda)
        poids = getInfoSommets(chemin, p, ph, d2)
        if poids <= S + 1e-7 
            a = lambda
            dist =  Dist(chemin, d1, d, D)
            if dist < borne_sup
                borne_sup = dist
                chemin_sup = chemin
                iter += 1
                #println("borne sup = ", borne_sup)
            end
        else
            b = lambda
        end        
    end
    return(borne_sup, chemin_sup, iter)

end


function main()
    name_instance="200_USA-road-d.COL.gr"
    bsup, chemin, _ = @time RechDich(name_instance)


    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam = initDelta(d, n)

    poids = getInfoSommets(chemin, p, ph, d2)
    valeur = Dist(chemin, d1, d, D)

    println("chemin  = ", chemin)
    println("est ce bien un chemin ?", isChemin(chemin, deltap, s, t))
    println("S = ", S)
    println("poids de notre chemin = ", poids, "\n")
    println("valeur obj = ", valeur)
    println(" [", ", ", bsup, "]")
end
