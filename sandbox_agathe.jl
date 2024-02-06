include("plans_coupants.jl")
include("heurVois.jl")
include("./utils/parsing.jl")
include("branch_cut.jl")


function main()
    name_instance="1000_USA-road-d.COL.gr"
    #n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    
    #results_ALG=plans_coupantsALG(n, s, t, S, d1, d2, p, ph, d, D, name_instance, "no_symmetry", "with initial values", 300)    


    #println("\n \n HEURISTIQUE \n")
    #nv_chemin = voisinages(name_instance)

    #=v_poids = @time getInfoSommets(nv_chemin, p, ph, d2)
    deltap, deltam = initDelta(d, n)
    println("Est bien un chemin ? ", isChemin(nv_chemin, deltap, s, t))
    println("Poids du chemin = ", nv_poids)
    println("valeur de S = ", S)
    println("Distance du chemin = ", Dist(nv_chemin, d1, d, D))=#

    branchAndCut(name_instance, "no_symmetry", "with initial values", 30)
end