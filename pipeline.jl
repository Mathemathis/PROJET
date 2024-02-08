using DataFrames
using CSV
using Distributed
include("heurVois.jl")
include("heurDich.jl")
include("branch_cut.jl")
include("PLNE_compacte.jl")


function pipelineHeurVois()
    csv_file_path = "results/heurVois/resultsHeurVois.csv"   # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers

    result_df = DataFrame(instance = [], time=[], ubound=[], nb_iter=[])
    for name_instance in files
        println("instance : ", name_instance)
        start = time()
        dist, _, iter = rechLoc(name_instance)  # changer ici nom fonction
        time_elapsed = time() - start


        push!(result_df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter))
        CSV.write(csv_file_path, DataFrame(result_df))
    end    
end

function pipelineHeurDich()
    csv_file_path = "results/heurVois/resultsHeurDich.csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers

    #result_df = DataFrame(instance = [], time=[], ubound=[], nb_iter=[])
    result_df = DataFrame(CSV.File(csv_file_path))
    for name_instance in files
        if !(any(x -> occursin(name_instance, string(x)), result_df.instance))
            println("instance : ", name_instance)
            start = time()
            dist, chemin, iter = RechDich(name_instance) # changer ici nom fonction
            time_elapsed = time() - start


            push!(result_df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter))
            CSV.write(csv_file_path, DataFrame(result_df))
        end
    end    
end

function mpBranchCut(csv_name)
    csv_file_path = "results/" * csv_name *".csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers
    result_df = DataFrame(CSV.File(csv_file_path))

    #Threads.@threads for name_instance in files
    for name_instance in files
        if !(any(x -> occursin(name_instance, string(x)), result_df.instance))
            println("instance : ", name_instance)
            println("processor : " , Threads.threadid())

            start = time()
            name_instance, computation_time, chemin, obj, bounds, iter23, iter24 = branchAndCut(name_instance, "no_symmetry", "with initial values", 600)
            vector_string = join(string.(chemin), ", ")
            modified_string = replace(vector_string, "," => "/")

            time_elapsed = time() - start
            nodenames = ["Instance", "Temps", "obj", "bounds", "iter23", "iter24", "chemin"]
            df =  DataFrame([[] for _ = nodenames] , nodenames)
            
            push!(df, (Instance = name_instance, Temps = computation_time, obj=obj, bounds = bounds, iter23 = iter23, iter24 = iter24, chemin = string(chemin)))
            println("on a mis dans result_df")
            CSV.write(csv_file_path, DataFrame(df), append=true)
            
        end
    end
end

function mpHeur(csv_name)
    eta = 1e-2
    csv_file_path = "results/heurDich/" * csv_name * string(eta) * ".csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers
    
    nodenames = ["instance", "time", "ubound", "nb_iter", "chemin"]
    result_df =  DataFrame([[] for _ = nodenames] , nodenames)



    Threads.@threads for name_instance in files
    #Threads.@threads for name_instance in files
        if !(any(x -> occursin(name_instance, string(x)), result_df.instance))
            println("instance : ", name_instance)
            println("processor : " , Threads.threadid())

            start = time()
            dist, chemin, iter = RechDich(name_instance, eta) # changer ici nom fonction
            
            vector_string = join(string.(chemin), ", ")
            modified_string = replace(vector_string, "," => "/")

            time_elapsed = time() - start
            
            df =  DataFrame([[] for _ = nodenames] , nodenames)
            
            push!(df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter, chemin = modified_string))
            println("on a mis dans result_df")
            CSV.write(csv_file_path, DataFrame(df), append=true)
            
        end
    end
end

function mpHeurVois(csv_name)
    k_max = 6
    csv_file_path = "results/heurVois/" * csv_name * string(k_max) * ".csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers
    
    nodenames = ["instance", "time", "ubound", "nb_iter", "chemin"]
    result_df =  DataFrame([[] for _ = nodenames] , nodenames)
    #result_df = DataFrame(CSV.File(csv_file_path))



    Threads.@threads for name_instance in files
    #for name_instance in files
        if !(any(x -> occursin(name_instance, string(x)), result_df.instance))
            println("instance : ", name_instance)
            println("processor : " , Threads.threadid())

            start = time()
            dist, chemin, iter = rechLoc(name_instance, k_max) # changer ici nom fonction
            
            vector_string = join(string.(chemin), ", ")
            modified_string = replace(vector_string, "," => "/")

            time_elapsed = time() - start
            
            df =  DataFrame([[] for _ = nodenames] , nodenames)
            
            push!(df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter, chemin = modified_string))
            println("on a mis dans result_df")
            CSV.write(csv_file_path, DataFrame(df), append=true)
        end
    end
end


function mpCompacte(csv_name)
    pertubated = "Yes"
    csv_file_path = "./results/" * csv_name * pertubated *".csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers
    result_df = DataFrame(CSV.File(csv_file_path))
    #nodenames = ["Instance", "Temps", "UB", "LB" ]
    #result_df =  DataFrame([[] for _ = nodenames] , nodenames)

    #Threads.@threads for name_instance in files
    instance_values = unique(result_df."Instance")
    filtered_elements = setdiff(files, instance_values)
    println("filetered elem = ", filtered_elements)

    Threads.@threads for name_instance in filtered_elements
       
            println("instance : ", name_instance)
            println("processor : " , Threads.threadid())

                
            name_instance, computation_time, UB, LB = plne_compacte(name_instance, pertubated, 600)
                
            nodenames = ["Instance", "Temps", "UB", "LB" ]
            df =  DataFrame([[] for _ = nodenames] , nodenames)
                
            push!(df, (Instance = name_instance, Temps = computation_time, UB=UB, LB = LB))
            println("on a mis dans result_df")
            CSV.write(csv_file_path, DataFrame(df), append=true)
        
    end
end

function main()
    mpCompacte("Compacte")
    #mpBranchCut("BranchCut")

end