using DataFrames
using CSV
using Distributed
include("heurVois.jl")
include("heurDich.jl")


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
            dist, _, iter = RechDich(name_instance) # changer ici nom fonction
            time_elapsed = time() - start


            push!(result_df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter))
            CSV.write(csv_file_path, DataFrame(result_df))
        end
    end    
end


function multiprocess(csv_name)
    csv_file_path = "results/" * csv_name *".csv"  # changer ici nom fichier
    files = readdir("data/")[2:end] # noms des fichiers
    result_df = DataFrame(CSV.File(csv_file_path))

    Threads.@threads for name_instance in files
        if !(any(x -> occursin(name_instance, string(x)), result_df.instance))
            println("instance : ", name_instance)
            println("processor : " , Threads.threadid())
            start = time()
            dist, _, iter = rechLoc(name_instance) # changer ici nom fonction
            time_elapsed = time() - start
            nodenames = ["instance", "time", "ubound", "nb_iter"]
            df =  DataFrame([[] for _ = nodenames] , nodenames)
            
            push!(df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter))
            println("on a mis dans result_df")
            CSV.write(csv_file_path, DataFrame(df), append=true)
            
        end
    end
end

function main()
    multiprocess("heurVois")

end