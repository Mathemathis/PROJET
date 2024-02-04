using DataFrames
using CSV
include("heurVois.jl")


function main()
    csv_file_path = "results/heurVois/resultsHeurVois.csv"  
    files = readdir("data/")[2:end] # noms des fichiers

    result_df = DataFrame(instance = [], time=[], ubound=[], nb_iter=[])
    for name_instance in files
        println("instance : ", name_instance)
        start = time()
        dist, _, iter = rechLoc(name_instance)
        time_elapsed = time() - start


        push!(result_df, (instance = name_instance, time = time_elapsed, ubound=dist, nb_iter = iter))
        CSV.write(csv_file_path, DataFrame(result_df))
    end    
end