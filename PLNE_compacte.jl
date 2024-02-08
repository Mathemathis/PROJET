using JuMP
using CPLEX
include("./utils/parsing.jl")
include("./utils/utils_heuristic.jl")
include("./constrSol.jl")

function plne_compacte(name_instance, is_perturbated, timelimit)
    """Crée un modèle et résout le problème compact (sans incertitudes)"""
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    deltap, deltam=initDelta(d, n)

    Keys=collect(keys(D))
    d_set, Kd_set, A_d_set=calcul_d_k_set(p, d, S, 100000)
    # p_set, Kp_set, A_p_set=calcul_p_k_set(p, S)
    
    initial_values=get_init_sol(name_instance)

    if is_perturbated == "Yes"
        for key in keys(d)
            d[key]+=rand(1:1000)/100000
        end
    end

    m = Model(CPLEX.Optimizer)

    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    @constraint(m,  sum(a[i]*p[i] for i in 1:n)<=S)

    if is_perturbated=="No"

        @variable(m, y[d_i in d_set, k in 1:Kd_set[d_i]], Bin)
        @constraint(m,  [d_i in d_set], sum(k*y[d_i,k] for k in 1:Kd_set[d_i]) ==sum(x[i,j] for (i,j) in A_d_set[d_i]))
        @constraint(m,  [d_i in d_set], sum(y[d_i,k] for k in 1:Kd_set[d_i]) <=1)
    end
    
    # @variable(m, wp[p_i in p_set, k in 1:Kp_set[p_i]], Bin)
    # @constraint(m,  [p_i in p_set], sum(k*wp[p_i,k] for k in 1:Kp_set[p_i]) ==sum(a[i] for i in A_p_set[p_i]))
    # @constraint(m,  [p_i in p_set], sum(wp[p_i,k] for k in 1:Kp_set[p_i]) <=1)

    # @constraint(m, sum(p_i*k*wp[p_i,k] for p_i in p_set for k in 1:Kp_set[p_i]) <= S)
    x_last=[]

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
    
            CPLEX.load_callback_variable_primal(cb_data, context_id) # on regardera les valeurs des variables dans la solution courante

            x_val=Dict(key => callback_value(cb_data, x[key]) for key in Keys) # on stocke x
            a_val = [callback_value(cb_data, a[i]) for i in 1:n]
            nodes_visited= [i for i in 1:n if a_val[i]>1e-5]
            edges_visited= [key for key in Keys if x_val[key]>1e-5]
            if edges_visited!=x_last
                println(nodes_visited, sum(x_val[key]*d[key] for key in Keys))
                y_val=Dict((d_i,k) => callback_value(cb_data, y[d_i, k]) for d_i in d_set for k in 1:Kd_set[d_i])
                K_set_star=Dict(d_i => k for d_i in d_set for k in 1:Kd_set[d_i] if y_val[d_i,k]>1e-5)
                cstr3 = @build_constraint(sum(sum(y[d_i,k] for k in K_set_star[d_i]:Kd_set[d_i]) for d_i in keys(K_set_star))<= sum(y_val[d_i, k] for d_i in d_set for k in 1:Kd_set[d_i] if y_val[d_i,k]>1e-5)-1+sum(x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr3)
                vars= vcat([x[i,j] for i in 1:n for j in deltap[i]], [a[i] for i in 1:n], [y[d_i, k] for d_i in d_set for k in 1:Kd_set[d_i]])
    
                vals= vcat([callback_value(cb_data,x[i,j]) for i in 1:n for j in deltap[i]], [callback_value(cb_data,a[i]) for i in 1:n], [callback_value(cb_data,y[d_i, k]) for d_i in d_set for k in 1:Kd_set[d_i]])
                
                MOI.submit(m, MOI.HeuristicSolution(cb_data), vars, vals)
                x_last=edges_visited
            end
        end
    end

    JuMP.set_start_value.(x, JuMP.Containers.SparseAxisArray(initial_values[1]))
    JuMP.set_start_value.(a, initial_values[2])


    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_HeuristicEffort", 0)
    
    if is_perturbated=="No"
        MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    end

    if is_perturbated == "Yes"
        logfile_name = "txtFiles/plne_compacte/perturbation/$name_instance.txt"
    else
        logfile_name = "txtFiles/plne_compacte/no_perturbation/$name_instance.txt"
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    #logfile_path = abspath(logfile_name)
    #logfile = open(logfile_path, "w")
    #redirect_stdout(logfile)

    @objective(m, Min, sum(x[i,j]*d[i,j] for i in 1:n for j in deltap[i])) 
    # @objective(m, Min, sum(d_i*k*y[d_i,k] for d_i in d_set for k in 1:Kd_set[d_i]))

    # Résolution d’un modèle
    start = time()
    optimize!(m)
    computation_time = time() - start

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println(JuMP.objective_value(m))
        return name_instance, computation_time, JuMP.objective_value(m), objective_bound(m)
    end
    #close(logfile)
end

# function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

#     # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
#     # appelé dans un des deux cas suivants :
#     # cas 1 - une solution entière a été obtenue; ou
#     # cas 2 - une relaxation non bornée a été obtenue
#     if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
#         return false
#     end

#     # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
#     # solution entière courante
#     ispoint_p = Ref{Cint}()
#     ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

#     # S'il n'y a pas de solution entière
#     if ret != 0 || ispoint_p[] == 0
#         return false
#     else
#         return true
#     end
# end