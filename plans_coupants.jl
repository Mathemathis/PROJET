using JuMP
using CPLEX

function plans_coupants(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
    deltap=Dict()
    deltam=Dict()
    for i in 1:n
        deltap[i]=[]
        deltam[i]=[]
    end
    for (i,j) in keys(d)
        push!(deltap[i],j)
        push!(deltam[j],i)
    end

    m = Model(CPLEX.Optimizer)
    set_silent(m)
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)

    MOI.set(m, MOI.NumberOfThreads(), 1)


    # contraintes du PLNE maitre

    @variable(m, z >=0 )    
    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    #@variable(m, a[1:n]>=0)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    @objective(m, Min, z)

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
    
            CPLEX.load_callback_variable_primal(cb_data, context_id)


            x_val = Matrix{Float64}(undef,n,n) # matrice

            for i in 1:n
                for j in deltap[i]
                    x_val[i,j] = callback_value(cb_data, x[i,j])
                end
            end

            a_val = [callback_value(cb_data, a[i]) for i in 1:n]

            z_etoile = callback_value(cb_data, z)
            println("current z* = ", z_etoile)
              
            # poids du chemin

            SP2 = Model(CPLEX.Optimizer)
            set_silent(SP2)

            @variable(SP2, delta_2[1:n]>=0)    

            @constraint(SP2,  [i in 1:n], delta_2[i]<=2)
            @constraint(SP2,  sum(delta_2[i] for i in 1:n)<=d2)
            
            set_optimizer_attribute(SP2, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(SP2, "CPXPARAM_TimeLimit", timelimit)
            set_optimizer_attribute(SP2, "CPX_PARAM_THREADS", 1)
            set_optimizer_attribute(SP2, "CPXPARAM_MIP_Display", 4)

            @objective(SP2, Max, sum(a_val[i]*(p[i]+delta_2[i]*ph[i]) for i in 1:n))
            optimize!(SP2)
            delta_2_val = JuMP.value.(delta_2)

            if JuMP.objective_value(SP2) >= S + 1e-5
                cstr = @build_constraint(sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("Ajout d'une contrainte de type 23")
            
            end

            # duree du chemin

            SP1 = Model(CPLEX.Optimizer)
            set_silent(SP1)
            @variable(SP1, delta_1[i in 1:n, j in deltap[i]]>=0)    

            @constraint(SP1,  [i in 1:n, j in deltap[i]], delta_1[i,j]<=d[i,j])

            @constraint(SP1,  sum(delta_1[i,j] for i in 1:n for j in deltap[i])<=d1)
            
            set_optimizer_attribute(SP1, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(SP1, "CPXPARAM_TimeLimit", timelimit)
            set_optimizer_attribute(SP1, "CPX_PARAM_THREADS", 1)
            set_optimizer_attribute(SP1, "CPXPARAM_MIP_Display", 4)

            @objective(SP1, Max, sum(x_val[i,j]*d[i,j]*(1+delta_1[i,j]) for i in 1:n for j in deltap[i]))
            optimize!(SP1)
            delta_1_val = JuMP.value.(delta_1)

            if JuMP.objective_value(SP1) >=  z_etoile + 1e-5
                println("solution SP1 =  ", JuMP.objective_value(SP1))
                cstr2 = @build_constraint(sum(x[i,j]*d[i,j]*(1+delta_1_val[i,j]) for i in 1:n for j in deltap[i]) <= z)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                println("Ajout d'une contrainte de type 24")
            end
        end
    end

    #logfile_name = "plans_coupants.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    #logfile_path = abspath(logfile_name)
    #logfile = open(logfile_path, "w")
    #redirect_stdout(logfile)

    # On précise que le modèle doit utiliser notre fonction de callback
    MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    optimize!(m)

    println(primal_status(m))
    println(termination_status(m))
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end
