using JuMP
using CPLEX


function branchAndCut(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
    """Résout la méthode par plans coupants (sans callback)"""
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
    #set_silent(m)
    
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

    while true
        optimize!(m)

        # solution courant
        x_val = JuMP.value.(x)
        a_val = JuMP.value.(a)
        z_val = JuMP.objective_value(m)

        # sous problème poids du chemin

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
            @constraint(m,  sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
            println("Ajout d'une contrainte de type 23")
        else
            # sous probleme duree du chemin

            SP1 = Model(CPLEX.Optimizer)
            set_silent(SP1)
            @variable(SP1, delta_1[i in 1:n, j in deltap[i]]>=0)    

            @constraint(SP1,  [i in 1:n, j in deltap[i]], delta_1[i,j]<=D[i,j])

            @constraint(SP1,  sum(delta_1[i,j] for i in 1:n for j in deltap[i])<=d1)
            
            set_optimizer_attribute(SP1, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(SP1, "CPXPARAM_TimeLimit", timelimit)
            set_optimizer_attribute(SP1, "CPX_PARAM_THREADS", 1)
            set_optimizer_attribute(SP1, "CPXPARAM_MIP_Display", 4)

            @objective(SP1, Max, sum(x_val[i,j]*d[i,j]*(1+delta_1[i,j]) for i in 1:n for j in deltap[i]))
            optimize!(SP1)
            delta_1_val = JuMP.value.(delta_1)

            if JuMP.objective_value(SP1) >=  z_val + 1e-5
                println("solution SP1 =  ", JuMP.objective_value(SP1))
                @constraint(m,  sum(x[i,j]*d[i,j]*(1+delta_1_val[i,j]) for i in 1:n for j in deltap[i]) <= z)
            else
                break
            end
        end
    end

    logfile_name = "txtFiles/branchAndCut.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    println("objectif = ", JuMP.objective_value(m))

end