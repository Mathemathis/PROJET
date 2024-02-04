using JuMP
using CPLEX
include("PLNE_compacte.jl")
include("./utils/utils_heuristic.jl")
include("./utils/parsing.jl")

function plne_dual(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, name_instance, formulation_agregee, with_initial_values, timelimit)
    """Résout le PLNE par la méthode duale (question 4)"""

    # initialisation des voisinages entrants et sortants pour chaque sommet
    if with_initial_values=="with initial values"
        Keys=collect(keys(D))
        d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S, 100000)
        p_set, ph_set, Kp_set, Kph_set, A_p_set, A_ph_set=calcul_p_ph_k_set(p, ph, S)
        # initial_values=plne_dual_poids(n, s, t, S, d1, d2, p, ph, d, D)
        initial_values=get_init_sol(name_instance)

        borne_sup=sum([initial_values[1][key]*d[key] for key in Keys])
        capa=0
        for d_D in reverse(d_D_set)
            if capa>d1-1e-5
                break
            end
            for (i,j) in A_d_D_set[d_D]
                if initial_values[1][i,j]>1e-5
                    if capa+D[i,j]<=d1
                        capa+=D[i,j]
                        borne_sup+=d[i,j]*D[i,j]
                    else
                        borne_sup+=d[i,j]*(d1-capa)
                        capa=d1
                        break
                    end
                end
            end
        end
        # println(borne_sup)
        n, s, t, S, d1, d2, p, ph, d, D, initial_values = simplify_instance(name_instance, borne_sup, initial_values)
    end

    deltap, deltam=initDelta(d, n)
    Keys=collect(keys(D))
    d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S, borne_sup)
    p_set, ph_set, Kp_set, Kph_set, A_p_set, A_ph_set=calcul_p_ph_k_set(p, ph, S)

    m = Model(CPLEX.Optimizer)

    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    # @variable(m, a[1:n]>=0)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    #Dualisation Contrainte poids


    if formulation_agregee=="formulation agregee"

        # @variable(m, wph[ph_i in ph_set, k in 1:Kph_set[ph_i]], Bin)
        # @constraint(m,  [ph_i in ph_set], sum(k*wph[ph_i,k] for k in 1:Kph_set[ph_i]) ==sum(a[i] for i in A_ph_set[ph_i]))
        # @constraint(m,  [ph_i in ph_set], sum(wph[ph_i,k] for k in 1:Kph_set[ph_i]) <=1)

        # @variable(m, gamma>=0)
        # @variable(m, eta[ph_i in ph_set, k in 1:Kph_set[ph_i]]>=0)
    
        # @constraint(m,  [ph_i in ph_set, k in 1:Kph_set[ph_i]], gamma+eta[ph_i]>=wph[ph_i]*ph_i)
        # @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(k*eta[ph_i,k] for ph_i in ph_set for k in 1:Kph_set[ph_i])+d2*gamma<=S)

        @variable(m, gamma>=0)
        @variable(m, eta[i in 1:n]>=0)
    
        @constraint(m,  [i in 1:n], gamma+eta[i]>=a[i]*ph[i])
        @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma<=S)

        @variable(m, y[d_D in d_D_set, k in 1:K_set[d_D]], Bin)
        @constraint(m,  [d_D in d_D_set], sum(k*y[d_D,k] for k in 1:K_set[d_D]) ==sum(x[i,j] for (i,j) in A_d_D_set[d_D]))
        @constraint(m,  [d_D in d_D_set], sum(y[d_D,k] for k in 1:K_set[d_D]) <=1)

        @variable(m, alpha>=0)
        @variable(m, beta[d_D in d_D_set, k in 1:K_set[d_D]]>=0)

        @constraint(m,  [d_D in d_D_set, k in 1:K_set[d_D]], alpha+beta[d_D, k]>=d_D[1]*y[d_D,k])

        @objective(m, Min, sum(k*d_D[1]*y[d_D,k] for d_D in d_D_set for k in 1:K_set[d_D])+alpha*d1+sum(beta[d_D,k]*k*d_D[2] for d_D in d_D_set for k in 1:K_set[d_D]))
    else
        @variable(m, gamma>=0)
        @variable(m, eta[i in 1:n]>=0)
    
        @constraint(m,  [i in 1:n], gamma+eta[i]>=a[i]*ph[i])
        @constraint(m,  sum(a[i]*p[i] for i in 1:n)+2*sum(eta[i] for i in 1:n)+d2*gamma<=S)

        @variable(m, alpha>=0)
        @variable(m, beta[i in 1:n, j in deltap[i]]>=0)

        @constraint(m,  [i in 1:n, j in deltap[i]], alpha+beta[i,j]>=x[i,j]*d[i,j])
        @objective(m, Min, sum(x[i,j]*d[i,j] for i in 1:n for j in deltap[i])+alpha*d1+sum(beta[i,j]*D[i,j] for i in 1:n for j in deltap[i]))
    end

    x_last=[]
    best_integer=-1
    i=0

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
    
            CPLEX.load_callback_variable_primal(cb_data, context_id) # on regardera les valeurs des variables dans la solution courante

            x_val=Dict(key => callback_value(cb_data, x[key]) for key in Keys) # on stocke x
            beta_val=Dict((d_D,k) => callback_value(cb_data, beta[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D]) # on stocke x
            a_val = [callback_value(cb_data, a[i]) for i in 1:n]
            alpha_val=callback_value(cb_data, alpha)
            nodes_visited= [i for i in 1:n if a_val[i]>1e-5]
            edges_visited= [key for key in Keys if x_val[key]>1e-5]
            # println(edges_visited, x_last)
            obj_value=sum(x_val[key]*d[key] for key in Keys)+alpha_val*d1+sum(beta_val[d_D,k]*d_D[2]*k for d_D in d_D_set for k in 1:K_set[d_D])
            if edges_visited!=x_last
                # println(nodes_visited, " ", obj_value, " ", best_integer)
                y_val=Dict((d_D,k) => callback_value(cb_data, y[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D])
                K_set_star=Dict(d_D => k for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)
                cstr3 = @build_constraint(sum(sum(y[d_D,k] for k in K_set_star[d_D]:K_set[d_D]) for d_D in keys(K_set_star))<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-sum(x_val[i,j]-x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr3)

                vars= vcat([x[i,j] for i in 1:n for j in deltap[i]], [a[i] for i in 1:n], [y[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D]],
                            [alpha], [beta[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D]], [gamma], [eta[i] for i in 1:n])
                
                vals= vcat([callback_value(cb_data,x[i,j]) for i in 1:n for j in deltap[i]], [callback_value(cb_data,a[i]) for i in 1:n], [callback_value(cb_data,y[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D]],
                            [callback_value(cb_data,alpha)], [callback_value(cb_data,beta[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D]], [callback_value(cb_data,gamma)], [callback_value(cb_data,eta[i]) for i in 1:n])
                
                MOI.submit(m, MOI.HeuristicSolution(cb_data), vars, vals)
                x_last=edges_visited
                if best_integer>=0
                    best_integer=min(best_integer, obj_value)
                else
                    best_integer=obj_value
                end
            end
        end
    end

    if with_initial_values=="with initial values"
        JuMP.set_start_value.(x, JuMP.Containers.SparseAxisArray(initial_values[1]))
        JuMP.set_start_value.(a, initial_values[2])
    end
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_HeuristicEffort", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)
    
    if formulation_agregee=="formulation agregee"
        MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    end

    if formulation_agregee=="formulation agregee"
        if with_initial_values=="with initial values"
            logfile_name = "txtFiles/plne_dual/formulation_agregee/with_initial_values/$name_instance.txt"
        else
            logfile_name = "txtFiles/plne_dual/formulation_agregee/without_initial_values/$name_instance.txt"
        end
    else
        if with_initial_values=="with initial values"
            logfile_name = "txtFiles/plne_dual/formulation_non_agregee/with_initial_values/$name_instance.txt"
        else
            logfile_name = "txtFiles/plne_dual/formulation_non_agregee/without_initial_values/$name_instance.txt"
        end
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    # Résolution d’un modèle
    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    println(primal_status(m))
    println(termination_status(m))

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        close(logfile)
        return JuMP.value.(x), JuMP.value.(a), JuMP.value.(eta), JuMP.value.(gamma)
    end
    close(logfile)
end