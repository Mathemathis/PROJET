using JuMP
using CPLEX


function branchAndCutPLNE(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
    """Résout la méthode par plans coupants (sans callback) avec des PLNE pour la résolution des sous-problèmes"""
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

    logfile_name = "txtFiles/branchAndCutPLNE.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    println("objectif = ", JuMP.objective_value(m))

end

function branchAndCutALG(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
    """Résout la méthode par plans coupants (sans callback) avec une méthode algorithmique pour la résolution des sous-problèmes"""
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

    order_ph=sortperm(ph, rev=true) # ph = p chapeau (incertitude pour le poids de chaque sommet) -> donne les indices dans l'ordre décroissant des ph
    Dh=[D[key] for key in collect(keys(D))] # stocke les distances dans un vecteur
    order_dh=sortperm(Dh, rev=true) # ordre décroissant des distances
    nD=size(order_dh)[1]
    Keys=collect(keys(D))

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

        # solution courante
        x_val = JuMP.value.(x)
        a_val = JuMP.value.(a)
        z_val = JuMP.objective_value(m)

        nodes_visited= [i for i in 1:n if a_val[i]>1e-5]

        # sous problème poids du chemin

        # poids du chemin
        """Idée : on ajoute une contrainte sur le poids des sommets en trouvant un augmentation des poids très élevée"""

        res=sum([a_val[i]*p[i] for i in 1:n])
        capa=0 # budget pour augmenter les delta^2
        I=1 # on regarde les indices dans l'ordre dans l'ordre des p^hat qui permettent d'augmenter le plus le poids du chemin
        delta_2_val = [0 for i in 1:n] # on stocke les delta^2

        while capa+2<=d2 && I<=n # tant que j'ai encore du budget pour augmenter le poids de mes sommets
            if a_val[order_ph[I]]>1e-5 # c'est bien un sommet que j'ai sélectionné 
                capa+=2 # augmentation du budget
                res+=2*ph[order_ph[I]] # augmentation du poids total des sommets
                delta_2_val[order_ph[I]]=2 # stockage valeur delta^2
            end
            I+=1
        end

        while I<=n && a_val[order_ph[I]]<1e-5 # on passe les sommets restants non selectionnés
            I+=1
        end

        if I<=n # on augmente avec le budget qu'il nous reste
            res+=(d2-capa)*ph[order_ph[I]] 
            delta_2_val[order_ph[I]]=d2-capa
        end

        if res >= S + 1e-5 # si notre poids est plus grand que S -> ajout de la contrainte
            @constraint(m,  sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
            println("Ajout d'une contrainte de type 23")
        else
            # sous probleme duree du chemin
            """Même idée sur les chemins : on augmente le poids d'un chemin en trouvant les arcs qui augmente le plus le coût"""

            res=sum([x_val[key]*d[key] for key in Keys])
            capa=0
            I=1
            delta_1_val = Dict(key => 0.0 for key in collect(keys(D))) # on initialise tous les delta^1 des arêtes à 0
            while capa+D[Keys[order_dh[I]]]<=d1 && I<nD # exactement même idée que pour les sommets
                if x_val[Keys[order_dh[I]]]>1e-5
                    capa+=D[Keys[order_dh[I]]]
                    res+=d[Keys[order_dh[I]]]*D[Keys[order_dh[I]]]
                    delta_1_val[Keys[order_dh[I]]]=D[Keys[order_dh[I]]]
                end
                I+=1
            end
            while I<=nD && x_val[Keys[order_dh[I]]]<1e-5
                I+=1
            end
            if I<=nD
                res+=d[Keys[order_dh[I]]]*(d1-capa)
                delta_1_val[Keys[order_dh[I]]]=d1-capa
            end

            if res >=  z_val + 1e-5
                @constraint(m,  sum(x[key]*d[key]*(1+delta_1_val[key]) for key in Keys) <= z)
            else
                break
            end
        end
    end

    logfile_name = "txtFiles/branchAndCutALG.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    println("objectif = ", JuMP.objective_value(m))

end