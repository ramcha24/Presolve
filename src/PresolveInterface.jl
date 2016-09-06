importall MathProgBase.SolverInterface

export Presolver

type Presolver <: AbstractMathProgSolver
    realsolver :: AbstractMathProgSolver
    # Presolver() gives an error. says get_defualt_solver undefined, why?
    function Presolver()
        realsolver = get_defualt_solver()
        new(realsolver)
    end
    function Presolver(solver::AbstractMathProgSolver)
        realsolver = solver
        new(realsolver)
    end
end

type PresolveMathProgModel <: AbstractMathProgModel
    m::Int
    n::Int
    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    c::Vector{Float64}
    p::Presolve_Problem
    orig_sense::Symbol
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Array{Float64,1}
    innermodel :: AbstractMathProgModel
    realsolver #gives an error if I specify the type. Can't converr from SCSSolver to AbstractMathProgSolver
    function PresolveMathProgModel()
        realsolver=get_default_solver()
        innermodel = LinearQuadraticModel(realsolver)
        p = Presolve_Problem()
        new(0,0,spzeros(0,0),Int[],Int[],:Min,p,:NotSolved,0.0,Int[],model,realsolver)
    end
    function PresolveMathProgModel(solver::AbstractMathProgSolver)
        realsolver=solver
        innermodel = LinearQuadraticModel(realsolver)
        p = Presolve_Problem()
        new(0,0,spzeros(0,0),Int[],Int[],:Min,p,:NotSolved,0.0,Int[],model,realsolver)
    end
end

function ConicModel(s::Presolver)
    return PresolveMathProgModel(s.realsolver)
end

LinearQuadraticModel(s::Presolver) = ConicToLPQPBridge(ConicModel(s))

status(m::PresolveMathProgModel) = m.solve_stat
getobjval(m::PresolveMathProgModel) = m.obj_val
getsolution(m::PresolveMathProgModel) = copy(m.primal_sol)
numvar(model::PresolveMathProgModel) = model.m  # my m doesn't change. in future it might.
numconstr(model::PresolveMathProgModel) = model.n
supportedcones(s::Presolver) = [:Free, :Zero, :NonNeg, :NonPos]

loadproblem!(model::PresolveMathProgModel,c,A,b,constr_cones,var_cones)
    = loadproblem!(model,c,sparse(A),b,constr_cones,var_cones)

function loadproblem!(model::PresolveMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    bad_cones = [:SOC,:SOCRotated,:SDP,:ExpPrimal,:ExpDual]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) is currently not supported")
        cone_vars[1] == :Free && error("There is a free constraint")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) is currently not supported")
    end

    # going to be used later. c_cones, and v_cones?
    cone_dict = Dict{Symbol,Int64}()
    cone_dict[:Zero] = 1
    cone_dict[:NonNeg] = 2
    cone_dict[:NonPos] = 3
    cone_dict[:Free] = 4

    c_cones = (Array{Int,1}(),Array{Int,1}(),Array{Int,1}(),Array{Int,1}())
    for (cone,id) in constr_cones
        if(isa(id,Number))
           push!(c_cones[cone_dict[c]],id)
        else
           append!(c_cones[cone_dict[c]],id)
        end
        push!(c_cones[cone_dict[cone]],id)
    end

    v_cones = (Array{Int,1}(),Array{Int,1}(),Array{Int,1}(),Array{Int,1}())
    for (cone,id) in var_cones
        if(isa(id,Number))
           push!(v_cones[cone_dict[c]],id)
        else
           append!(v_cones[cone_dict[c]],id)
        end
        push!(v_cones[cone_dict[cone]],id)
    end

    # Column or Variable bounds
    collb = fill(-Inf, length(c))
    colub = fill(Inf, length(c))
    for (cone,idxs) in var_cones
        l = (cone == :Free || cone == :NonPos) ? -Inf : 0.0
        u = (cone == :Free || cone == :NonNeg) ?  Inf : 0.0
        for idx in idxs
            collb[idx] = l
            colub[idx] = u
        end
    end

    # Row or Constraint bounds
    rowlb = Array(Float64,length(b))
    rowub = Array(Float64,length(b))
    for (cone,idxs) in constr_cones
        for idx in idxs
            rowlb[idx] = (cone == :Zero || cone == :NonPos) ? b[idx] : -Inf
            rowub[idx] = (cone == :Zero || cone == :NonNeg) ? b[idx] :  Inf
        end
    end

    model.m,model.n = size(A)
    model.A = A
    model.b = b
    model.c = c

    #TODO.. change the variable orders in presolve routines.

    make_presolve!(false,model.p,A,collb,colub,c,rowlb,rowub)

    newA,newcollb,newcolub,newc,newrowlb,newrowub = presolver!(false,model.p,A,collb,colub,c,rowlb,rowub)

    loadproblem!(model.innermodel,A,collb,colub,c,rowlb,rowub,model.orig_sense)

    return model
end

function optimize!(m::PresolveMathProgModel)
    ans = Array{Float64,1}()
    if(length(find(m.p.independentvar))!=0)
        optimize!(m.innermodel)
        if(m.innermodel.solve_stat!= :Optimal)
            error("Problem status is not optimal, solution might be inaccurate")
        end
        m.solve_stat = m.innermodel.solve_stat
        ans = m.innermodel.getsolution()
    end
    m.primal_sol = return_postsolved(ans,m.p.independentvar,m.p.pstack)
    m.obj_val = dot(m.c, m.primal_sol)
end


# TODO: getdual and getvardaul, can be done once add the dualsolution support for all the presolve routines.
