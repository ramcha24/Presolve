importall MathProgBase.SolverInterface

export Presolver

immutable Presolver <: AbstractMathProgSolver
    realsolver :: AbstractMathProgSolver
end

Presolver(;kwargs...) = Presolver(kwargs)


type PresolveMathProgModel <: AbstractMathProgModel
    #variables that will be needed
    m::Int
    n::Int
    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    c::Vector{Float64}
    p::Presolve_Problem
    #realsolver::AbstractMathProgSolver
    orig_sense::Symbol
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Float64
    realsolver :: AbstractMath
end

# I have to neccesitate one keyword argument no arguemnt implies I take the default solver.

PresolveMathProgModel(;kwargs...) = PresolveMathProgModel(0,0,spzeros(0,0),Int[],Int[],Presolve_Problem(false,0,0),
:Min,:0,0.0,Int[],kwargs)

PresolveMathProgModel() = PresolveMathProgModel(0,0,spzeros(0,0),Int[],Int[],Presolve_Problem(false,0,0),
:Min,:0,0.0,Int[],AbstractMathProgSolver())

# i have to resolve how I take in the options to the real solver.

ConicModel(s::Presolver) = PresolveMathProgModel(;s.options...)
LinearQuadraticModel(s::Presolver) = ConicToLPQPBridge(ConicModel(s))

function optimize!(m::PresolveMathProgModel)
    solution = options.

    m.presolve_sol = solution
end

status(m::PresolveMathProgModel) = m.solve_stat
getobjval(m::PresolveMathProgModel) = m.obj_val
getsolution(m::PresolveMathProgModel) = copy(m.primal_sol)

#orderconesforscs

loadproblem!(model::PresolveMathProgModel,c,A,b,constr_cones,var_cones)
    = loadproblem!(model,c,sparse(A),b,constr_cones,var_cones)

function loadproblem!(model::PresolveMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    bad_cones = [:SOC,:SOCRotated,:SDP,:ExpPrimal,:ExpDual]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) is currently not supported")
    end
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) is currently not supported")
    end

    c_cones = [(cone,[idxs..]) for (cone,idxs) in constr_cones]
    v_cones = [(cone,[idxs..]) for (cone,idxs) in var_cones]

    #here I need to generate col lb, col ub, row lb, row ub and then make my presolve problem

    model.m,model.n = size(A)
    model.A = A
    model.b = b
    model.c = c

    make_presolve!(false,model.p,c,A,b,collb,colub,rowlb,rowub)

    # this is where the conic model is built.
    # fill in the model variables.

    return model
end

numvar(model::PresolveMathProgModel) = model.input_numvar
numconstr(model::PresolveMathProgModel) = model.input_numconstr

supportedcones(s::Presolver) = [:Free, :Zero, :NonNeg, :NonPos]

# TODO: getdual and getvardaul
