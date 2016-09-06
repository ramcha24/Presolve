import MathProgBase
export can_solve_mip, can_solve_socp, can_solve_sdp, can_solve_exp
export set_default_solver, get_default_solver

function set_default_solver(solver::MathProgBase.AbstractMathProgSolver)
  global DEFAULT_SOLVER
  DEFAULT_SOLVER = solver
end

function get_default_solver()
  if DEFAULT_SOLVER == nothing
    error("The default solver is set to `nothing`
         You must have at least one solver installed to use Convex.
         You can install a solver such as SCS by running:
         Pkg.add(\"SCS\").
         You will have to restart Julia after that.")
  end
  return DEFAULT_SOLVER
end

solvers = [("SCS", "SCSSolver"), ("ECOS", "ECOSSolver"), ("Gurobi", "GurobiSolver"), ("Mosek", "MosekSolver"),
          ("GLPKMathProgInterface", "GLPKSolverMIP")]

for (dir, solver) in solvers
  if isdir(Pkg.dir(dir)) && DEFAULT_SOLVER == nothing
    eval(parse("using "*dir))
    eval(parse("set_default_solver("*solver*"())"))
  end
end


if get_default_solver() == nothing
  packages = ""
  for (dir, solver) in solvers
    packages = packages*dir*" | "
  end
  warn("***********************************************************************************************
       You don't have any of
       "*packages*" installed.
       You must have at least one of these solvers. You can install a solver such as SCS by running:
       Pkg.add(\"SCS\")
       You will have to restart Julia after that.
       ***********************************************************************************************")
end
