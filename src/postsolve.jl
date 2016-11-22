using MathProgBase

"
--- PostSolving Utilities ---
add_to_stack!       : function that will add the Linear_Dependency element to the stack.
post_solve!         : function that will post solve one Linear_Dependency element.
return_postsolved   : function that will take in the solution from solver for reduced problem and returns solution for original problem
"

function add_to_stack!(l::Linear_Dependency, independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    if (length(l.vec1) != length(l.vec2))
        error("vector1 size not equal to vector 2 size for LD element")
    end
    if (l.flag == 1 || l.flag == 2)
        independentvar[l.index] = false # the variable at this index is not independent anymore
    elseif (l.flag == 3 || l.flag == 4)
        active_constr[l.index] = false # ??

    push!(pstack,l)
end

function post_solve!(post_solvedX::Array{Float64,1}, l::Linear_Dependency)
    post_solvedX[l.index] = l.value

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
        #println("made postsolved at $(l.index) to value $(post_solvedX[l.index])")
    end
end

function return_postsolved(cpsol::Array{Float64,1},cdsol::Array{Float64,1},vpsol::Array{Float64,1},vdsol::Array{Float64,1}, independent_var::BitArray{1}, active_constr::BitArray{1}, pstack :: Array{Presolve_Stack,1})
    n = length(independent_var)
    m = length(active_constr)

    var_primalsol = zeros(n)
    var_dualsol = zeros(n)
    constr_primalsol = zeros(m)
    constr_dualsol = zeros(m)

    sol = [var_primalsol,var_dualsol,constr_primalsol,constr_dualsol]

    newcols = find(independent_var)
    newrows = find(active_constr)

    for i in 1:length(newcols)
        var_primalsol[newcols[i]] = vpsol[i]
        var_dualsol[newcols[i]] = vdsol[i]
    end

    for i in 1:length(newrows)
        constr_primalsol[newrows[i]] = cpsol[i]
        constr_dualsol[newrows[i]] = cdsol[i]
    end

    for i in reverse(collect(1:length(pstack)))
        post_solve!(sol[pstack[i].flag],pstack[i])
    end
    return sol
end

"
--- Miscellaneous Utilities ---
"
function is_zero(i::Float64)
    if(abs(i-0.0) <= 1e-3)
        return true
    else
        return false
    end
end

function is_equal(a::Array{Float64,1}, b::Array{Float64,1})
    (length(a) != length(b)) && error("trying to determine equality of arrays of different sizes")
    for i in 1:length(a)
        if(is_zero(a[i]-b[i]) == false)
            return false
        end
    end
    return true
end

function is_lb_Unbounded(lb::Float64)
    if(lb == -Inf)
        return true
    else
        return false
    end
end


function is_ub_Unbounded(ub::Float64)
    if(ub == +Inf)
        return true
    else
        return false
    end
end

# generates the unique key from row,col index for creating the dictionary
function rc(x::Int, y::Int, M::Int)
    return (x-1)*M + y
end

function roughly(x::Float64, y::Float64)
    if(abs(x-y) < 1e-3)
        return true
    else
        return false
    end
end

# TODO.. cleanup these utility functions, why is there one of zero and one for roughly. make it into one function. no need unbounded functions.
