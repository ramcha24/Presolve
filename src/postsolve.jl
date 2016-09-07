
using MathProgBase

"
--- PostSolving Utilities ---
add_to_stack!       : function that will add the Linear_Dependency element to the stack.
post_solve!         : function that will post solve one Linear_Dependency element.
return_postsolved   : function that will take in the solution from solver for reduced problem and returns solution for original problem
"
function add_to_stack!(l::Linear_Dependency, independentvar::BitArray{1}, pstack::Array{Presolve_Stack,1})
    if (length(l.vec1) != length(l.vec2))
        error("vector1 size not equal to vector 2 size for LD element")
    end
    independentvar[l.index] = false # the variable at this index is not independent anymore
    push!(pstack,l)
end

function post_solve!(post_solvedX::Array{Float64,1}, l::Linear_Dependency)
    post_solvedX[l.index] = l.value

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
        #println("made postsolved at $(l.index) to value $(post_solvedX[l.index])")
    end
end

function return_postsolved(x::Array{Float64,1}, independentvar::BitArray{1}, pstack :: Array{Presolve_Stack,1})
    postsolvedX = zeros(length(independentvar))
    newcols = find(independentvar)

    for i in 1:length(newcols)
        postsolvedX[newcols[i]] = x[i]
    end

    for i in reverse(collect(1:length(pstack)))
        post_solve!(postsolvedX,pstack[i])
    end
    return postsolvedX
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
