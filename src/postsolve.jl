using MathProgBase

"
--- PostSolving Utilities ---
add_to_stack!           : function that will add one type of a presolve stack element to the stack.
postsolve_stack_element!: function that will post solve one presolve stack element using the corresponding function.
postsolve_lindep        : function that will post solve one Linear_Dependency element.
return_postsolved       : function that will take in the solution from solver for reduced problem and returns solution for original problem
"

# EMPTY ROW
function add_to_stack!(r::Empty_Row,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    active_constr[r.i] = false
    push!(pstack,r)
end

function postsolve_stack_element!(r::Empty_Row, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    constr_primalsol[r.i] = 0.0
    constr_dualsol[r.i] = 0.0
end


# FREE ROW
function add_to_stack!(r::Free_Row,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    active_constr[r.i] = false
    push!(pstack,r)
end

function postsolve_stack_element!(r::Free_Row, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    postsolve_lindep!(constr_primalsol,l)
    constr_dualsol[r.i] = 0.0
end

# SINGLETON ROW EQUALITY
function add_to_stack!(r::Singleton_Row_Equality,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    active_constr[r.i] = false
    push!(pstack,r)
end

function postsolve_stack_element!(r::Singleton_Row_Equality, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    var_primalsol[r.j] = r.x_j
    constr_primalsol[r.i] = r.y_i
    constr_dualsol[r.row] = (var_dualsol[r.j] * r.x_j)/ r.y_i
    var_dualsol[r.j] = 0.0
end

# SINGLETON ROW INEQUALITY
function add_to_stack!(r::Singleton_Row_Inequality,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    #active_constr[r.i] = false
    push!(pstack,r)
end

function postsolve_stack_element!(r::Singleton_Row_Inequality, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    l = r.l
    constr_primalsol[r.i] = l.value
    for j in 1:length(l.vec1)
        constr_primalsol[r.i] += l.vec2[j]*constr_dualsol[l.vec1[j]]
    end

    # col bound
    if(var_primalsol[r.j] = r.l_new)
        col_bound = 1
    elseif (l_new <= var_primalsol[r.j] && var_primalsol[r.j] <= u_new)
        col_bound = 2
    elseif (var_primalsol[r.j] == r.u_new)
        col_bound = 3
    end

    row_bound = r.row_bound

    if (col_bound == 2)
        constr_dualsol[r.i] = 0.0
    end

    if(matval > 0)
        if (col_bound == 1 && (row_bound == 1 || row_bound == 2))
            constr_dualsol[r.i] = 0.0
        elseif (col_bound == 1 && (row_bound == 3 || row_bound == 4))
            constr_dualsol[r.i] = var_dualsol[r.j] / r.matval
        elseif (col_bound == 3 && (row_bound == 1 || row_bound == 3))
            constr_dualsol[r.i] = 0.0
        elseif (col_bound == 3 && (row_bound == 2 || row_bound == 4))
            constr_dualsol[r.i] = var_dualsol[r.j] / r.matval
        end
    else
        if (col_bound == 1 && (row_bound == 1 || row_bound == 2))
            constr_dualsol[r.i] = 0.0
        elseif (col_bound == 1 && (row_bound == 3 || row_bound == 4))
            constr_dualsol[r.i] = -var_dualsol[r.j] / r.matval
        elseif (col_bound == 3 && (row_bound == 1 || row_bound == 3))
            constr_dualsol[r.i] = 0.0
        elseif (col_bound == 3 && (row_bound == 2 || row_bound == 4))
            constr_dualsol[r.i] = -var_dualsol[r.j] / r.matval
        end
    end
end

# Forcing ROW
function add_to_stack!(r::Forcing_Row,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    active_constr[r.i] = false
    push!(pstack,r)
end

function postsolve_stack_element!(r::Forcing_Row, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    # lastcol , val
    lastcol = -1
    val = 0
    for(j in 1:length(r.col_ind))
        dual = var_dualsol[r.col_ind[j]]
        temp = dual/r.mat_val[r.col_ind[j]]

        if((r.col_bound == -1) && (dual < 0 && val < temp))
            lastcol = r.col_ind[j]
            val = temp
        elseif((r.col_bound == 1) && (dual > 0 && val < temp))
            lastcol = r.col_ind[j]
            val = temp
        end
    end

    constr_primalsol[r.i] = r.fixed_val

    if(lastcol == -1)
        constr_dualsol[r.i] = 0.0
    else
        #constr_dualsol[r.i] = var_dualsol[lastcol] /r.mat_val[lastcol]
        constr_dualsol[r.i] = val
        for(j in 1:length(r.col_ind))
            if(j!=lastcol)
                var_dualsol[j] -= (constr_dualsol[r.i])*r.mat_val[j]
            end
        end
        var_dualsol[lastcol] = 0.0
    end

end


# EMPTY COL
function add_to_stack!(c::Empty_Col,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    independentvar[c.j] = false
    push!(pstack,c)
end

function postsolve_stack_element!(c::Empty_Col, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    var_primalsol[c.j] = x_j
    var_dualsol[c.j] = c_j
end

# FIXED COL
function add_to_stack!(c::Fixed_Col,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    independentvar[c.j] = false
    push!(pstack,c)
end

function postsolve_stack_element!(c::Fixed_Col, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    var_primalsol[c.j] = c.x_j

    l = c.l
    var_dualsol[c.j] = l.value
    for i in 1:length(l.vec1)
        var_dualsol[c.j] += l.vec2[i]*constr_dualsol[l.vec1[i]]
        constr_primalsol[l.vec1[i]] += l.vec2[i]*c.x_j
    end
end

"
# SINGLETON COL
function add_to_stack!(c::Fixed_Col,  independentvar::BitArray{1}, active_constr::BitArray{1}, pstack::Array{Presolve_Stack,1})
    independentvar[c.j] = false
    push!(pstack,c)
end

function postsolve_stack_element!(c::Fixed_Col, var_primalsol::Array{Float64,1}, var_dualsol::Array{Float64,1}, constr_primalsol::Array{Float64,1}, constr_dualsol::Array{Float64,1})
    var_primalsol[c.j] = c.x_j

    l = c.l
    var_dualsol[c.j] = l.value
    for i in 1:length(l.vec1)
        var_dualsol[c.j] += l.vec2[i]*constr_dualsol[l.vec1[i]]
        constr_primalsol[l.vec1[i]] += l.vec2[i]*c.x_j
    end
end
"

# General Linear Dependency postsolve
function postsolve_lindep!(post_solvedX::Array{Float64,1}, l::Linear_Dependency)
    post_solvedX[l.index] = l.value

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
        #println("made postsolved at $(l.index) to value $(post_solvedX[l.index])")
    end
end

function return_postsolved(vpsol::Array{Float64,1}, vdsol::Array{Float64,1}, cpsol::Array{Float64,1}, cdsol::Array{Float64,1}, independent_var::BitArray{1}, active_constr::BitArray{1}, pstack :: Array{Presolve_Stack,1})
    n = length(independent_var)
    m = length(active_constr)

    var_primalsol = zeros(n)
    constr_primalsol = zeros(m)
    var_dualsol = zeros(n)
    constr_dualsol = zeros(m)

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


    sol = [var_primalsol,var_dualsol,constr_primalsol,constr_dualsol]

    #NOTE. What is the mapping? from different sols to flags. comment it.
    # 1 - var_primalsol
    # 2 - var_dualsol
    # 3 - constr_primalsol
    # 4 - constr_dualsol

    for i in reverse(collect(1:length(pstack)))
        postsolve_stack_element!(pstack[i],sol[1],sol[2],sol[3],sol[4])
    end
    return sol
end
