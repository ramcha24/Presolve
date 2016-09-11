# Routines and Algorithms that are needed for Presolving.

"
--- Presolver Core ---
empty_row!          : Processes the empty row. Removes it or reports an Infeasibility
presolver!          : Traverses the active list of rows and detects if redundancies are found. Call approporiate functions to handle redundancies.
singleton_row!      : Processes the singleton row. Deletes the row and makes changes to the constraint matrix appropriately
other functions will be added here in the future."
function presolver!(verbose::Bool, p::Presolve_Problem, A::SparseMatrixCSC{Float64,Int64}, collb::Array{Float64,1},colub::Array{Float64,1}, c::Array{Float64,1}, rowlb::Array{Float64,1},rowub::Array{Float64,1})
    v = verbose
    v && println("PRESOLVE ROUTINES...............................")
    row = Presolve_Row()
    col = Presolve_Col()
    tmp = p.rowque

    while(tmp != nothing)
        row = tmp
        deque_row!(v,p,row)
        if(row.aij == nothing)
            empty_row!(v,p,row)
        else
            if(row.lb == row.ub)
                if(row.aij.row_next == row.aij)
                    singleton_row!(v,p,row,true)
                else
                    forcing_constraints(v,p,row,true)
                end
            else
                if(row.lb == -Inf && row.ub == Inf)
                    free_row!(v,p,row)
                else
                    if(row.aij.row_next == row.aij)
                        singleton_row!(v,p,row,false)
                    else
                        forcing_constraints(v,p,row,false)
                    end
                end
            end
        end
        if(tmp.next == tmp)
            tmp = nothing
        else
            tmp = tmp.next
        end
    end

    @show p.colque.j
    tmp = p.colque
    while(tmp != nothing)
        col = tmp
        #@show col.j
        deque_col!(v,p,col)
        if(col.aij == nothing)
            empty_col!(v,p,col)
        else
            if(col.l == col.u)
                fixed_col!(v,p,col)
            elseif(col.aij.col_next == col.aij)
                singleton_col!(v,p,row,true)
            end
        end
        if(tmp.next == tmp)
            tmp = nothing
        else
            tmp = tmp.next
        end
    end

    v && println("Making the reduced Problem")
    A,collb,colub,c,rowlb,rowub = make_new(v,p)
    return A,collb,colub,c,rowlb,rowub
end

function empty_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("EMPTY ROW FOUND AT $(row.i)")

    if(!roughly(row.b_val,0.0))
        error("Empty Row Infeasibility at row $row.i and b[i] is - $(row.b_val)")
    else
        remove_row!(v,p,row)
        p.activeconstr[row.i] = false
    end

    v && println("Exiting Empty Row")
end

function free_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    # TODO.. dual variable
    remove_row!(v,p,row)
end

function singleton_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row, has_bval::Bool)
    v = verbose
    v && println("SINGETON ROW FOUND AT $(row.i)")

    # TODO.. make the changes for has_bval to be true and for it to be false

    col = row.aij.col
    i = row.i
    j = col.j
    matval = row.aij.val

    if(has_bval)
        b_val = row.aij.row.b_val
        xj = b_val/matval

        p.activeconstr[row.i] = false
        col.l = col.u = xj

        fixed_col!(v,p,col)
        remove_row!(v,p,row)
    else
        #there is only one variable and it is of the form , lb <= aij*xj <= ub , check if it can be made tighter.
        l_val = (matval > 0) ? row.lb/matval : row.ub/matval
        u_val = (matval > 0) ? row.ub/matval : row.lb/matval

        (l_val > col.u || u_val < col.l) && error("Infeasibility in singleton_row $(row.i) with double bound")

        col.l = max(col.l,l_val)
        col.u = min(col.u,u_val)
        if(col.l == col.u)
            fixed_col!(v,p,col)
        end
        free_row!(v,p,row)
    end
end

function forcing_constraints!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row, has_bval::Bool)
    tmp = row.aij
    # make g and h. traverse through each aij element in the row.
    g = 0
    while(tmp != nothing)
        if(g == -Inf || g == +Inf)
            break
        end
        if(tmp.val < 0)
            g += tmp.val * tmp.col.u
        else
            g += tmp.val * tmp.col.u
        end

        if(tmp.row_next != tmp)
            tmp = tmp.row_next
        else
            tmp = nothing
        end
    end

    tmp = row.aij
    h = 0
    while(tmp != nothing)
        if(h == -Inf || h == +Inf)
            break
        end
        if(tmp.val < 0)
            h += tmp.val * tmp.col.l
        else
            h += tmp.val * tmp.col.l
        end

        if(tmp.row_next != tmp)
            tmp = tmp.row_next
        else
            tmp = nothing
        end
    end

    (h < row.lb || g > row.ub) && error("Infeasible problem in forcing constraints procedure")

    if(has_bval)
        b_val = row.b_val
        # analysis
        if(is_zero(g-b_val))
            v && println("Found lower bound Forcing constraint at row $i")
            #we can fix all variables in this row to their lowerbound or upperbound.
            tmp = row.aij
            while(tmp != nothing)
                col = tmp.col
                if(tmp.val > 0)
                    col.u = col.l
                else
                    col.l = col.u
                end
                fixed_col!(v,p,col)

                if(tmp.row_next != tmp)
                    tmp = tmp.row_next
                else
                    tmp = nothing
                end
            end
            remove_row!(v,p,row)
        elseif(is_zero(h-b_val))
            v && println("Found upper bound Forcing constraint at row $i")
            #we can fix all variables in this row to their lowerbound or upperbound.
            tmp = row.aij
            while(tmp != nothing)
                col = tmp.col
                if(tmp.val > 0)
                    col.l = col.u
                else
                    col.u = col.l
                end
                fixed_col!(v,p,col)

                if(tmp.row_next != tmp)
                    tmp = tmp.row_next
                else
                    tmp = nothing
                end
            end
            remove_row!(v,p,row)
        end
    else
        # here we have g <= ub and lb <= h, If g is tighter than lb, update. If h is tighter than ub, update
        if(g > row.lb)
            row.lb = g
        end

        if(h < row.ub)
            row.ub = h
        end
    end
end

function empty_col!(v::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("EMPTY COL FOUND AT $(col.j)")

    if(col.c_val == 0.0)
        # nothing to be done here. the variable has no effect on the problem. can just remove the column. fix it to value 0 (default)
        add_to_stack!(LinearDependency(col.j,0),p.independentvar,p.pstack)
    elseif (col.c_val > 0)
        col.l == -Inf && error("Problem is unbounded due to variable $(col.j)")
        add_to_stack!(LinearDependency(col.j,col.l),p.independentvar,p.pstack)
    else
        col.u == Inf && error("Problem is unbounded due to variable $(col.j)")
        add_to_stack!(LinearDependency(col.j,col.u),p.independentvar,p.pstack)
    end
    remove_col!(v,p,p.dictcol[col.j])
end

function fixed_col!(v::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("FIXED COL FOUND AT $(col.j)")

    (col.l == -Inf || col.l == Inf) && error("Problem is unbounded due to variable $(col.j) in fixed col")

    # need to substitute in the matrix whenever this variable occurs.
    f_val = col.l
    tmp = col.aij
    while(tmp != nothing)
        diff = f_val*tmp.val
        tmp.row.b_val -= diff
        tmp.row.lb -= diff
        tmp.row.ub -= diff
        if(tmp.col_next != tmp)
            tmp = tmp.col_next
        else
            tmp = nothing
        end
    end
    add_to_stack!(LinearDependency(col.j,col.l),p.independentvar,p.pstack)
    remove_col!(v,p,p.dictcol[col.j])
end

function singleton_col!(p::Presolve_Problem, v::Bool)

end


# TODO.. MANAGE THE return status to proper ones wrt MPB interface.
