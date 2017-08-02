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
                    singleton_row_eq!(v,p,row)
                else
                    forcing_constraints_eq!(v,p,row)
                end
            else
                if(row.lb == -Inf && row.ub == Inf)
                    free_row!(v,p,row)
                else
                    if(row.aij.row_next == row.aij)
                        singleton_row_ineq!(v,p,row)
                    else
                        forcing_constraints_ineq!(v,p,row)
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
    # output notif
    v = verbose
    v && println("EMPTY ROW FOUND AT $(row.i)")

    # feasibility check
    if(row.lb > 0 || row.ub < 0 || row.lb > row.ub)
        error("Empty Row Infeasibility at row $row.i where the bounds are lb - $(row.lb), ub - $(row.ub)")
    else
        # fixing constraint primal y_i and dual (alpha-beta)_i.
        add_to_stack!(Empty_Row(row.i),p.independent_var,p.active_constr,p.pstack)
        # row is redundant now.
        remove_row!(v,p,row)
    end

    # output notif
    v && println("Exiting Empty Row")
end

function free_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    # making lin dep for constraint primal y_i
    v = verbose

    vec1 = Array{Int,1}()
    vec2 = Array{Float64,1}()

    tmp = row.aij
    while(tmp != nothing)
        push!(vec1,tmp.col.j)
        push!(vec2,tmp.aij)
        if(tmp.row_next != tmp)
            tmp = tmp.row_next
        else
            tmp = nothing
        end
    end
    l = Linear_Dependency(row.i,0.0,vec1,vec2)
    # fixing y_i and (alpha-beta)_i
    add_to_stack!(Free_Row(row.i,l),p.independent_var,p.active_constr,p.pstack)
    # row is redundant now.
    remove_row!(v,p,row)
end

function singleton_row_eq!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row, has_bval::Bool)
    v = verbose

    matval = row.aij.val
    col = row.aij.col
    i = row.i
    j = col.j
    y_i = row.lb
    x_j = y_i/matval

    (x_j > col.u || x-j < col.l) && error("Primal Infeasibility in singleton_row $(row.i) with equality")

    add_to_stack!(Singleton_Row_Equality(i,j,x_j,y_i),p.independent_var,p.active_constr,p.pstack)
    remove_row!(v,p,row)
    col.l = col.u = xj
    fixed_col!(v,p,col)
end

function singleton_row_ineq!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row, has_bval::Bool)
    v = verbose
    v && println("SINGlETON ROW FOUND AT $(row.i)")

    matval = row.aij.val
    col = row.aij.col
    i = row.i
    j = col.j

    if(matval > 0)
        l_new = row.lb / matval
        u_new = row.ub / matval
    else
        l_new = row.ub / matval
        u_new = row.lb / matval
    end

    (l_new > col.u || u_new < col.l) && error("Primal Infeasibility in singleton_row $(row.i) with double bound")

    if(l_new <= col.l && col.u <= u_new)
        row_bound = 1
    elseif(l_new <= col.l && col.u > u_new)
        row_bound = 2
    elseif(l_new > col.l && col.u <= u_new)
        row_bound = 3
    else
        row_bound = 4
    end

    vec1 = Array{Int,1}()
    vec2 = Array{Float64,1}()

    tmp = row.aij
    while(tmp != nothing)
        push!(vec1,tmp.col.j)
        push!(vec2,tmp.aij)
        if(tmp.row_next != tmp)
            tmp = tmp.row_next
        else
            tmp = nothing
        end
    end
    l = Linear_Dependency(row.i,0.0,vec1,vec2)

    col.l = max(col.l,l_new)
    col.u = min(col.u,u_new)
    add_to_stack!(Singleton_Row_Inequality(i,j,l,l_new,u_new,row_bound,matval),p.independent_var,p.active_constr,p.pstack)
    remove_row!(v,p,row)
end

function forcing_constraints_eq!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row, has_bval::Bool)
    v = verbose

    tmp = row.aij
    # make g and h. traverse through each aij element in the row.
    g = 0
    h = 0
    while(tmp != nothing)
        if ((g == -Inf || g == +Inf) || (h == -Inf || h == +Inf))
            return
        end
        if(tmp.val < 0)
            g += tmp.val * tmp.col.u
            h += tmp.val * tmp.col.l
        else
            g += tmp.val * tmp.col.l
            h += tmp.val * tmp.col.u
        end

        if(tmp.row_next != tmp)
            tmp = tmp.row_next
        else
            tmp = nothing
        end
    end

    (h < row.lb || g > row.ub) && error("Infeasible problem in forcing constraints procedure")

    # analysis

    col_ind = Array{Int,1}()
    col_bound = Array{Int,1}()
    mat_val = Array{Float64,1}()
    if( g == row.ub )
        v && println("Found lower bound Forcing constraint at row $i")
        #we can fix all variables in this row to their lowerbound or upperbound.
        tmp = row.aij
        while(tmp != nothing)
            col = tmp.col
            push!(col_ind,col.j)
            push!(mat_val,tmp.val)
            if(tmp.val > 0)
                col.u = col.l
                push!(col_bound,-1)
            else
                col.l = col.u
                push!(col_bound,1)
            end
            fixed_col!(v,p,col)

            if(tmp.row_next != tmp)
                tmp = tmp.row_next
            else
                tmp = nothing
            end
        end
        add_to_stack!(Forcing_Row(row.i,col_ind,col_bound,mat_val,g),p.independent_var,p.active_constr,p.pstack)
        remove_row!(v,p,row)
    elseif( h == row.lb )
        v && println("Found upper bound Forcing constraint at row $i")
        #we can fix all variables in this row to their lowerbound or upperbound.
        tmp = row.aij
        while(tmp != nothing)
            col = tmp.col
            push!(col_ind,col.j)
            push!(mat_val,tmp.val)
            if(tmp.val > 0)
                col.l = col.u
                push!(col_bound,1)
            else
                col.u = col.l
                push!(col_bound,-1)
            end
            fixed_col!(v,p,col)
            if(tmp.row_next != tmp)
                tmp = tmp.row_next
            else
                tmp = nothing
            end
        end
        add_to_stack!(Forcing_Row(row.i,col_ind,col_bound,mat_val,h),p.independent_var,p.active_constr,p.pstack)
        remove_row!(v,p,row)
    elseif( g >= row.lb && h <= row.ub )
        # here we have g < ub and lb < h, If g is tighter than lb, update. If h is tighter than ub, update
        row.lb = -Inf
        row.ub = Inf
        free_row!(v,p,row)
    end
end

function empty_col!(v::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("EMPTY COL FOUND AT $(col.j)")

    if((col.c_val > 0 && col.l == -Inf) || (col.c_val < 0 && col.u == Inf))
        error("Dual Infeasibility from empty col at $(col.j)")
    end

    if(col.l == -Inf && col.u == Inf)
        x_j = 0
    elseif (col.l == -Inf)
        x_j = col.u
    elseif (col.u == Inf)
        x_j = col.l
    else
        if col.c_val >= 0
            x_j = col.l
        else
            x_j = col.u
        end
    end

    add_to_stack!(Empty_Col(col.j,x_j,col.c_val),p.independent_var,p.active_constr,p.pstack)
    # TODO.. accounting for constant term in the objective function
    remove_col!(v,p,p.dictcol[col.j])
end

function fixed_col!(v::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("FIXED COL FOUND AT $(col.j)")

    (col.l == -Inf || col.l == Inf) && error("Problem is unbounded due to variable $(col.j) in fixed col")

    vec1 = Array{Int,1}()
    vec2 = Array{Float64,1}()

    # need to substitute in the matrix whenever this variable occurs.
    fixed_val = col.l
    tmp = col.aij
    while(tmp != nothing)
        diff = fixed_val*tmp.val

        if(row.b_val != nothing)
            tmp.row.b_val -= diff
        end
        if(row.lb != -Inf)
            tmp.row.lb -= diff
        end
        if(row.ub != -Inf)
            tmp.row.ub -= diff
        end

        push!(vec1,tmp.row.i)
        push!(vec2,tmp.aij)

        if(tmp.col_next != tmp)
            tmp = tmp.col_next
        else
            tmp = nothing
        end
    end
    add_to_stack!(Fixed_Col(col.j,fixed_val,col.c,Linear_Dependency(col.j,0,vec1,vec2)),p.independent_var,p.active_constr,p.pstack)

    remove_col!(v,p,p.dictcol[col.j])
end


function singleton_col!(v::Bool, p::Presolve_Problem, col::Presolve_Col)
end


# TODO.. MANAGE THE return status to proper ones wrt MPB interface.
