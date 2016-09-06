# Routines and Algorithms that are needed for Presolving.

`
--- Presolver Core ---
empty_row!          : Processes the empty row. Removes it or reports an Infeasibility
presolver!          : Traverses the active list of rows and detects if redundancies are found. Call approporiate functions to handle redundancies.
singleton_row!      : Processes the singleton row. Deletes the row and makes changes to the constraint matrix appropriately
other functions will be added here in the future.
`

function presolver!(verbose::Bool,p::Presolve_Problem,c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
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
            if(row.aij.row_next == row.aij)
                singleton_row!(v,p,row)
            else
                v && println("happy for now")
                #forcing_constraints!(p,row,v)
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

function singleton_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("SINGETONE ROW FOUND AT $(row.i)")

    i = row.aij.row.i
    j = row.aij.col.j
    matval = row.aij.val
    b_val = row.aij.row.b_val

    xj = b_val/matval
    add_to_stack!(Linear_Dependency(j,xj),p.independentvar,p.pstack)
    remove_row!(v,p,row)
    p.activeconstr[row.i] = false
    if(!haskey(p.dictcol,j))
        col = p.dictcol[j]
        col.aij.col_next == col && error("Col Singleton also here")
        error("dictcol key error")
    end
    aij = p.dictcol[j].aij
    while(aij != nothing)
        r = aij.row
        r.b_val -= xj*aij.val
        if(aij.col_next != aij)
            aij = aij.col_next
        else
            aij = nothing
        end
    end
    remove_col!(v,p,p.dictcol[j])
end
