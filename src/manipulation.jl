# This file contains operations that modify the presolve problem data

"
--- Row and Column Operations ---
There are four functions for each - add, enque, deque and remove.
There are separate functions for enque and deque intentionally.
There are times we want to remove a row from the active list but not remove it from the problem.
Julia doesnt have NULL. Thus, self-referencing indicates end of the list in either direction.
Initially, the aij pointers are assigned the value nothing. They are updated later by the add_aij functions.
"

# Creates row i with b[i] value b_val and adds it to Presolve_Problem p.
function add_row!(verbose::Bool, p::Presolve_Problem, i::Int, lb::Float64, ub::Float64)
    v = verbose
    v && println("Adding row : $(i)")

    row = Presolve_Row()
    row.i = i
    row.lb = lb
    row.ub = ub
    if(lb == ub)
        row.b_val = lb
    else
        row.b_val = nothing
    end
    row.aij = nothing
    row.is_active = false
    if(i!=1)
        row.prev = p.dictrow[i-1]
        row.prev.next = row
    else
        row.prev = row
        p.rowptr = row
    end
    row.next = row
    enque_row!(v,p,row)
    p.dictrow[i] = row
end

# places the specified row in the active list.
function enque_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("Queueing row : $(row.i)")

    if(row.is_active == false)
        row.is_active = true
        row.active_prev = row.prev
        row.active_next = row
        if(p.rowque.i == -1)
            p.rowque = row
        end
        if(row.active_prev != row)
            row.active_prev.active_next = row
        end
    end
end

# removes the specified row from the active list.
function deque_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("Dequeueing row : $(row.i)")

    if(row.is_active == true)
        row.is_active = false
        if(row.active_prev == row)
            p.rowque = row.active_next
        else
            row.active_prev.active_next = row.active_next
        end
        if(row.active_next != row)
            row.active_next.active_prev = row.active_prev
        end
    end
end

# Creates col j with c[i] value c_val and adds it to Presolve_Problem p.
function add_col!(verbose::Bool, p::Presolve_Problem, j::Int, c_val::Float64, l::Float64, u::Float64)
    v = verbose
    v && println("Adding col : $(j)")

    col = Presolve_Col()
    col.j = j
    col.c_val = c_val
    col.l = l
    col.u = u
    col.aij = nothing
    col.is_independent = false
    if(j!=1)
        col.prev = p.dictcol[j-1]
        col.prev.next = col
    else
        col.prev = col
        p.colptr = col
    end
    col.next = col
    enque_col!(v,p,col)
    p.dictcol[j] = col
end

# places the specified col in the independent list.
function enque_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("Queueing col : $(col.j)")

    if(col.is_independent == false)
        col.ind_prev = col.prev
        col.is_independent = true
        col.ind_next = col
        if(p.colque.j == -1)
            p.colque = col
        end
        if(col.ind_prev != col)
            col.ind_prev.ind_next = col
        end
    end
end

# removes the specified col from the independent list.
function deque_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("Dequeueing col : $(col.j)")

    if(col.is_independent == true)
        col.is_independent = false
        if(col.ind_prev == col)
            p.colque = col.ind_next
        else
            col.ind_prev.ind_next = col.ind_next
        end
        if(col.ind_next != col)
            col.ind_next.ind_prev = col.ind_prev
        end
    end
end

function remove_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    deque_row!(v,p,row)

    while(row.aij != nothing)
        v && println("INSIDE------ and aij is $(row.aij.row.i),$(row.aij.col.j)")
        tmp = row.aij
        key = rc(tmp.row.i,tmp.col.j,p.originaln)
        #enque_col!(p,tmp.col)
        row.aij = tmp.row_next

        if(tmp.col_prev == tmp)
            if(tmp.col_next == tmp)
                tmp.col.aij = nothing
            else
                tmp.col.aij = tmp.col_next
                tmp.col_next.col_prev = tmp.col_next
            end
        else
            if(tmp.col_next == tmp)
                tmp.col_prev.col_next = tmp.col_prev
            else
                tmp.col_prev.col_next = tmp.col_next
            end
        end

        delete!(p.dictaij,key)

        if(row.aij == tmp)
            row.aij = nothing
        end
        tmp = nothing
    end

    if(row.prev == row)
        if(row.next != row)
            p.rowptr = row.next
            row.next.prev = row.next
        else
            p.rowptr = Presolve_Row()
        end
    else
        if(row.next != row)
            row.prev.next = row.next
            row.next.prev = row.prev
        else
            row.prev.next = row.prev
        end
    end

    delete!(p.dictrow,row.i)
end

function remove_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("REMOVE COLUMN CALL")
    deque_col!(v,p,col)

    while(col.aij != nothing)
        tmp = col.aij
        key = rc(tmp.row.i,tmp.col.j,p.originaln)
        #enque_row!(p,tmp.row)
        col.aij = tmp.col_next

        if(tmp.row_prev == tmp)
            if(tmp.row_next == tmp)
                tmp.row.aij = nothing
            else
        #        println("here 1 ")
                tmp.row.aij = tmp.row_next
                tmp.row_next.row_prev = tmp.row_next
            end
        else
            if(tmp.row_next == tmp)
                tmp.row_prev.row_next = tmp.row_prev
            else
                tmp.row_prev.row_next = tmp.row_next
            end
        end

        delete!(p.dictaij,key)

        if(col.aij == tmp)
            col.aij = nothing
        end
        tmp = nothing
    end

    if(col.prev == col)
        if(col.next != col)
            p.colptr = col.next
            col.next.prev = col.next
        else
            p.colptr = Presolve_Col()
        end
    else
        if(col.next != col)
            col.prev.next = col.next
            col.next.prev = col.prev
        else
            col.prev.next = col.prev
        end
    end
    delete!(p.dictcol,col.j)
end

"
--- Matrix Element Operations ---
Matrix elements have only two functions - add_aij_normal and add_aij_transpose
They are never separately removed outside of removing rows or columns.
The input is a sparse matrix stored in the CSC format.
It can efficiently be accessed only in a column major order.
We need our matrix elements to act as two doubly linked lists.
One along the row and one along the column.
add_aij_normal creates the links along the column.
add_aij_transpose creates the links along the row.

We traverse the CSC matrix twice. First in the regular column major order and call add_aij_normal.
The second time we traverse the transpose(A) in column major order and call add_aij_transpose.
Note that the matrix element is already created by the time we call the second function.
"
function add_aij_normal!(verbose::Bool, p::Presolve_Problem, row_id::Int, col_id::Int, row_prev_id::Int, val::Float64)
    v = verbose
    v && println("Adding mat element (normal): $(row_id),$(col_id)")

    aij = Presolve_Matrix()
    aij.row = p.dictrow[row_id]
    aij.col = p.dictcol[col_id]
    aij.val = val
    aij.row_prev = aij
    aij.row_next = aij
    aij.col_next = aij
    aij.col_prev = aij

    if(row_prev_id != -1)
        prev_key = rc(row_prev_id,col_id,p.originaln)
        p.dictaij[prev_key].col_next = aij
        aij.col_prev = p.dictaij[prev_key]
    end
    #can also be done by checking if row_prev == -1
    if(p.dictcol[col_id].aij == nothing)
        p.dictcol[col_id].aij = aij
    end

    key = rc(row_id,col_id,p.originaln)
    p.dictaij[key] = aij
end

function add_aij_transpose!(verbose::Bool, p::Presolve_Problem, row_id::Int, col_id::Int, col_prev_id::Int, val::Float64)
    v = verbose
    v && println("Adding mat element (transpose): $(row_id),$(col_id)")

    aij = p.dictaij[rc(row_id,col_id,p.originaln)]
    if(col_prev_id != -1)
        prev_key = rc(row_id,col_prev_id,p.originaln)
        p.dictaij[prev_key].row_next = aij
        aij.row_prev = p.dictaij[prev_key]
    end

    if(p.dictrow[row_id].aij == nothing)
        p.dictrow[row_id].aij = aij
    end
end

"
--- Presolve Setup and Cleanup ---
make_presolve       : Sets up the linked list connections.
print_info          : Prints all the linked list information. Useful for debugging.
make_new            : Converts the final presolve problem data into the format of the original problem"
function make_presolve!(verbose::Bool, p::Presolve_Problem, A::SparseMatrixCSC{Float64,Int64}, collb::Array{Float64,1},colub::Array{Float64,1}, c::Array{Float64,1}, rowlb::Array{Float64,1},rowub::Array{Float64,1})
    v = verbose
    m,n = size(A)

    # TODO.. make changes everywhere with collb,colub and rowlb,rowub

    # checks to ensure input problem is valid.
    p.originalm != m && error("Wrong size of b wrt A")
    p.originalm != length(rowlb) && error("Wrong size of rowlb wrt A")
    p.originalm != length(rowub) && error("Wrong size of rowub wrt A")
    p.originaln != n && error("Wrong size of c wrt A")
    p.originaln != length(collb) && error("Wrong size of collb wrt A")
    p.originaln != length(colub) && error("Wrong size of colub wrt A")

    v && println("Row SETUP ----- ")
    for i in 1:p.originalm
        add_row!(v,p,i,rowlb[i],rowub[i])
    end

    v && println("COL SETUP -----")
    for j in 1:p.originaln
        add_col!(v,p,j,c[j],collb[j],colub[j])
    end

    # Iterating through the non-zeros of sparse matrix A to construct the dictionary
    Arows = rowvals(A)
    v && println("MAT ELEMENT SETUP -----")
    Avals = nonzeros(A)
    for j = 1:p.originaln
        tmp = -1
        for i in nzrange(A,j)
            r = Arows[i]
            rcval = rc(r,j,p.originaln)
            #dictA[rcval] = Avals[i]
            p.rowcounter[r] += 1
            p.colcounter[j] += 1
            add_aij_normal!(v,p,r,j,tmp,Avals[i])
            tmp = r
        end
    end

    B = A'
    Arows = rowvals(B)
    Avals = nonzeros(B)
    for i = 1:p.originalm
        tmp = -1
        for c in nzrange(B,i)
            j = Arows[c]
            rcval = rc(i,j,p.originaln)
            add_aij_transpose!(v,p,i,j,tmp,Avals[c])
            tmp = j
        end
    end
end

# For Debugging.
function print_info(p::Presolve_Problem)
    println("Row Information--------------------------------------")
    for key in keys(p.dictrow)
        println("-----------")
        row = p.dictrow[key]
        @show row.i
        @show row.lb
        @show row.ub
        @show row.b_val
        if(row.aij != nothing)
            @show row.aij.row.i
            @show row.aij.col.j
            @show row.aij.val
        end
        @show row.prev.i
        @show row.next.i
        @show row.is_active
        @show row.active_prev.i
        @show row.active_next.i
    end

    println("Col Information-------------------------------------------")
    for key in keys(p.dictcol)
        println("-----------")
        col = p.dictcol[key]
        @show col.j
        @show col.c_val
        @show col.l
        @show col.u
        if(col.aij != nothing)
            @show col.aij.row.i
            @show col.aij.col.j
            @show col.aij.val
        end
        @show col.prev.j
        @show col.next.j
        @show col.is_independent
        @show col.ind_prev.j
        @show col.ind_next.j
    end

    println("MAT ELEMENT Information---------------------------")
    for key in keys(p.dictaij)
        println("-----------")
        aij = p.dictaij[key]
        @show aij.row.i
        @show aij.col.j
        @show aij.val

        @show aij.row_prev.row.i, aij.row_prev.col.j, aij.row_prev.val
        @show aij.row_next.row.i, aij.row_next.col.j, aij.row_next.val
        @show aij.col_prev.row.i, aij.col_prev.col.j, aij.col_prev.val
        @show aij.col_next.row.i, aij.col_next.col.j, aij.col_next.val
    end
end

# Constructs the reduced problem
function make_new(verbose::Bool, p::Presolve_Problem)
    v = verbose

    currentn = 0
    newc = Array{Float64,1}()
    newcollb = Array{Float64,1}()
    newcolub = Array{Float64,1}()

    col = p.colptr
    if(col.j != -1)
        v && println("Constructing newc,newlb,newub")
        while(col != nothing)
            v && @show col.j
            push!(newc,col.c_val)
            push!(newlb,col.l)
            push!(newub,col.u)

            currentn = currentn + 1
            p.finalcols[col.j] = currentn

            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    v && @show currentn

    currentm = 0
    newrowlb = Array{Float64,1}()
    newrowub = Array{Float64,1}()
    row = p.rowptr
    if(row.i != -1)
        v && println("Constructing newb")
        while(row != nothing)
            push!(newrowlb,row.lb)
            push!(newrowub,row.ub)
            currentm = currentm + 1
            p.finalrows[row.i] = currentm
            if(row.next == row)
                row = nothing
            else
                row = row.next
            end
        end
    end
    v && @show currentm

    v && println(p.finalcols)
    v && println(p.finalrows)

    I = Array{Int64,1}()
    J = Array{Int64,1}()
    Val = Array{Float64,1}()
    col = p.colptr
    if(col.j != -1)
        v && println("Constructing the new A matrix")
        while(col != nothing)
            tmp = col.aij
            while(tmp != nothing)
                v && @show tmp.row.i , tmp.col.j
                v && @show tmp.val
                push!(J,p.finalcols[tmp.col.j])
                push!(I,p.finalrows[tmp.row.i])
                push!(Val,tmp.val)
                if(tmp.col_next == tmp)
                    tmp = nothing
                else
                    tmp = tmp.col_next
                end
            end
            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    newA = sparse(I,J,Val,currentm,currentn)
    return newA,newcollb,newcolub,newc,newrowlb,newrowub
end
