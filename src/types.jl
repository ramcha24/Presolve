export Presolve_Element, Presolve_Row, Presolve_Col, Presolve_Matrix
export Presolve_Stack, Linear_Dependency
export Presolve_Problem

`
--- Introduction ---
The optimization problem being solved is of the form -
    min c*x
    s.t A*x = b
        lb <= x <= ub

    where x, c, lb, ub are n-dimensional vectors. b is a m-dimensional vector.
    A is a m x n constraint co-efficient matrix (usually spase)

The information we store about the problem is logically divided into - rows (m-dims), columns (n-dims) and constraint matrix (m x n dims)

Presolve_Element is the overall abstract type for storing the information about the LP problem.
The three subtypes are Presolve_Row , Presolve_Col , Presolve_Matrix
`

`
--- Active Rows and Independent Columns ---
The Presolve_Row and Presolve_Col data types defined below holds two doubly linked lists within it.
First - the normal previous and next references for rows(cols) of constraint matrix.
Second - a doubly linked list that signifies "active" rows ("independent cols"). This is explained here.

Active row means that row constraint is still to be processed.
Independent col means that corresponding variabls xj is still independent of other variables and we havent tested it yet.

Initially all rows are made "active" and all columns are made "Independent". The following is a description for rows. A similar workflow is adopted for columns.
In the function presolver! every row which is initially active is traversed and we first "deactivate" or "deque" them
After they are dequed from the active queue, they are then analyzed for possible redundancie.
If we are able to detect a redundancy we will appropriately remove the row.
If no redundancies can be detected, they have only been dequed and can still be accessed by the normal prev and next references.
At the end of the presolver! call, there are no active rows.
We make the new optimization problem by using the prev and next references among the rows and columns that havent been deleted yet.
`

abstract Presolve_Element

type Presolve_Row <: Presolve_Element
    # Each row represents a constraint row i
    i :: Int                                # Row number in the original problem.
    lb :: Float64                        # row[i] lower bound in the original problem.
    ub :: Float64                        # row[i] upper bound in the original problem.
    b_val :: Union{Float64,Nothing}                        # b[i] in the original problem.
    aij :: Union{Presolve_Element,Nothing}  # Reference to the first non-zero matrix entry in row i. aij = nothing indicates empty row.
    prev :: Presolve_Row                    # Reference to previous row in the original problem. Row 5 holds ref to Row 4 etc
    next :: Presolve_Row                    # Reference to next row in the original problem. Row 5 holds ref to Row 6 etc
    is_active :: Bool                       # true if the row is in the active doubly linked list. false otherwise.
    active_prev :: Presolve_Row             # if in the active doubly linked list, reference to the previous row in the active doubly linked list.
    active_next :: Presolve_Row             # if in the active doubly linked list, reference to the next row in the active doubly linked list.
    function Presolve_Row()
        # Constructor that creates a default row which has an invalid row.i value.
        n = new()
        n.i = -1
        n.rowlb = 0.0
        n.rowub = 0.0
        n.b_val = 0.0
        n.is_active = true
        n.aij = nothing
        n.prev = n
        n.next = n
        n.active_prev = n
        n.active_next = n
        n
    end

end
type Presolve_Col <: Presolve_Element
    # Each column represents a variable xj
    j :: Int                                # Col number in the original problem.
    l :: Float64                           # lb[j] in the original problem
    u :: Float64                           # ub[j] in the original problem
    c_val :: Float64                        # c[j] in the original problem.
    aij :: Union{Presolve_Element,Nothing}  # Reference to the first non-zero matrix entry in col j. aij = nothing indicates empty col.
    prev :: Presolve_Col                    # Reference to previous Col in the original problem. Col 5 holds ref to Col 4 etc
    next :: Presolve_Col                    # Reference to next Col in the original problem. Col 5 holds ref to Col 6 etc
    is_independent :: Bool                  # true if the Col is in the independent doubly linked list. false otherwise.
    ind_prev :: Presolve_Col                # if in the independent doubly linked list, reference to the previous Col in the independent doubly linked list.
    ind_next :: Presolve_Col                # if in the independent doubly linked list, reference to the next Col in the independent doubly linked list.
    function Presolve_Col()
        # Constructor that creates a default row which has an invalid row.i value.
        n = new()
        n.j = -1
        n.l = -Inf
        n.u = +Inf
        n.c_val = 0.0
        n.is_independent = true
        n.aij = nothing
        n.prev = n
        n.next = n
        n.ind_prev = n
        n.ind_next = n
        n
    end
end

type Presolve_Matrix <: Presolve_Element
    row :: Presolve_Row                     # the row associated with element aij. for example, a(4,5) will hold row 4.
    col :: Presolve_Col                     # the col associated with element aij. for example, a(4,5) will hold col 5.
    val :: Float64                          # the matrix value in the constraint matrix
    row_prev :: Presolve_Matrix             # reference to the previous aij element along the same row.
    row_next :: Presolve_Matrix             # reference to the next aij element along the same row.
    col_prev :: Presolve_Matrix             # reference to the previous aij element along the same col.
    col_next :: Presolve_Matrix             # reference to the next aij element along the same col.
    function Presolve_Matrix()
        n = new()
        row = Presolve_Row()
        col = Presolve_Col()
        val = 0.0
        row_prev = n
        row_next = n
        col_prev = n
        col_next = n
        n
    end
end

`
--- Presolve Stack ---
There are multiple redundancies possible , we combine those that are similar and can be resolved together
into a subtype of the abstract type Presolve_Stack.
As of now we have implemented Linear_Dependency which can help resolve - singleton rows, singleton columns and forcing constraints

The overall abstract type is called "Presolve_Stack" as we need to do the post-solving in the reverse order
and hence we refer to it as a stack for the LIFO logic.

Each subtype contains only as much information as is required for the postsolving.
`

abstract Presolve_Stack

`
--- Linear Dependency ---
Here we detect that a variable xj is linearly  dependent on some other variables by the equation -
xj = constant + ∑ xvalues * linear co-efficients

Consider a column singleton in column 5 and the corresponding element being a(4,5) (let dimensions be m,n = 6,6)
Suppose the fourth row of the constraint matrix looks like this -
                                b
4th row :   1 4 0 2 3 -1        10

This represents the equation -
1*x_1 + 4*x_2 + 0*x_3 + 2*x_4 + 3*x_5 + (-1)*x_6 = 10.

We can substitute x5 out of the problem as x5  = 10/3 - (1/3)*x_1 - (4/3)*x2 - (0/3)*x3 - (2/3)*x4 - (-1/3)*x6

We view this as
x[index] = value + ∑ x[vec1[i]]*vec2[i]
where,
index here is 5 for x_5
value here is 10/3
vec1 here is [1,2,4,6] (not 3 as coefficient of 3 is 0)
vec2 here is [-1/3, -4/3, -2/3, 1/3]

For a row singleton vec1 and vec2 are empty.
`

type Linear_Dependency <: Presolve_Stack
    index :: Int                            # index of the variable that is being presolved wrt original problem.
    vec1 :: Vector{Int}                     # indices of the variables that x[index] is dependent on
    vec2 :: Vector{Float64}                 # corresponding co-efficients of dependency. See above for an example
    value :: Float64                        # constant value of the Linear Dependence equation

    function Linear_Dependency(ind::Int, val::Number)
        vec1 = Array{Int,1}()
        vec2 = Array{Float64,1}()
        new(ind,vec1,vec2,val)
    end

    function Linear_Dependency(ind::Int, vec1::Vector{Int}, vec2::Vector{Float64}, val::Number)
        new(ind,vec1,vec2,val)
    end

end

`
--- Presolve Problem ---
We take the original problem and work with internally before reporting back the smaller resultant problem.
our internal workspace consist of problem type Presolve_Problem which holds information in the way we want for our internal functions.
This is never accessed by the user and its scope is the presolver! function call.
We have a constructor which initializes the variables to default values.
`

type Presolve_Problem
    # The dimensions of the original problem
    originalm :: Int64                      # number of rows in original problem
    originaln :: Int64                      # number of cols in original problem

    # Linked List storage
    dictrow :: Dict{Int64,Presolve_Row}     # Rows of the original problem.
    dictcol :: Dict{Int64,Presolve_Col}     # Cols of the original problem.
    dictaij :: Dict{Int64,Presolve_Matrix}  # Non-zero entries of the constraint matrix
    rowptr :: Presolve_Row                  # Reference to the first valid row. Will be updated as rows are deleted.
    colptr :: Presolve_Col                  # Reference to the first valid col. Will be updated as cols are deleted.
    rowque :: Presolve_Row                  # Reference to the first active row. Will be updated as rows are dequed.
    colque :: Presolve_Col                  # Reference to the first independent col. Will be updated as cols are dequed.

    # Boolean status fields
    independentvar :: BitArray{1}
    activeconstr :: BitArray{1}

    # counter variables for aij elements in row/col. Can be done away with probably.
    rowcounter :: Array{Float64,1}
    colcounter :: Array{Float64,1}

    # the stack that will be fed into the postsolver.
    pstack :: Array{Presolve_Stack,1}

    # map from the index of the inital rows to final rows. -1 if they are deleted.
    finalrows :: Array{Int,1}
    finalcols :: Array{Int,1}

    # Constructor that creates a Presolve_Problem
    function Presolve_Problem(verbose::Bool,m::Int,n::Int)
        verbose && println("-----------INSIDE PRESOLVE CONSTRUCTOR------------")

        originalm,originaln = m,n
        dictrow = Dict{Int64,Presolve_Row}()
        dictcol = Dict{Int64,Presolve_Col}()
        dictaij = Dict{Int64,Presolve_Matrix}()
        rowptr = Presolve_Row()
        colptr = Presolve_Col()
        rowque = Presolve_Row()
        colque = Presolve_Col()
        independentvar = trues(originaln)
        activeconstr = trues(originalm)
        rowcounter = zeros(originalm)
        colcounter = zeros(originaln)
        pstack = Array{Presolve_Stack,1}()
        finalrows = fill(-1,originalm)              # Initially everything is -1.
        finalcols = fill(-1,originaln)

        new(originalm,originaln,dictrow,dictcol,dictaij,rowptr,colptr,rowque,colque,independentvar,activeconstr,rowcounter,colcounter,pstack,finalrows,finalcols)
    end
    function Presolve_Problem()
        return Presolve_Problem(false,0,0)
    end
end
