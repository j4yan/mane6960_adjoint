#
# Create Jacobian matrix
#
function createJac{Tsol}(q::AbstractArray{Tsol, 3})
    n_field = size(q, 1)
    n_dof_per_elem = size(q, 2)
    n_elem = size(q, 3)

    m = length(q)
    n = m
    blk_size = n_field * n_dof_per_elem
    nnz = blk_size * blk_size * (n_elem * 3 - 2)
    colptr = zeros(Int32, m+1)
    rowval = zeros(Int32, nnz)
    nzval = zeros(Float64, nnz)

    # special treatment for 1st and last element
    colptr[1] = 1
    colptr[2 : blk_size+1] = blk_size * 2
    colptr[m-blk_size+1 : m+1] = blk_size * 2
    for el = 2 : n_elem - 1
        for j = 1 : blk_size
            col = (el - 1) * blk_size + j
            colptr[col+1] = blk_size * 3
        end
    end

    for j = 1 : m
        colptr[j+1] += colptr[j]
    end

    @assert(colptr[m+1] == nnz + 1)

    # special treatment for 1st and last element
    for j = 1 : blk_size    # loop over local columns
        # first element
        row_offset = 0
        col_offset = 0
        col = col_offset + j
        for i = 1 : 2*blk_size # loop over rows
            loc = colptr[col] + i - 1
            rowval[loc] = i + row_offset
        end
        # last element
        row_offset = m - 2*blk_size
        col_offset = m - blk_size
        col = col_offset + j
        for i = 1 : 2*blk_size    # loop over rows 
            loc = colptr[col] + i - 1
            rowval[loc] = i + row_offset
        end
    end

    for el = 2 : n_elem - 1
        col0 = (el - 1) * blk_size + 1 
        col1 = col0 + blk_size - 1
        row0 = col0 - blk_size
        row1 = row0 + 3*blk_size - 1
        for col = col0 : col1
            for i = 1 : blk_size * 3
                rowval[colptr[col] + i - 1] = row0 + i -1
            end
        end
    end
    
    return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end


