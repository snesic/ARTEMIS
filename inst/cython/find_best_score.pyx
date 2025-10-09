# distutils: extra_compile_args=-O3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
import numpy as np
cimport numpy as np
import math

# Define C types for better performance
ctypedef np.float64_t DTYPE_t
ctypedef np.int32_t ITYPE_t

def find_best_score(
    np.ndarray[DTYPE_t, ndim=2] H,
    ITYPE_t s1_len,
    ITYPE_t s2_len,
    int mem,  # mem is often used as a control flow variable, keeping as 'int' is fine
    int verbose
):
    """
    Cythonized version of find_best_score.
    Calculates best score and a list of potential alignment endpoints.
    """
    cdef np.ndarray[DTYPE_t, ndim=1] mem_score
    cdef np.ndarray[ITYPE_t, ndim=2] mem_index
    cdef DTYPE_t max_score
    cdef ITYPE_t max_idx, i
    cdef tuple max_index
    cdef np.ndarray mem_array  # Use generic array for the zip/object array
    cdef DTYPE_t mem_min
    cdef np.ndarray mem_score_final
    cdef np.ndarray mem_index_final
    cdef DTYPE_t finalScore

    # 1. Get score and index arrays
    mem_score = H[:, s1_len].copy()
    
    mem_index = np.empty((s2_len + 1, 2), dtype=np.int32)
    mem_index[:, 0] = np.arange(s2_len + 1, dtype=np.int32)
    mem_index[:, 1] = s1_len

    # 2. Find max score
    max_idx = np.argmax(mem_score)
    max_score = mem_score[max_idx]
    finalScore = max_score
    max_index = (max_idx, s1_len)
    finalIndex = max_index


    # 3. Sort arrays based on score
    sort_indices = mem_score.argsort()
    # Apply the sort indices to both arrays
    mem_score_sorted = mem_score[sort_indices]
    mem_index_sorted = mem_index[sort_indices]
    # 4. Filtering Logic

    if mem == -1:
        mem = max(1, math.floor(s2_len / s1_len))
        
        # Get the score of the N-th best alignment (using negative index on sorted array)
        mem_min = mem_score_sorted[-mem]
        
        # Create a boolean mask for filtering
        mask = mem_score_sorted >= mem_min * 0.9

        # Apply the mask to both arrays (these slices remain 1D/2D typed arrays)
        mem_score_final = mem_score_sorted[mask]
        mem_index_final = mem_index_sorted[mask]
        
        if verbose == 2:
            print("Calculated mem: ")
            print(mem)

    elif mem == 0:
        # Return empty 2D int array and 1D float array
        mem_index_final = np.empty((0, 2), dtype=np.int32)
        mem_score_final = np.empty((0,), dtype=np.float64)

    else: # mem >= 1
        mem_min = mem_score_sorted[-mem]
        
        mask = mem_score_sorted >= mem_min * 0.9
        
        mem_score_final = mem_score_sorted[mask]
        mem_index_final = mem_index_sorted[mask]


    # 5. Reverse the order (highest score first)
    # The arrays are currently sorted ascending by score; reverse to get descending.
    mem_score_final = mem_score_final[::-1]
    mem_index_final = mem_index_final[::-1]

    return finalScore, finalIndex, mem_index_final, mem_score_final
