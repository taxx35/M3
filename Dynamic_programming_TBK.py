# Importing required packages
import numpy as np
import time  # Importing the time module for measuring execution time

# Defining seq1 and seq2 and user prompt for inputing sequences
seq1 = input("Enter the first sequence:")
seq2 = input("Enter the second sequence:")

s1 = len(seq1)
s2 = len(seq2)
if seq1 and seq2:
    # Time measurement for Global Alignment
    start_time_global = time.time()

    # Defining the rewards and penalties
    match = 1
    mismatch = -1
    gap = -2

    # Lengths of seq1 and seq2
    m = len(seq1)
    n = len(seq2)

    # Initializing the matrix with zeros, size (m+1) x (n+1) -> for global and local matrix
    matrix_global = [[0] * (n + 1) for _ in range(m + 1)]
    matrix_local = [[0] * (n + 1) for _ in range(m + 1)]

    # Filling the first row and column with gap penalties for global alignment
    for c in range(m + 1):
        matrix_global[c][0] = c * gap  # Sets value in the first column of each row to c * gap.
    for r in range(n + 1):
        matrix_global[0][r] = r * gap  # Sets value in the first row of each column to r * gap.

    # Filling the first row and column with zeros for local alignment
    for i in range(m + 1):
        matrix_local[i][0] = 0  # Sets value in the first column of each row to 0
    for j in range(n + 1):
        matrix_local[0][j] = 0  # Sets value in the first row of each column to 0

    # Filling in the rest of the matrix for global alignment
    # Iterating over the columns of the matrix (starting second column - index 1)
    for c in range(1, m + 1):
        # Iterating over the rows of the matrix (starting second row - index 1)
        for r in range(1, n + 1):

            # If characters match, use match score; otherwise, use mismatch score
            score = match if seq1[c - 1] == seq2[r - 1] else mismatch

            # Fill rest of the matrix with max score
            matrix_global[c][r] = max(
                matrix_global[c - 1][r - 1] + score,  # Diagonal (match or mismatch)
                matrix_global[c][r - 1] + gap,         # Horizontal (gap in seq2)
                matrix_global[c - 1][r] + gap          # Vertical (gap in seq1)
            )

    # Traceback process for global alignment
    # Empty list to store the aligned sequences
    aligned1_g = []
    aligned2_g = []
    alignment_score_global = 0  # Initialize the alignment score

    # Starting traceback from the bottom-right corner (matrix[m][n])
    c = m
    r = n

    # Loop until both c and r reach 0 -> top-left corner of matrix
    while c > 0 or r > 0:
        score = match if seq1[c - 1] == seq2[r - 1] else mismatch
        if matrix_global[c][r] == matrix_global[c - 1][r - 1] + score:
                # Match or mismatch (diagonal move)
                aligned1_g.append(seq1[c - 1])  # Append the characters to aligned1_g
                aligned2_g.append(seq2[r - 1])
                alignment_score_global += score
                c -= 1  # Both c and r are decremented to move diagonally up-left.
                r -= 1
        elif c > 0 and matrix_global[c][r] == matrix_global[c - 1][r] + gap:
            # Gap in seq2 (up move)
            aligned1_g.append(seq1[c - 1])  # Append the characters to aligned1_g
            aligned2_g.append('-')  # Append "-" to aligned2_g for that character in aligned1_g
            alignment_score_global += gap
            c -= 1
        elif r > 0 and matrix_global[c][r] == matrix_global[c][r - 1] + gap:
            # Gap in seq1 (left move)
            aligned1_g.append('-')
            aligned2_g.append(seq2[r - 1])
            alignment_score_global += gap
            r -= 1

    # Reverse the alignments (since we built them backwards)
    aligned1_g.reverse()
    aligned2_g.reverse()

    # Convert lists into strings
    aligned_seq1_global = ''.join(aligned1_g)
    aligned_seq2_global = ''.join(aligned2_g)

    # End time measurement for Global Alignment
    end_time_global = time.time()
    elapsed_time_global = end_time_global - start_time_global  # Time taken for global alignment

    # Time measurement for Local Alignment
    start_time_local = time.time()

    # Filling in the rest of the matrix for local alignment (similar to global but with 0 initiation)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # If characters match, use match score; otherwise, use mismatch score
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            # For local alignment, we take the maximum score from diagonal, horizontal, or vertical,
            # but we don't allow negative scores, so we use 0 instead of gap penalties if the value is negative
            matrix_local[i][j] = max(
                0,  # For excluding the negative scores
                matrix_local[i - 1][j - 1] + score,  # Diagonal (match or mismatch)
                matrix_local[i][j - 1] + gap,         # Horizontal (gap in seq2)
                matrix_local[i - 1][j] + gap          # Vertical (gap in seq1)
            )

    # Traceback process for local alignment
    # Empty lists to store the aligned sequences
    aligned1_l = []
    aligned2_l = []
    alignment_score_local = 0  # Initialize the alignment score

    # Start traceback from the highest score in the matrix
    # Find the maximum score and its position
    max_score = 0
    max_pos = (0, 0)

    # Create loop to iterate over the matrix for max score
    for i in range(m + 1):
        for j in range(n + 1):
            if matrix_local[i][j] > max_score:
                max_score = matrix_local[i][j]
                max_pos = (i, j)

    i, j = max_pos

    # Loop until both i and j reach 0 -> top-left corner of matrix
    while i > 0 and j > 0 and matrix_local[i][j] > 0:
        score = match if seq1[i - 1] == seq2[j - 1] else mismatch

        if matrix_local[i][j] == matrix_local[i - 1][j - 1] + score:
            # Match or mismatch (diagonal move)
            aligned1_l.append(seq1[i - 1])
            aligned2_l.append(seq2[j - 1])
            alignment_score_local += score
            i -= 1
            j -= 1
        elif matrix_local[i][j] == matrix_local[i - 1][j] + gap:
            # Gap in seq2 (up move)
            aligned1_l.append(seq1[i - 1])
            aligned2_l.append('-')
            alignment_score_local += gap
            i -= 1
        elif matrix_local[i][j] == matrix_local[i][j - 1] + gap:
            # Gap in seq1 (left move)
            aligned1_l.append('-')
            aligned2_l.append(seq2[j - 1])
            alignment_score_local += gap
            j -= 1

    # Reverse the alignments (since we built them backwards)
    aligned1_l.reverse()
    aligned2_l.reverse()

    # Convert lists into strings
    aligned_seq1_local = ''.join(aligned1_l)
    aligned_seq2_local = ''.join(aligned2_l)

    # End time measurement for Local Alignment
    end_time_local = time.time()
    elapsed_time_local = end_time_local - start_time_local  # Time taken for local alignment

    # Display results of alignment, along with score, time taken and length 
    print("Global Alignment")
    print(aligned_seq1_global)
    print(aligned_seq2_global)
    print("Alignment score:", alignment_score_global)
    print("Time taken for global alignment:", elapsed_time_global, " seconds")  # Display the execution time
    print("The length of the sequence 1 is:", s1, "Bases")  # Number of bases for seq1 and seq2 
    print("The length of the sequence 2 is:", s2, "Bases")

    # Display the aligned sequences and alignment score 
    print("Local Alignment")
    print(aligned_seq1_local)
    print(aligned_seq2_local)
    print("Alignment score:", alignment_score_local)
    print("Time taken for local alignment:", elapsed_time_local, " seconds")
   
