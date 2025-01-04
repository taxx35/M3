# pip install streamlit -> install the required streamlit package in terminal
# Importing required packages
import numpy as np
import streamlit as st
import time  # Importing the time module for measuring execution time

# Streamlit UI
# Title
st.title("Sequence Alignment Tools")

# Side Header -> Global, Local alignment
st.sidebar.header("Alignment Options")
align_type = st.sidebar.radio("Choose Alignment Type", options=["Global", "Local"])

# Header
st.header("Input Sequences")
seq1 = st.text_input("Enter Sequence 1:", placeholder="e.g., GATTACA")
seq2 = st.text_input("Enter Sequence 2:", placeholder="e.g., ATTACATTAC")

# If align button pressed, then perform the following commands
if st.button("Align"):
    if seq1 and seq2:
        start_time = time.time()

        # Defining the rewards and penalties
        match = 1
        mismatch = -1
        gap = -1

        # Lengths of seq1 and seq2
        m = len(seq1)
        n = len(seq2)

        # Initializing the matrix with zeros, size (m+1) x (n+1) -> for global and local
        matrix_global = [[0] * (n + 1) for _ in range(m + 1)]
        matrix_local = [[0] * (n + 1) for _ in range(m + 1)]
        
       # Fill the first row and column with gap penalties for global alignment 
        for c in range(m + 1):
            matrix_global[c][0] = c * gap
        for r in range(n + 1):
            matrix_global[0][r] = r * gap

        # Filling the first row and column with zeros for local alignment 
        for i in range(m + 1):
            matrix_local[i][0] = 0  # First column of matrix for local alignment

        for j in range(n + 1):
            matrix_local[0][j] = 0  # First row of matrix for local alignment

        if align_type == "Global":
            # Filling in the rest of the matrix for global alignment
            for c in range(1, m + 1):
                for r in range(1, n + 1):
                    # If characters match, use match score; otherwise, use mismatch score
                    if seq1[c - 1] == seq2[r - 1]:
                        score = match
                    else:
                        score = mismatch

                    # Fill rest of the matrix with max score
                    matrix_global[c][r] = max(
                        matrix_global[c - 1][r - 1] + score,  # Diagonal (match or mismatch)
                        matrix_global[c][r - 1] + gap,         # Horizontal (gap in seq2)
                        matrix_global[c - 1][r] + gap          # Vertical (gap in seq1)
                    )

            # Traceback process for global alignment
            aligned1 = []
            aligned2 = []
            match_line = []  # For the match/mismatch line

            # Starting traceback from the bottom-right corner (matrix[m][n])
            c = m
            r = n
            alignment_score = 0  # Initialize the alignment score

            while c > 0 or r > 0:
                if c > 0 and r > 0 and matrix_global[c][r] == matrix_global[c - 1][r - 1] + (match if seq1[c - 1] == seq2[r - 1] else mismatch):
                    # Match or mismatch (diagonal move)
                    aligned1.append(seq1[c - 1])
                    aligned2.append(seq2[r - 1])
                    alignment_score += match if seq1[c - 1] == seq2[r - 1] else mismatch
                    c -= 1
                    r -= 1
                elif c > 0 and matrix_global[c][r] == matrix_global[c - 1][r] + gap:
                    # Gap in seq2 (up move)
                    aligned1.append(seq1[c - 1])
                    aligned2.append('-')
                    alignment_score += gap
                    c -= 1
                elif r > 0 and matrix_global[c][r] == matrix_global[c][r - 1] + gap:
                    # Gap in seq1 (left move)
                    aligned1.append('-')
                    aligned2.append(seq2[r - 1])
                    alignment_score += gap
                    r -= 1

            # Reverse the alignments (since we built them backwards)
            aligned1.reverse()
            aligned2.reverse()

            # Display the aligned sequences and alignment score
            st.subheader("Global Alignment Result")

            # Convert lists into strings
            aligned_seq1 = ''.join(aligned1)
            aligned_seq2 = ''.join(aligned2)

            # End time measurement
            end_time = time.time()
            elapsed_time = end_time - start_time  # Time taken for the alignment

            # Display results of alignment, along with score and time taken
            st.text(f"{aligned_seq1}")
            st.text(f"{aligned_seq2}")
            st.text(f"Alignment score: {alignment_score}")
            st.text(f"Time taken for alignment: {elapsed_time:.4f} seconds")  # Display the execution time
            st.text(f"The length of the sequence 1 is:   {m}  Bases") # Number of bases for seq1 and seq2 
            st.text(f"The length of the sequence 2 is:   {n}  Bases")

        elif align_type == "Local":
            # Filling in the rest of the matrix for local alignment (similar to global but with 0 initiation)
            for i in range(1, m + 1):
                for j in range(1, n + 1):
                    # If characters match, use match score; otherwise, use mismatch score
                    if seq1[i - 1] == seq2[j - 1]:
                        score = match
                    else:
                        score = mismatch

                    # For local alignment, we take the maximum score from diagonal, horizontal, or vertical,
                    # but we don't allow negative scores, so we use 0 instead of gap penalties if the value is negative
                    matrix_local[i][j] = max(
                        0,  # For excluding the negative scores
                        matrix_local[i - 1][j - 1] + score,  # Diagonal (match or mismatch)
                        matrix_local[i][j - 1] + gap,         # Horizontal (gap in seq2)
                        matrix_local[i - 1][j] + gap          # Vertical (gap in seq1)
                    )

            # Traceback process for local alignment
            aligned1 = []
            aligned2 = []
            alignment_score = 0  # Initialize the alignment score

            # Start traceback from the highest score in the matrix
            # Find the maximum score and its position
            max_score = 0
            max_pos = (0, 0)
            for i in range(m + 1):
                for j in range(n + 1):
                    if matrix_local[i][j] > max_score:
                        max_score = matrix_local[i][j]
                        max_pos = (i, j)

            i, j = max_pos

            while i > 0 and j > 0 and matrix_local[i][j] > 0:
                if matrix_local[i][j] == matrix_local[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
                    # Match or mismatch (diagonal move)
                    aligned1.append(seq1[i - 1])
                    aligned2.append(seq2[j - 1])
                    alignment_score += match if seq1[i - 1] == seq2[j - 1] else mismatch
                    i -= 1
                    j -= 1
                elif matrix_local[i][j] == matrix_local[i - 1][j] + gap:
                    # Gap in seq2 (up move)
                    aligned1.append(seq1[i - 1])
                    aligned2.append('-')
                    alignment_score += gap
                    i -= 1
                elif matrix_local[i][j] == matrix_local[i][j - 1] + gap:
                    # Gap in seq1 (left move)
                    aligned1.append('-')
                    aligned2.append(seq2[j - 1])
                    alignment_score += gap
                    j -= 1

            # Reverse the alignments (since we built them backwards)
            aligned1.reverse()
            aligned2.reverse()

            # Display the aligned sequences and alignment score
            st.subheader("Local Alignment Result")

            # Convert lists into strings
            aligned_seq1 = ''.join(aligned1)
            aligned_seq2 = ''.join(aligned2)

            # End time measurement
            end_time = time.time()
            elapsed_time = end_time - start_time  # Time taken for the alignment

            # Display results of alignment, along with score and time taken
            st.text(f"Score: {alignment_score}")
            st.text(f"Aligned Sequences:")
            st.text(f"{aligned_seq1}")
            st.text(f"{aligned_seq2}")
            st.text(f"Time taken for alignment: {elapsed_time:.4f} seconds")  # Display the execution time
            st.text(f"The length of the sequence 1 is:   {m}  Bases") # Number of bases for seq1 and seq2 
            st.text(f"The length of the sequence 2 is:   {n}  Bases")
    else:
        st.warning("Please enter both sequences.")
