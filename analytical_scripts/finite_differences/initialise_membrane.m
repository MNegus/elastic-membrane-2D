function w = initialise_membrane(w_previous, A_mat, B_mat)
        rhs = 0.5 * B_mat * w_previous;
        w = A_mat \ rhs;
end