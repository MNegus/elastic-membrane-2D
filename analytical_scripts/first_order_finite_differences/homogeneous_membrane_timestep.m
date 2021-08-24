function w_next = homogeneous_membrane_timestep(w, w_previous, A_mat)

        
        %% Determine right-hand-side and solve for w_next
        w_next = A_mat \ (2 * w - w_previous);
end