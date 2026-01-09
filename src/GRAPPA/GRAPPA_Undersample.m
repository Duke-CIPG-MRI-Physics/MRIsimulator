function [Sampling_Table_GRAPPA] = GRAPPA_Undersample(Complete_Matrix_Size,Sampling_Table,R)
%This function removes entries from a sampling table, in accordance with
%GRAPPA undersampling
%  
%The strategy is to zero lines according to the R factor

arguments
    %Complete_Matrix_Size is the final image size
    Complete_Matrix_Size (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}

    %Sampling_Table is the output from TWIST.m, it can have already been
    %passed through PF_Undersample.m
    Sampling_Table (:,4) 

    %GRAPPA acceleration factors: [phase (rows), slice (columns)]
    R (1,2) {mustBeNumeric,mustBeInteger,mustBePositive}
   
end


%% --- Deleting GRAPPA Lines ---

% 1. Define the regularly-spaced rows and columns to keep (no change needed)
Rows_to_keep = R(1):R(1):Complete_Matrix_Size(2);
Columns_to_keep = R(2):R(2):Complete_Matrix_Size(3);

% 2. Use ismember() to find all rows that match the criteria in a single operation.
% This creates a logical vector (true/false) for each condition.
is_row_kept = ismember(Sampling_Table.("Row (phase)"), Rows_to_keep);
is_col_kept = ismember(Sampling_Table.("Column (slice)"), Columns_to_keep);

% 3. Combine the two logical vectors using an element-wise AND (&).
% The final mask is true only for rows that satisfy BOTH conditions.
final_mask = is_row_kept & is_col_kept;

% 4. Use the final logical mask to select all desired rows at once.
% This is extremely fast and avoids loops and resizing arrays.
Sampling_Table_GRAPPA = Sampling_Table(final_mask, :);
end


