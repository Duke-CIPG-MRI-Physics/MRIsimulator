function [Sampling_Table_PF] = PF_Undersample(Complete_Matrix_Size,Sampling_Table,PF_Factor)
%This function removes entries from a sampling table, in accordance with
%Partial Fourier undersampling

arguments (Input)
    %Complete_Matrix_Size is the final image size
    Complete_Matrix_Size (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}

    %Sampling_Table is the output from TWIST.m, it can have already been
    %passed through GRAPPA_Undersample.m
    Sampling_Table (:,4) 

    %PF_Factor is the percentage of k-space which will be acquired
    PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor,1), mustBeGreaterThan(PF_Factor,.5)}

end

% %% --- Deleting PF areas
% Partial_Fourier_Number = [round(Complete_Matrix_Size(2)*PF_Factor(1)), round(Complete_Matrix_Size(3)*PF_Factor(2))];
% 
% Rows_to_keep = 1:Partial_Fourier_Number(1);
% Columns_to_keep = (Complete_Matrix_Size(3)-Partial_Fourier_Number(2)):Complete_Matrix_Size(3);
% 
% Sampling_Table_Accelerated = [];
% 
% for ii = 1:height(Sampling_Table)
%     if any(Sampling_Table{ii,"Row (phase)"} == Rows_to_keep) && any(Sampling_Table{ii,"Column (slice)"} == Columns_to_keep)
% 
%         Sampling_Table_Accelerated = [Sampling_Table_Accelerated;Sampling_Table(ii,:)];
% 
%     end
% 
% end
% 
% Sampling_Table_PF = Sampling_Table_Accelerated;

%% --- Deleting PF areas (Vectorized) ---

% 1. Define the rows and columns to keep (no change needed)
Partial_Fourier_Number = [round(Complete_Matrix_Size(2) * PF_Factor(1)), round(Complete_Matrix_Size(3) * PF_Factor(2))];
Rows_to_keep = 1:Partial_Fourier_Number(1);
Columns_to_keep = (Complete_Matrix_Size(3) - Partial_Fourier_Number(2)):Complete_Matrix_Size(3);

% 2. Use ismember() to find all rows that match the criteria in a single step
is_row_kept = ismember(Sampling_Table.("Row (phase)"), Rows_to_keep);
is_col_kept = ismember(Sampling_Table.("Column (slice)"), Columns_to_keep);

% 3. Combine the conditions using an element-wise AND (&)
final_mask = is_row_kept & is_col_kept;

% 4. Use the final logical mask to select all desired rows at once
Sampling_Table_PF = Sampling_Table(final_mask, :);
end