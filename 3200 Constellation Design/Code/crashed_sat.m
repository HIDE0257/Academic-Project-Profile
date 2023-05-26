%% Crashed Satellite Check
function crashed = crashed_sat(Nsc,r_ACI)
% Inputs:   Nsc = Number of Satellites 
%           r_ACI = Satellite position in the ACI frame
% Outputs:  crashed = cell array {#ofCombination x 1}(#ofIndex x 3)[ith sat, jth sat, 1/0] 

% Ncomb = factorial(Nsc)/(factorial(2)*factorial(Nsc-2));
% crashed = zeros(length(r_ACI{1,1}),3);
comb = 2;
col_comb = 0;
for i = 1:Nsc-1
    for j = comb:Nsc
        col_comb = col_comb + 1;
        for k = 1:length(r_ACI{1,1})
            if r_ACI{i}(k,:) == r_ACI{j}(k,:)
                crashed{col_comb,1}(k,:) = {i, j, 1};
                fprintf('Our satellites (%d & %d) are crashed, baby!\n', i,j);
            else
                crashed{col_comb,1}(k,:) = {i, j, 0};
            end
        end
    end
    comb = comb + 1; 
end

end

