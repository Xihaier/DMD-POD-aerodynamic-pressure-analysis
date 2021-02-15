%%
% assemble Hankel matrix

function [output_data] = Hankel_matrix(data, num)

output_data = zeros(125*num.delay, num.snapshots - num.delay + 1);

for kk = 1:num.delay
    output_data((125*(kk-1)+1):(125*kk), :) = data(:, kk:(num.snapshots-num.delay+kk));
end

end


