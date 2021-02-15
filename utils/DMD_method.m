%%
% perform dynamic mode decomposition

function [DMD_infor] = DMD_method(X1, X2, num, dt)

% singular value decomposition
[U, S, V] = svd(X1, 'econ');

pre_select = 0;
if pre_select == 1
    sigma = diag(S);
    beta = size(X1,1)/size(X1,2);
    if beta>1
        beta=0.99;
    end
    thresh = optimal_SVHT_coef(beta,0) * median(sigma);
    r = length(sigma(sigma>thresh));
    r = min(num.truncate,r);
    num.truncate = r;
end

U_r = U(:, 1:num.truncate);
S_r = S(1:num.truncate, 1:num.truncate);
V_r = V(:, 1:num.truncate);

% build low rank approximation
Atilde = U_r' * X2 * V_r / S_r;

% compute eigenvalues and eigenmodes
[eig_vecs, eig_vals] = eig(Atilde);
eig_vals = diag(eig_vals);

% DMD modes
dmd_modes_p = U_r*eig_vecs;
dmd_modes_e = X2*V_r/S_r*eig_vecs;

% amplitudes
amplitudes_p = pinv(dmd_modes_p)*X1(:,1);
amplitudes_e = pinv(dmd_modes_e)*X1(:,1);

% frequencies
frequencies = imag(log(eig_vals))/(2*pi*dt);

% growth rates
growth_rates = real(log(eig_vals))/dt;

% output data
DMD_infor.val = eig_vals;
DMD_infor.vec = eig_vecs;

DMD_infor.p_modes = dmd_modes_p;
DMD_infor.e_modes = dmd_modes_e;

DMD_infor.p_amp = amplitudes_p;
DMD_infor.e_amp = amplitudes_e;

DMD_infor.fre = frequencies;

DMD_infor.growth = growth_rates;

end


