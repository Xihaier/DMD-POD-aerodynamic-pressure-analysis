%%
% summarize DMD results


function results_summary(DMD_infor, X1, dt)

resolution = [125 575];

num_x = resolution(1);
num_y = resolution(2);
xx = -0.04:0.02:0.04; xx = repmat(xx, 25, 1);
yy = 0.01:0.02:0.49; yy = flip(repmat(yy, 5, 1)');
Xq = linspace(-0.04,0.04,num_x); Xq = repmat(Xq, num_y, 1);
Yq = linspace(0.01,0.49,num_y); Yq = flip(repmat(Yq, num_x, 1)');

fre = log(DMD_infor.val)/dt;
time = (0:(size(X1, 2)-1))*dt;
[fre, time] = meshgrid(fre, time);

vander = exp(fre'.*time');
dynamics = repmat(DMD_infor.e_amp, 1, size(X1, 2)) .* vander;

[~, index] = sort(abs(imag(log(DMD_infor.val)/dt)));
dynamics = dynamics(index, :);

real_e = real(dynamics);
for kk = 1:125
    real_e(kk,:) = real_e(kk,:)/max(abs(real_e(kk,:)));
end

imag_e = imag(dynamics);
for kk = 1:125
    imag_e(kk,:) = imag_e(kk,:)/max(abs(imag_e(kk,:)));
end

eig_vec = real(DMD_infor.e_modes);

aa = [15 57 117];

for jj=1:3
    
    kk = aa(jj);
    
    figure;
    data_integral = interp2(xx, yy, reshape(eig_vec(1:125,kk), 25, 5), Xq, Yq, 'linear');
    imagesc(data_integral);
    set(gcf,'Position',[1 1 round(290) round(870)]);
    colormap(curl);
    colorbar;
    axis equal;
    axis off;
    set(gca, 'FontSize',25, 'FontName','Times', 'LineWidth',1.2);
    set(gcf,'Position',[1 1 round(330) round(990)]);
    
    figure;
    plot(1:1000, real_e(kk,:), 'LineWidth',3);
    xlabel('time step');
    grid on;
    set(gcf,'Position',[1 1 round(870) round(290)]);
    set(gca, 'FontSize',25, 'FontName','Times', 'LineWidth',1.2);
    
    figure;
    [pxx1,w1] = periodogram(real_e(kk,:));
    plot(w1/pi,10*log10(pxx1), 'r', 'LineWidth',3);
    xlabel('normalized frequency');
    grid on;
    set(gcf,'Position',[1 1 round(870) round(290)]);
    set(gca, 'FontSize',25, 'FontName','Times', 'LineWidth',1.2);
    
end

end


