clc
clear
dat = load('data\wind_speed.txt')';
dmax = max(dat);
dmin = min(dat);
v = dmin : 0.01 : dmax;
%% Distribution of original data
dat_sort = sort(dat);
Possi = 1:length(dat_sort);
Possi = Possi / length(dat_sort);
tic
for j = 1:length(v)
    ori_F(j) = length(find(dat <= v(j))) / length(dat);
end
toc
% ori_F_L is the original CDF low density
ori_F_L = ori_F(1:35:length(ori_F));
v_L = v(1:35:length(v));

for i = 2:length(v_L)
    ori_f(i) = (ori_F_L(i) - ori_F_L(i-1)) ...
                / (v_L(i) - v_L(i-1));
end

figure

bar (v_L,ori_f,'FaceColor','none','EdgeColor','b');

%% Matlab Weibull Toolbox
tic
[p] = wblfit(dat);

MW_c = p(1,1);
MW_k = p(1,2);
MW_f = wblpdf(v,MW_c,MW_k);
MW_F = wblcdf(v,MW_c,MW_k);
toc
hold on;
plot (v,MW_f,'g','LineWidth',1);
grid on
% figure
% plot (v,MW_F,'LineWidth',2);
%% Graphic method
tic
logF = log( -log(1 - ori_F));
logV = log(v);

p = polyfit(logV(1:length(logV)-1),logF(1:length(logV)-1),1)

GM_k = p(1)
GM_c = exp(-p(2) / p(1))
GM_f = ((GM_k / GM_c) .* (v / GM_c).^(GM_k - 1)) ...
        .* exp(-(v / GM_c).^GM_k);
toc
    hold on;
plot (v,GM_f,':r','LineWidth',1);

%% Maximum likelihood method

tic
ML_k = 2;
ML_k_T = 0;
while (abs(ML_k - ML_k_T) >0.001)
   ML_k_T = ML_k;
   ML_k = (sum(dat.^ML_k .* log(dat)) / sum(dat.^ML_k) - ...
            sum(log(dat)/length(dat)))^(-1);
end
ML_c = (sum(dat.^ML_k) / length(dat))^(1 / ML_k);

ML_f = ((ML_k / ML_c) .* (v / ML_c).^(ML_k - 1)) ...
        .* exp(-(v / ML_c).^ML_k);
toc
hold on;
plot (v,ML_f,'-.b','LineWidth',1);


%% Moment method
tic
MM_k = (std(dat) / mean(dat))^(-1.086)
MM_c = mean(dat) / gamma(1 + 1/MM_k)

MM_f = ((MM_k / MM_c) .* (v / MM_c).^(MM_k - 1)) ...
        .* exp(-(v / MM_c).^MM_k);
toc
plot (v,MM_f,'--c','LineWidth',1);

grid on;
xlabel('Wind Speed (m/s)');
ylabel('Density');
title('Wind Speed Distribution');
legend('Original Distribution','Weibull Fit By Matlab', ...
        'Weibull Fit By GM','Weibull Fit By ML', ...
        'Weibull Fit By ML','Location','NorthEast')
