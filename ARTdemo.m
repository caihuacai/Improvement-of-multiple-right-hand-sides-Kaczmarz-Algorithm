%ARTdemo (script) Demonstrates the use of, and the results from, the ART methods.
%
% This script illustrates the use of the ART methods kaczmarz, symmetric
% kaczmarz, and randomized kaczmarz.
%
% The script creates a parallel-beam test problem, adds noise, and solves
% the problems with the ART methods.  The exact solution and the results
% from the methods are shown.
%
% See also: nonnegdemo, SIRTdemo, trainingdemo.

% Maria Saxild-Hansen and Per Chr. Hansen, Mar 11, 2011, DTU Compute.
clear;
close all;
fprintf(1,'\nStarting ARTdemo:\n\n');

% Set the parameters for the test problem.
%N = 50;           % The discretization points.
%theta = 0:5:179;  % No. of used angles
N = 60;           % The discretization points.
theta = 0:179;
P = 180;
%p = 75;           % No. of parallel rays.
% eta = 0.05;       % Relative noise level. z

fprintf(1,'Creating a parallel-bema tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),P]);

% Create the test problem.
[A,b_ex,x_ex] = paralleltomo(N,theta,P);

% k = 3600;
% 
% % 计算奇异值的子集
% sigma = svds(A, k);
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
% 
% % 计算条件数
% cond_A = max(sigma) / min(sigma)

% Noise level.
% delta = eta*norm(b_ex);
% 
% % Add noise to the rhs.
% randn('state',0);
% e = randn(size(b_ex));
% e = delta*e/norm(e); b = b_ex+e;
% 2Qa4e

% scaleFactor = 0.6;  
% noise = randn(size(b_ex)) * scaleFactor;
% b = b_ex + noise;
% 
 b=b_ex;


%  noise = 0.008 * randn(size(b_ex)); 
%  b = b_ex + noise;
 % noise = 0.016 * randn(size(b_ex));
 
 
 
%Show the exact solution.
% figure
% imagesc(reshape(x_ex,N,N)), colormap gray,
% axis image off
% c = caxis;


title('Exact phantom')

% No. of iterations.
% % %%--------------------------------------------------------李顺昌2-----------------------%%
% k = 70000;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% % [Xkacz_max_li2,psnr_li2,ssim_li2,p_li2,normalized_variance_kaczmarz_li2,time_li2] = lishunchang2(A,b,k,x_ex);
% [Xkacz_max_li2,psnr_li2,ssim_li2,p_li2,normalized_variance_kaczmarz_li2,time_li2] = lishunchang2(A,b,k,x_ex);
% figure
% original_image = reshape(Xkacz_max_aver,N,N);
% normalized_image = double(original_image)/255;
% imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_prop,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('kaczmarz_prop reconstruction')

% % %%--------------------------------------------------------李顺昌1-----------------------%%
% k = 70000;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% [Xkacz_max_li1,psnr_li1,ssim_li1,p_li1,normalized_variance_kaczmarz_li1,time_li1] = lishunchang(A,b,k,x_ex);

% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_prop,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('kaczmarz_prop reconstruction')
% %%--------------------------------------------------------changqiaodanyouduan-----------------------%%
k = 5000;
fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
fprintf(1,'\nThis takes a moment ...');
%Perform the kaczmarz iterations.
[Xkacz_max_chang,psnr_chang,ssim_chang,p_chang,normalized_variance_kaczmarz_chang,time_chang] = changqiaodanyouduan(A,b,k,x_ex);
earfigure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_prop,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('kaczmarz_prop reconstruction')
% %%--------------------------------------------------------kaczmarz-----------------------%%
% k = 7e7;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% [Xkacz_max_FGBK,psnr,ssim,p,normalized_variance_kaczmarz,time] = kaczmarz_block(A,b,k,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_FGBK,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('kaczmarz reconstruction')


% %%--------------------------------------------------------kaczmarz_prop-----------------------%%
% k = 70000;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% [Xkacz_max_prop,psnr_prop,ssim_prop,p_prop,normalized_variance_kaczmarz_prop,time_prop] = kaczmarz_prop(A,b,k,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_prop,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('kaczmarz_prop reconstruction')
% %%--------------------------------------------------------cancha_tuiguang-----------------------%%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 700000;
% [Xkacz_max_cancha,psnr_cancha,ssim_cancha,p_cancha,normalized_variance_kaczmarz_max_cancha,time_cancha] = cancha(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_cancha,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('cancha reconstruction')
% k = toc
% % %%--------------------------------------------------------GBK_tuiguang-----------------------%%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 70000;
% [Xkacz_max_GBK,psnr_GBK,ssim_GBK,p_GBK,normalized_variance_kaczmarz_max_GBK,time_GBK] = GBK(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_GBK,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('GBK reconstruction')
% k = toc

% %%--------------------------------------------------------FGBK_tuiguang-----------------------%%
%FGBK
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_FGBK_tuiguang''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 700000;
% [Xkacz_max_FGBK,psnr_FGBK,ssim_FGBK,p_FGBK,normalized_variance_kaczmarz_max_FGBK,time_FGBK] = FGBK_tuiguang(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_FGBK,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('FGBK reconstruction')
% k = toc
% % %%--------------------------------------------------------kaczmarz_max_ave_nostep-----------------------%%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_ave_nostep''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 10000;
% [Xkacz_max_aver_nostep,normalized_variance_kaczmarz_max_aver_nostep,time_nostep] = kaczmarz_max_ave_nostep(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_aver_nostep,N,N)), colormap gray,
% 
% axis image off
% %caxis(c);
% title('FFGBK reconstruction')
% k = toc

% % %%--------------------------------------------------------kaczmarz_max_ave_newweight-----------------------%%
% tic

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_ave_newweight''s method.',k,'LineWidth', 2);
fprintf(1,'\nThis takes a moment ...');
%Perform the kaczmarz iterations.
k_max_aver = 10000;

[Xkacz_max_FFGBK,psnr_FFGBK,ssim_FFGBK,p_FFGBK,normalized_variance_kaczmarz_max_FFGBK,time_FFGBK] = kaczmarz_max_ave_newweigh(A,b,k_max_aver,x_ex);
figure
% original_image = reshape(Xkacz_max_aver,N,N);
% normalized_image = double(original_image)/255;
% imshow(normalized_image,[0,1]);
imagesc(reshape(Xkacz_max_FFGBK,N,N)), colormap gray,
axis image off
%caxis(c);
title('FFGBK reconstruction');
% k = toc
% %%--------------------------------------------------------new_set-----------------------%%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_new_set''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 80000;
% [Xkacz_max_K_means_FGBK,psnr_K_means_FGBK,ssim_K_means_FGBK,p_K_means_FGBK,normalized_variance_kaczmarz_max_K_means_FGBK,time_K_means_FGBK] = kaczmarz_max_ave_new_set_2(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_K_means_FGBK,N,N)), colormap gray,
% 
% axis image off
% %caxis(c);
% title('K-means FGBK reconstruction')
% k = toc
% %%--------------------------------------------------------combine_two_1-----------------------%%
tic
fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_combine_two_1''s method.',k,'LineWidth', 2);
fprintf(1,'\nThis takes a moment ...');
%Perform the kaczmarz iterations.
k_max_aver = 20000;
[Xkacz_K_means_FFGBK,psnr_K_means_FFGBK,ssim_K_means_FFGBK,p_K_means_FFGBK,normalized_variance_kaczmarz_K_means_FFGBK,time_K_means_FFGBK] = combine_two_2(A,b,k_max_aver,x_ex);
figure
imagesc(reshape(Xkacz_K_means_FFGBK,N,N)), colormap gray,
axis image off
%caxis(c);
title('K-means FFGBK reconstruction')
k = toc






%original_image = reshape(Xkacz_max_aver,N,N);
% normalized_image = double(original_image)/255;
% imshow(normalized_image,[0,1]);

% a1 = 1:100:length(time_FGBK);
%  plot(a1,time_FGBK);
%  hold on；
% a2 = 1:100:length(time_nostep);
%  plot(a2,time_nostep);
% hold on；
% a3 = 1:100:length(time_newweight);
%  plot(a3,time_newweight);
% hold on；
% a4 = 1:100:length(time_new_set);
%  plot(a4,time_new_set);
%  hold on；
% a5 = 1:100:length(time_combine);
%  plot(a5,time_combine);
% hold on；
%  plot(1:length(time_FGBK),time_FGBK);
%  hold on;
%   plot(1:length(time_nostep),time_nostep);
%    hold on;
%   plot(1:length(time_newweight),time_newweight);
%    hold on;
%   plot(1:length(time_new_set),time_new_set);
%    hold on;
%   plot(1:length(time_combine),time_combine);
% figure;
% colors = get(gca,'ColorOrder')
%  semilogy(1:length(time_FGBK),time_FGBK,'Color',colors(3,:));
%  hold on;
% %   semilogy(1:length(time_nostep),time_nostep,'Color',colors(5,:));
% % % % %    hold on;
%   semilogy(1:length(time_newweight),time_newweight','Color',colors(6,:));
%    hold on;
%   semilogy(1:length(time_new_set),time_new_set,'Color',colors(7,:));
%    hold on;
%   semilogy(1:length(time_combine),time_combine,'Color',colors(1,:));
%  
%   legend('FGBK','FFGBK','NWFGBK','K-means FGBK','K-means FFGBK');
 


% semilogy(1:length(time_FGBK),time_FGBK,'Color', [0.6 0.8 1]);
%  hold on;
%   semilogy(1:length(time_nostep),time_nostep,'Color', [0.8 0 0]);
%    hold on;
%   semilogy(1:length(time_newweight),time_newweight,'Color', [0.7 1 0.7]');
%    hold on;
%   semilogy(1:length(time_new_set),time_new_set,'Color', [1 0.8 0]);
%    hold on;
%   semilogy(1:length(time_combine),time_combine,'Color', [0.7 0.7 1]);
%   
% legend('FGBK','FFGBK','NWFGBK','K-means FGBK','K-means FFGBK');

% 
% plot(1:length(time), time, 'b-', 'DisplayName', '时间');
% hold on;
% plot(1:length(normalized_variance_kaczmarz_combine_two_1),normalized_variance_kaczmarz_combine_two_1,'DisplayName','精度');
% % hold on；
% legend on('时间'，’精度‘)
% time_sum = cumsum(time)
% plot(time_sum, normalized_variance_kaczmarz_combine_two_1, '-o');
% hold on;
% 添加标题和标签
% title('时间与精度');
% xlabel('时间');
% ylabel('精度');




%%-------------------------------------绘图------------------------------------------------------%%

% a1=1:100:length(normalized_variance_kaczmarz_max_FGBK);a1(1)=1;
% plot(a1,normalized_variance_kaczmarz_max_FGBK(a1),'p--m','MarkerSize',5,'MarkerFaceColor','m','MarkerEdgeColor','m','LineWidth',2)
% hold on;
% a2=1:100:length(normalized_variance_kaczmarz_max_aver_nostep);a2(1)=1;
% plot(a2,normalized_variance_kaczmarz_max_aver_nostep(a2),'d--g','MarkerSize',5,'MarkerFaceColor','g','MarkerEdgeColor','g','LineWidth',2)
% hold on;
% a3=1:100:length(normalized_variance_kaczmarz_max_aver_newweight);a3(1)=1;
% plot(a3,normalized_variance_kaczmarz_max_aver_newweight(a3),'s--r','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2)
% hold on;
% legend('FontSize', 10);
% legend('FGBK','aver_nostep','aver_newweight');




% plot(1:length(normalized_variance_kaczmarz_prop), normalized_variance_kaczmarz_prop, 'b-','LineWidth', 2);
% plot(a1,normalized_variance_kaczmarz_prop(a1),'p--m','MarkerSize',5,'MarkerFaceColor','m','MarkerEdgeColor','m','LineWidth',2)
% hold on;
% plot(a1,normalized_variance_kaczmarz_prop_non(a1),'d--g','MarkerSize',5,'MarkerFaceColor','g','MarkerEdgeColor','g','LineWidth',2)
% hold on;
% plot(a1,normalized_variance_kaczmarz(a1),'s--r','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2)
% hold on;
% plot(a1,normalized_variance_kaczmarz_max_3(a1),'*--b','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2)
% hold on;
% plot(a1,normalized_variance_kaczmarz_max_4(a1),'h--k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
% hold on;
%%-------------------------------------------------kaczmarz_max_aver-----------------------------------------------%%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_aver''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max_aver = 40000;
% [Xkacz_max_aver,normalized_variance_kaczmarz_max_aver] = kaczmarz_max_aver(A,b,k_max_aver,x_ex);
% figure
% % original_image = reshape(Xkacz_max_aver,N,N);
% % normalized_image = double(original_image)/255;
% % imshow(normalized_image,[0,1]);
% imagesc(reshape(Xkacz_max_aver,N,N)), colormap gray,
% 
% axis image off
% %caxis(c);
% title('Kaczmarz max average reconstruction')
% k = toc
% %%-----------------------------------------------cancha两种概率选择方式都是残差----------------------------%%
% tic
% k_prop_cancha = 40000;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_prop_non''s method.',k);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% %Xkacz = kaczmarz(A,b,k);
% [Xkacz_prop_cancha,normalized_variance_kaczmarz_constrain_cancha] = untitled(A,b,k_prop_cancha,x_ex);
% %Show the kaczmarz solution.
% figure
% imagesc(reshape(Xkacz_prop_cancha,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('Kaczmarz cancha reconstruction')
% hold on
% 
% k = toc;
% k
% plot(1:length(normalized_variance_kaczmarz_constrain), normalized_variance_kaczmarz_constrain, 'b-', 'DisplayName', '正则化');
% hold on;
%---------------------------------------------constrain_kaczmarz_nonzero----------------------------------------------------%
% tic
% k_prop = 500;
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_prop_non''s method.',k);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% %Xkacz = kaczmarz(A,b,k);
% [Xkacz_prop_non,normalized_variance_kaczmarz_constrain] = constrain_kaczmarz_cancha(A,b,k_prop,x_ex);
% %Show the kaczmarz solution.
% figure
% imagesc(reshape(Xkacz_prop_non,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('Kaczmarz constrain reconstruction')
% hold on
% 
% k = toc;
% k
% plot(1:length(normalized_variance_kaczmarz_constrain), normalized_variance_kaczmarz_constrain, 'b-', 'DisplayName', '正则化');
% hold on;
%-----------------------------------------------kaczmarz_max_3-----------------------------%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_3''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max = 10001;
% [Xkacz_max_5,normalized_variance_kaczmarz_max_5] = greedy(A,b,k_max,x_ex);
% figure
% imagesc(reshape(Xkacz_max_5,N,N)), colormap gray,
% axis image off
% %caxis(c);
% title('Kaczmarz greedy reconstruction')


% plot(1:length(normalized_variance_kaczmarz_max_5),normalized_variance_kaczmarz_max_5);
% hold on;
% [Xkacz_max_3,normalized_variance_kaczmarz_max_3] = kaczmarz_max_3(A,b,k_max,x_ex,Xkacz_max_5);
% % Show the kaczmarz solution.
% figure
% imagesc(reshape(Xkacz_max_3,N,N)), colormap gray,
% axis image off
% % caxis(c);
% title('Kaczmarz max 3 reconstruction')
% toc
% k = toc;
% k
% plot(1:length(normalized_variance_kaczmarz_max_3),normalized_variance_kaczmarz_max_3);
% hold on;
%-----------------------------------------------kaczmarz_max_4-----------------------------%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_4''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max = 10000;
% [Xkacz_max_4,normalized_variance_kaczmarz_max_4] = kaczmarz_max_4(A,b,k_max,x_ex);
% % Show the kaczmarz solution.
% figure
% imagesc(reshape(Xkacz_max_4,N,N)), colormap gray,
% axis image off
% % caxis(c);
% title('Kaczmarz max 4 reconstruction')
% toc
% k = toc;
% k
% 
% %-----------------------------------------------kaczmarz_max_5-----------------------------%
% tic
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with Kaczmarz_max_5''s method.',k,'LineWidth', 2);
% fprintf(1,'\nThis takes a moment ...');
% %Perform the kaczmarz iterations.
% k_max = 40000;
% [Xkacz_max_5,normalized_variance_kaczmarz_max_5] = kaczmarz_max_5(A,b,k_max,x_ex);
% %Show the kaczmarz solution.
% figure
% imagesc(reshape(Xkacz_max_5,N,N)), colormap gray,
% axis image off
% % caxis(c);
% title('Kaczmarz max 5 reconstruction')
% toc
% k = toc;
% k
% plot(1:length(normalized_variance_kaczmarz_max_4),normalized_variance_kaczmarz_max_4);
% hold on;
% plot(1:length(normalized_variance_kaczmarz_max_5),normalized_variance_kaczmarz_max_5);
% hold on;
%  legend('kaczmarz max 4','kaczmarz max 5');
