%列残差最大行+A :K-means
%只选择一列残差最大行+k-means



function [s] = zhizuidak_means(A,B,x,x0,tol,h,k_fenlei,kb,maxit)
[m,n]=size(A);
idx = new_set(A,k_fenlei);
uniqueValues = unique(idx);
% 初始化元胞数组
 % 遍历唯一值，并将相同值的索引存放在元胞数组中
for i = 1:numel(uniqueValues)
    indexCellArray{i} = find(idx == uniqueValues(i));
end



ta=clock;
x0=zeros(n,kb);%kb个右端项
a=(sum(A.^2,2)).^(1/2);%计算系数矩阵A每行的二范数
maxit=1e+7;
res=1;
iter=0;
while res>=tol&&iter<=maxit
    tic
    % 遍历元胞数组，随机选择一个值
for i = 1:numel(indexCellArray)
    % 如果当前元胞数组非空，则随机选择一个值
    if ~isempty(indexCellArray{i})
        selectedIndex = randi(length(indexCellArray{i}));
        ik(i) = indexCellArray{i}(selectedIndex);
    end
end 
    g =abs(B -A*x0);%求残差
    [~, top_three_indices] = maxk(g(:), h);
    [index_1, ~] = ind2sub(size(g), top_three_indices);
    [index_1, ~, ic] = unique(index_1);
    
    index = unique([ik'; index_1]);%并集
    %index = intersect(ik, top_three_indices_1);%交集
    Ak=A(index,:);%将指标构成块进行迭代
    Bk=B(index,:);
    iter=iter+1;
    ak=(sum(Ak.^2,2));

%     C=(Bk-Ak*x0)./repmat(ak,1, size(Ak,2));
%     z = (1/h)*Ak'* C ;
%     a4 = Ak';
%     x0=x0+z;%平均块的更新，x=x+求和（w(bi-ai*x/f:ai)）
    C = Bk-Ak*x0;
    [U, S, V] = svd(Ak);
    S_plus = zeros(size(S));
for i = 1:rank(Ak)
    S_plus(i,i) = 1/S(i,i);
end

    % 计算伪逆
    A_plus = V * S_plus' * U';
    x0 = x0 + A_plus*C;
  time(:,iter) = toc;
    RES=sum((x-x0).^2,1)./sum(x0.^2,1);
   
    %RES = norm(B-A*x);
   %RES=(sum((x-x0).^2,1)./sum(x0.^2,1)).^1/2;
    %for i=1:kb
   % RES(i)=norm(x(:,i)-x0(:,i))/norm(x0(:,i));
    %end
  
%              errors_matrix = (x - x0).^2;

        % 计算原始数据矩阵的方差
%         original_variance = var(x(:));
%         original_variance =original_variance.^2;
%         % 计算误差矩阵的方差
%         errors_variance = var(errors_matrix(:));
%        
%         % 计算归一化方差
%         facha = (sum(errors_variance) / sum(original_variance))^(0.5);
%         pre(:,iter) = facha;
%         ssimValue(:,iter) = ssim(x, x0);
%         psnrValue(:,iter) = psnr(x0,x,255);
        res=max(RES);
       
     error(:,iter)= res;
        
    
%     threshold = 1e-5;
%     index_set = find(max(A, [], 1) < threshold);%选择每列残差最大值小于某值的列索引
end
tb=clock;
cup=etime(tb,ta);
% x_uint8 = uint8(x0);
% figure; % 创建一个新的图形窗口
% 
% subplot(1, 2, 1); % 创建第一个子图，占一行两列，在第一列
% imshow(img);
% title('原图');
% 
% subplot(1, 2, 2); % 创建第二个子图，在第二列
% imshow(x_uint8);
% title('重建后');
s(:,1)=cup;
s(:,2)=iter;
% s(:,3)=res;
% s(:,4)=size(index,1);
end

function [idx] = new_set(A,k_fenlei)
normalizedA = normalize(A, 2, 'norm');
uniqueFirstValues = unique(normalizedA(:, 1));
selectedRows = [];
for i = 1:length(uniqueFirstValues)
    index = find(normalizedA(:, 1) == uniqueFirstValues(i), 1);
    selectedRows = [selectedRows; A(index, :)];
    if size(selectedRows, 1) >= k_fenlei
        break; % 当找到十行时退出循环
    end
end
% 使用 selectedRows 作为初始中心
selectedRows = full(selectedRows);

[idx, centers] = kmeans(A, k_fenlei, 'Start', selectedRows,'MaxIter',50);

% 显示结果
% disp('K-means 聚类结果:');
% disp(idx);
% disp('聚类中心:');
% disp(centers);
end