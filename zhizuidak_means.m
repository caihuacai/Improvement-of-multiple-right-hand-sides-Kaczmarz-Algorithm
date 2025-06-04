%�вв������+A :K-means
%ֻѡ��һ�вв������+k-means



function [s] = zhizuidak_means(A,B,x,x0,tol,h,k_fenlei,kb,maxit)
[m,n]=size(A);
idx = new_set(A,k_fenlei);
uniqueValues = unique(idx);
% ��ʼ��Ԫ������
 % ����Ψһֵ��������ֵͬ�����������Ԫ��������
for i = 1:numel(uniqueValues)
    indexCellArray{i} = find(idx == uniqueValues(i));
end



ta=clock;
x0=zeros(n,kb);%kb���Ҷ���
a=(sum(A.^2,2)).^(1/2);%����ϵ������Aÿ�еĶ�����
maxit=1e+7;
res=1;
iter=0;
while res>=tol&&iter<=maxit
    tic
    % ����Ԫ�����飬���ѡ��һ��ֵ
for i = 1:numel(indexCellArray)
    % �����ǰԪ������ǿգ������ѡ��һ��ֵ
    if ~isempty(indexCellArray{i})
        selectedIndex = randi(length(indexCellArray{i}));
        ik(i) = indexCellArray{i}(selectedIndex);
    end
end 
    g =abs(B -A*x0);%��в�
    [~, top_three_indices] = maxk(g(:), h);
    [index_1, ~] = ind2sub(size(g), top_three_indices);
    [index_1, ~, ic] = unique(index_1);
    
    index = unique([ik'; index_1]);%����
    %index = intersect(ik, top_three_indices_1);%����
    Ak=A(index,:);%��ָ�깹�ɿ���е���
    Bk=B(index,:);
    iter=iter+1;
    ak=(sum(Ak.^2,2));

%     C=(Bk-Ak*x0)./repmat(ak,1, size(Ak,2));
%     z = (1/h)*Ak'* C ;
%     a4 = Ak';
%     x0=x0+z;%ƽ����ĸ��£�x=x+��ͣ�w(bi-ai*x/f:ai)��
    C = Bk-Ak*x0;
    [U, S, V] = svd(Ak);
    S_plus = zeros(size(S));
for i = 1:rank(Ak)
    S_plus(i,i) = 1/S(i,i);
end

    % ����α��
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

        % ����ԭʼ���ݾ���ķ���
%         original_variance = var(x(:));
%         original_variance =original_variance.^2;
%         % ����������ķ���
%         errors_variance = var(errors_matrix(:));
%        
%         % �����һ������
%         facha = (sum(errors_variance) / sum(original_variance))^(0.5);
%         pre(:,iter) = facha;
%         ssimValue(:,iter) = ssim(x, x0);
%         psnrValue(:,iter) = psnr(x0,x,255);
        res=max(RES);
       
     error(:,iter)= res;
        
    
%     threshold = 1e-5;
%     index_set = find(max(A, [], 1) < threshold);%ѡ��ÿ�вв����ֵС��ĳֵ��������
end
tb=clock;
cup=etime(tb,ta);
% x_uint8 = uint8(x0);
% figure; % ����һ���µ�ͼ�δ���
% 
% subplot(1, 2, 1); % ������һ����ͼ��ռһ�����У��ڵ�һ��
% imshow(img);
% title('ԭͼ');
% 
% subplot(1, 2, 2); % �����ڶ�����ͼ���ڵڶ���
% imshow(x_uint8);
% title('�ؽ���');
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
        break; % ���ҵ�ʮ��ʱ�˳�ѭ��
    end
end
% ʹ�� selectedRows ��Ϊ��ʼ����
selectedRows = full(selectedRows);

[idx, centers] = kmeans(A, k_fenlei, 'Start', selectedRows,'MaxIter',50);

% ��ʾ���
% disp('K-means ������:');
% disp(idx);
% disp('��������:');
% disp(centers);
end