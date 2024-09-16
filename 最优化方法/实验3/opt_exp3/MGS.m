load ./MatrixA.mat;
load ./MatrixB.mat;
% A = [1, 2, 1, 2; 1, 6, -1,3; -1, -6, 0,0;-1,-2,0,-1];
[Q, R] = schmidt_qr(A);
I1 = Q' *  Q;
h = heatmap(I1);
h.CellLabelFormat = '%0.2f';
h.CellLabelColor = 'none';
column_sum_abs = sum(abs(I1));
figure;
plot(column_sum_abs, '-o'); % 使用'-o'表示用圆点连接数据点
title('每列绝对值求和结果');
xlabel('列索引');
ylabel('绝对值求和');
grid on; % 显示网格


disp('Q:');
disp(Q);
disp('R:');
disp(R);
disp(I1);
check(Q);

function [Q, R] = schmidt_qr(A)
    A1 = A;
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    Q(:,1)=A(:,1)/norm(A(:,1)); %定义第一个向量基准
    R(1,1)=norm(A(:,1));
    for k=2:n
         for j = 1: k - 1
             R(j, k) = Q(:, j)' * A1(:, k);
         end

         for i=k:n
             A(:,i)=A(:,i)-A(:,i)'*Q(:,k-1)*Q(:,k-1); 
         end


         R(k,k)=norm(A(:,k)); %求出R（K,K）
         Q(:,k)=A(:,k)/R(k,k); %求出标准正交化向量qk
    end
end


function check(Q)
    %%正交性偏差验证
    n = size(Q,1)
    E = zeros(1,n);
     
    for k=2:n
        max = 0;
        for i=1:k-1
            temp = abs(Q(:,i)' *  Q(:,k));
            if temp > max
                max = temp;
            end
        end
        E(1,k)=max;
    end    
    plot(E)

end
    
