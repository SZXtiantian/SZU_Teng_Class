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
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for j = 1:n
        v = A(:, j);

        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end

function check(Q)
    %%正交性偏差验证
    n = size(Q,1);
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
