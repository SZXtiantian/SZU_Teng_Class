load ./MatrixA.mat;
load ./MatrixB.mat;
% A = [1, 2, 1, 2; 1, 6, -1,3; -1, -6, 0,0;-1,-2,0,-1];
[Q, R] = qr_householder(A);
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

function [Q, R] = qr_householder(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for k = 1:min(m-1, n)
        v = R(k:m, k);
        v(1) = v(1) + sign(v(1)) * norm(v);
        v = v / norm(v);

        H_k = eye(m);
        H_k(k:end, k:end) = H_k(k:end, k:end) - 2 * v * v';

        Q = Q * H_k';
        R = H_k * R;
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
    figure;
    plot(E);

end
