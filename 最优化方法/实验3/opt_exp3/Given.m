load ./MatrixA.mat;
load ./MatrixB.mat;
% A = [1, 2, 1, 2; 1, 6, -1,3; -1, -6, 0,0;-1,-2,0,-1];
[Q, R] = qr_givens(A);
I1 = Q' *  Q;
h = heatmap(I1);
h.CellLabelFormat = '%0.2f';
h.CellLabelColor = 'none';

check(Q);

function [Q, R] = qr_givens(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for k = 1:n
        for i = m:-1:k+1
            % 计算Givens旋转矩阵
            [c, s] = givens(R(i-1, k), R(i, k));

            % 更新R
            G = [c -s; s c];
            R(i-1:i, k:end) = G' * R(i-1:i, k:end);

            % 更新Q
            Q(:, i-1:i) = Q(:, i-1:i) * G;
        end
    end
end

function [c, s] = givens(a, b)
    % 计算Givens旋转参数
    if b == 0
        c = 1;
        s = 0;
    else
        r = hypot(a, b);
        c = a / r;
        s = b / r;
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