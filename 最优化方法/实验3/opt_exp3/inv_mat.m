load ./MatrixB.mat;
inv_matrix = inv_matr(B);
% temp = inv_matrix * B;
% I = abs(eye(size(B, 1)) - inv_matrix * B);
% figure;
% h = heatmap(temp);
% h.CellLabelFormat = '%0.2f';
% h.CellLabelColor = 'none';
% B_inv = inv(B);
% bias = abs(inv_matrix - B_inv);
% figure;
% h = heatmap(bias);
% h.CellLabelFormat = '%0.2f';
% h.CellLabelColor = 'none';


function m = inv_matr(A)
    tic;
    m = [];
    [Q, R] = qr_householder(A);
    % [Q, R] = schmidt_qr(A);
    diagonal_elements = diag(A);
    product_of_diagonal = prod(diagonal_elements);
    if product_of_diagonal == 0
        disp("矩阵不可逆！")
    else
        Q = Q';
        for i = 1 : size(R, 1)
            [x] = Gback(R, Q(:, i), size(R, 1));
            m = [m, x];
        end     
    end
    elapsedTime = toc;
    disp(['代码执行时间：', num2str(elapsedTime), ' 秒']);
end


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
function [x]=Gback(A1, b1, n)
    for i = n: -1: 1
        if(A1(i, i) ~= 0)
            x(i) = b1(i);
            for j=n: -1 :i + 1
                x(i) = x(i) - A1(i, j) * x(j);
            end
            x(i) = x(i) / A1(i, i);
        else
            disp("error!");
            return;
        end
    end
    x = x';
end


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