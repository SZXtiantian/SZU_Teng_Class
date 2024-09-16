load ./FittingData.mat;
% figure;
% plot(x, y);
% title('a示例');
% xlabel('X轴标签');
% ylabel('Y轴标签');
% grid on;  % 添加网格线
A = getA(x, 2);
a = kktfun(A, y', [1,1,1], 4);
a1 = kktQR(A, y', [1,1,1], 4);
a2 = lsqlin(A, y',[],[],[1,1,1],4);
disp(a);
disp(a1);
disp(a2);
disp(calerror(x, y, a));
disp(calerror(x, y, a1));
show(x, y, a);
show(x, y, a1);

function show(x, y, a)
    figure;
    y1 = a(3,1)*x.^2 + a(2,1)*x + a(1,1);
    plot(x, y1);
    hold on;
    scatter(x, y, 'filled', 'MarkerFaceColor', 'r');
    % 添加标题和轴标签
    title('Graph');
    xlabel('x');
    ylabel('y');
end

function [A] = getA(x, mi)
    col = size(x, 2);
    A = zeros(col, (mi+1));
    for i = 1 : col
        A(i, 1) = 1;
        for j = 2 : (mi+1) 
            A(i, j) = A(i, j - 1) * x(i);
        end
    end
end
 
function [x] = kktfun(A, b, C, d)
    tic;
    gram = 2 * (A' * A);
    ATb = 2 * (A' * b);
    row = size(C, 1);
    zeroMat = zeros(row, row);
    matrix = [gram, C';C, zeroMat];
    vec = [ATb; d];
    x = matrix \ vec;
    elapsed_time = toc;
    fprintf('time: %.4f seconds\n', elapsed_time)
end


function [x] = kktQR(A, b, C, d)
    
    rowA = size(A, 1);
    colA = size(A, 2);
    AC=[A; C];
    tic;

    [~, R] = qr(AC);
    R = R(1: colA,1: colA);
    I= eye(size(R, 1));
    Q = AC * (R \ I);
    Q1 = Q(1:rowA, :);
    Q2 = Q(rowA + 1:size(Q, 1), :);
    [Qtuta,Rtuta] = qr(Q2');
    temp = Rtuta' \ d;
    w = Rtuta \(2 * (Qtuta' * Q1' *b)-2 *temp);
    x = R \ (Q1' * b - 0.5 * Q2' *w);
    elapsed_time = toc;
    fprintf('time: %.4f seconds\n', elapsed_time)
end

function [norm_v] = calerror(x, y, a)
    num = size(x, 2);
    error = zeros(num, 1);
    for i = 1 : num
        error(i, 1) = a(3,1)*x(1, i)^2 + a(2,1)*x(1, i) + a(1,1) - y(1, i);
    end
    norm_v = norm(error, 2);
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