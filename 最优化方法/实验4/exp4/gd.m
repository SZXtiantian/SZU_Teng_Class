load ./Matrix_A_b.mat;
row = size(A,1);
col = size(A,2);

a = normal_linear(A, b);
figure;
plot(1:40, a,'-or');
title('精确解');
xlabel('第i个解');
ylabel('大小');
grid on;  % 添加网格线
% b1 = normal_linear_1(A, b);

% total = 0;
% n = 1;
% for i = 1: n
%     [cnt, c] = gd_faseter(A, b, a);
%     % total = total + cnt;
% end


% total = total/n;
% disp("平均："+total);
% [d, p] = adaGrad(A, b, a);
e = temp(A,b);
% disp("p")
% disp(p);
% disp([d,a,e,c]);
% disp([A*p - b,A*a - b,A*e - b]);

function [beta]=normal_linear(X, y)
    beta = (X' * X) \ (X' * y);

end

function [x_least]=normal_linear_1(X, y)
    [Q,R] = qr(X);
    x_least=R\(Q'*y); 
end

function [cnt, x]=gd_faseter(A, b, gt)
    n = size(A, 2);
    % x = randn(n, 1);
    x = zeros(n, 1);
    cnt1 = x;
    n = 600;
    for j = 1 : 6
        % 
        cnt = 0;
        % times1 = n * j;
        disp(j);
        for i = 1 : 100
            % f1(1,i) = norm(gt - x, 2);
            p = A' * (A * x - b);
            norm_p = norm(p, 2);
            norm_p = norm_p * norm_p;
            norm_Ap = norm(A * p, 2);
            norm_Ap = norm_Ap *norm_Ap;
            Yita = norm_p / norm_Ap;
            x1 = x - Yita * p;
            % f2(1, i) = Yita;
            % f(1,i) = norm(x - x1, 2)/norm(x, 2);
            x = x1;
            % if f(1,i) < 0.001
            %     % disp("迭代次数：" + i);
            %     cnt = i;
            %     n = i;
            %     break;
            % end
            % if f1(1,i) < 0.001
            %    cnt = i;
            %    n = i;
            %    break;
            % end
        end
        disp(size(x));
        disp(size(cnt1))
        cnt1 = [cnt1, x];
    end
    % figure;
    % plot(1:n, f);
    % title('迭代次数与迭代解间的差异的关系');
    % xlabel('迭代次数');
    % ylabel('迭代解间的差异');
    % grid on;  % 添加网格线
    % figure;
    % plot(1:n, f1);
    % title('最速下降法迭代次数与误差关系');
    % xlabel('迭代次数');
    % ylabel('误差关系');
    % grid on;  % 添加网格线
    % figure;
    % plot(1:n, f2);
    % title('Yita值与迭代次数的关系');
    % xlabel('X轴标签');
    % ylabel('Y轴标签');
    % grid on;  % 添加网格线
    % figure;
    % plot(1:40, [x, gt]);
    % title('近似解');
    % xlabel('第i个');
    % ylabel('大小');
    % grid on;  % 添加网格线
    Dim = 1: 40;
    for k = 1 : 6
        t = 100 *k ;
        str = ['近似解与精确解(迭代次数',num2str(t),')'];
        subplot(2,3,k);
        plot(Dim,cnt1(:, k + 1),'-*b',Dim,gt,'-or'); %线性，颜色，标记
        title(str);
        axis( [0,41,-0.8,0.8])  %确定x轴与y轴框图大小
        set(gca,'XTick',[0:40]) %x轴范围
        set(gca,'YTick',[-1:0.2:1]) %y轴范围
        legend('近似值','准确值','Location','SouthEast');   %右下角标注
        xlabel('第i个解')  %x轴坐标描述
        ylabel('大小') %y轴坐标描述
    end
    

    

end

function [x, p]= adaGrad(A, b, gt)
    ep = 1e-6;
    alpha = 0.01;
    n = size(A, 2);
    x = zeros(n, 1);
    g = zeros(n,1);
    epsilon = ones(n,1) *ep;
    

    for i = 1 : 15000
        f1(1,i) = norm(gt - x, 2);
        p = A' * (A * x - b);
        g = g + p.^2;
        g1 = g ./ i;
        Yita = alpha ./ sqrt(g1 + epsilon);
        x1 = x - Yita .* p;
        f(1,i) = norm(x - x1, 2)/norm(x, 2);
        x = x1;
    end
    figure;
    plot(1:15000, f);
    title('adaGrad折线图示例');
    xlabel('X轴标签');
    ylabel('Y轴标签');
    grid on;  % 添加网格线
    figure;
    plot(1:15000, f1);
    title('adaGrad迭代次数与误差关系');
    xlabel('X轴标签');
    ylabel('Y轴标签');
    grid on;  % 添加网格线
end

function [x]=temp(A,b)
    
    min=0.01;
    x=zeros(40,1);
    for k = 1:30 %或指定迭代次数
        f(1,k)=0.5*norm(A*x-b,2)^2; % 目标函数值
        p = A'*(A*x-b);
        a = norm(p,2)^2 / norm(A*p,2)^2;
        y = x - a * p; %y为x（k+1）
        % temp1(1,k) = norm((x-y),2)/norm(x,2); %迭代解间的相对接近程度
        % error(1,k) = norm((x_least - x),2); %误差迭代
        % 
        %   if norm((x-y),2)/norm(x,2) < min
        %       break
        %   end
        x = y; %迭代
    end

end
