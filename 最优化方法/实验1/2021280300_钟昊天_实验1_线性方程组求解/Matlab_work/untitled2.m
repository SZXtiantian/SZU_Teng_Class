% A = input("输入系数矩阵：");
%b = input("输入常数向量b:");
Test = load("A_b.mat");
%main(A, b');
% GSel(A, b',size(A, 1));
main(Test.A, Test.b);
%inv_A(A);
% Y = round(0.998,2);
% disp(Y);

% 检查A，b
% 符合条件   marker返回1和系数矩阵维数为n
% 不符合条件 marker返回0，n为-1
function [marker, n]=check(A, b)
    % 1.系数矩阵是否为方阵
    if size(A, 1) == size(A, 2)

        % 2.常数向量b输入是否为n行1列
        if size(b, 1) == size(A, 1) && size(b, 2) == 1 

            % 3.行列式是否为0
            if(det(A) ~= 0)
                n = size(A, 1);
                marker = 1;
                return;
            else 
                disp("行列式等于0，无解！");
                marker = 0;
                n = -1;
                return;
            end
        else
            disp("常数向量b输入不合法！");
            marker = 0;
            n = -1;
            return;
        end
    else
        disp("输入的系数矩阵A不是方阵！");
        marker = 0;
        n = -1;
        return;
    end

end

% 主元需要是较大值
function [A_temp, b_temp] = findMax(A, b, k, n)
    temp = A(k, k);
    mark = k;
    for i = (k + 1):n
        if(A(i, k) ~= 0 && abs(A(i, k) ) > abs(temp) )
            temp = A(i, k);
            mark = i;
        end
    end
     A([k,mark],:) = A([mark,k],:);
     b([k,mark],:) = b([mark,k],:);
     A_temp = A;
     b_temp = b;
end

% 高斯消元法，下三角
function [A1, b1]=GSel(A, b, n)

    for k=1:(n - 1)
       % 主元非0
       if(A(k, k) ~= 0)
           [A, b] = findMax(A, b, k, n);
           % disp(A);
           % 在第k列中将k+1到n行的元素置0
           for i=(k+1):n

                % 计算系数c
                c=(-1 * A(i, k)/A(k,k));

                % 行变换
                for j=1:n
                    A(i, j) = A(i, j) + c * A(k, j);
                end
                % 理论上能看作增广矩阵，但这里分开运算
                b(i) = b(i) + c * b(k);
                % disp(A);
           end
       
       % 主元为0，在该列向下找到第一个非零元素的行，然后两行交换
       else

           % m = -1 说明下面元素全为0，
           % 否则 m 为向下第一个非0元素所在的行
           m = FindFirstNonZero(A, k, n);

           % 找不到说明无解
           if(m == -1)
               disp("error!");            
               return; 

           % 找到了交换两行然后消元
           else
               A([k,m],:) = A([m,k],:);
               b([k,m],:) = b([m,k],:);
               for i=(k+1):n                 
                   c=(-1 * A(i, k)/A(k,k));                 
                   for j=1:n                     
                       A(i, j) = A(i, j) + c * A(k, j);
                   end                 
                   b(i) = b(i) + c * b(k); 
               end
           end     
       end

    end

    % 返回A1, b1
    A1 = A;
    b1 = b;
end

function [A1, b1]=GSel2(A, b, n)

    for k=1:(n - 1)
       [A, b] = findMax(A, b, k, n);
       for i=(k+1):n
                % 计算系数c
                c=(-1 * A(i, k)/A(k,k));
                % 行变换
                for j=1:n
                    A(i, j) = A(i, j) + c * A(k, j);
                end
                % 理论上能看作增广矩阵，但这里分开运算
                for j=1:n
                    b(i, j) = b(i, j) + c * b(k, j);
                end
                % disp(A);
       end
    end
    % 返回A1, b1
    A1 = A;
    b1 = b;
end

function [A1, b1]=GSel1(A, b, n)

    for k=1:(n - 1)
       [A, b] = findMax(A, b, k, n);
       for i=(k+1):n
                % 计算系数c
                c=(-1 * A(i, k)/A(k,k));
                % 行变换
                for j=1:n
                    A(i, j) = A(i, j) + c * A(k, j);
                end
                % 理论上能看作增广矩阵，但这里分开运算
                b(i) = b(i) + c * b(k);
                % disp(A);
       end
    end
    % 返回A1, b1
    A1 = A;
    b1 = b;
end


% 在该列向下找到第一个非零元素的行 m 并返回
function [m]=FindFirstNonZero(A, k, n)
    for i=k+1:n
        if(A(i, k)~=0)
            m = i;
            return;
        end
    end
    m = -1;
end


% 回代，从下至上
function [x]=Gback(A1, b1, n)
    for i = n: -1: 1
        if(A1(i, i) ~= 0)
            x(i) = b1(i);
            for j=n: -1 :i + 1
                x(i) = x(i) - A1(i, j) * x(j);
            end
            x(i) = x(i) / A1(i, i);
            %x(i) = round(x(i), 3);
        else
            disp("error!");
            return;
        end
    end
    x = x';
end


% 主函数
function main(A, b)     
    [marker, n] = check(A, b);     
    if marker == 1         
        [A1, b1] = GSel1(A, b, n);
        [x] = Gback(A1, b1, n);
        disp("高斯消元法求解");
        disp("x=");
        disp(x);
        x1 = A\b;
        disp("matlab函数求得精确解(基准)：");
        disp("a=");
        disp(x1);
        disp("绝对误差(x - a)：")
        temp = x - x1;
        disp(x - x1);
        % disp("MSE=");
        % disp(mse(x - x1))
        disp("MSE=");
        disp(mean((x - x1).^2))
        plot(1:100,temp,'-*b');
        xlabel('x-第i个解')  %x轴坐标描述
        ylabel('y-绝对误差') %y轴坐标描述


        %scatter(1:100, [x';x1']);
        % plot(1:100, [x';x1'],'o','g*')
        %plot(1:100 ,[x;x1],'bo',x(1:3),x1(1:3),'g*')
        %x2=-10:1:10; 
        %plot(x2,x2,'r') 
        %grid on 
        %hold on 
        %scatter(x,x1,[],"+")
        %xlabel("matlab函数求得精确解")
        %ylabel("高斯消元法求的解")
        %legend("y=x")
        %scatter(1:100,x1,100,"+")
        %hold on
        %scatter(1:100,x,100)
        %xlabel("x-第i个解")
        %ylabel("y-第i个解的值")
        %legend("高斯消元法求的解","matlab函数求得精确解")
        %set(gca, 'FontSize', 14);
        
    else
        disp("error!");
    end  
end

% 求逆矩阵
function inv_A(A)
    m = [];
    if size(A, 1) == size(A, 2)
        if det(A) ~= 0
            I = eye(size(A, 1));
            [A1, b1] = GSel2(A, I, size(A, 1));
            for i = 1 : size(A, 1)           
                [x] = Gback(A1, b1(:, i), size(A, 1));
                m = [m, x];
            end
            disp("高斯消元法求得的逆=")
            disp(m)
            disp("matlab函数求得逆(基准)：");
            disp(inv(A))
            disp("对应元素的绝对误差：")
            disp(m - inv(A))
        else
            disp("行列式为0，不可逆！");
            return;
        end
    else
        disp("不是方阵，不可逆！");
        return;
    end
end


