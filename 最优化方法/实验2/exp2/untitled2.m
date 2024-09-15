% data = rand(100, 2);
filename = 'example.xlsx';
sheet = 'Sheet1';
result = [];
data = [];
labels = [];

times = 10;
% 设置聚类数量
for type = 1 : 1
    disp("type:" + type);
    for t =1 : 6
        sum = 0;
        tic;
    
        numPic = 10000*t;
        k = 10;
        for i = 1 : times
        % randomOrtVec( );
            [data, labels] = DataRead(numPic);
            % % 运行 K-means 算法
            [~,~,centroids, label] = kmeans_algorithm(data, k, labels, type);
            % [centroids, label] = bikmeans_algorithm(data, k, labels);
            cnt = countRight(labels, label', k);
            sum = sum + cnt;
            % sum = sum + countRight1(label, k, numPic);
            result(i) = cnt;
        end
        elapsedTime = toc;
        disp(['代码执行时间：', num2str(elapsedTime/times), ' 秒']);
        disp(sum / times);
    end
end
% disp([labels', label])
% if exist(filename, 'file')
%     existingData = xlsread(filename, sheet);
% else
%     existingData = [];
% end
% 
% % 将新数据追加到现有数据后面
% combinedData = [existingData; result];
% 
% % 使用 xlswrite 函数将合并后的数据写入 Excel 文件
% xlswrite(filename, combinedData, sheet);
% plot(1:times, result,'o-','DisplayName', '相互最远向量','Color',"blue");
% % disp(result);
% legend('Location', 'northeast');


function [correct] = countRight1(label,k, numPic)
    totalNum = numPic;
    len = k;
    correctSum = 0;
    for i = 1 : len
        vector = label(i);
        most_frequent_element = mode(vector{1});
        count_most_frequent = sum(vector{1} == most_frequent_element);
        correctSum = correctSum + count_most_frequent;
        % disp(correctSum);
    end
    % disp(correctSum);
    % disp(correctSum/totalNum);
    correct = correctSum/totalNum;
end

function [correct] = countRight(labels, label, k)
    totalNum = size(label, 2);

    cnt = zeros(k, k + 1);
    for i = 1 : totalNum
        cnt(labels(i) + 1, label(i)) = cnt(labels(i) + 1, label(i)) + 1;  
    end
    correctSum = 0;
    for i = 1 : k
        maxValNum = max(cnt(i, :));
        correctSum = correctSum + maxValNum;
    end
    % disp(correctSum);
    % disp(correctSum/totalNum);
    correct = correctSum/totalNum;
end

function [orthogonal_vectors] = randomOrtVec()
    num_vectors = 10;
    vector_dim = 784;
    vectors = randn(vector_dim, num_vectors);
    orthogonal_vectors = zeros(vector_dim, 10);

    for i = 1:10
        vector = vectors(:, i);
        for j = 1:i - 1
            coefficient = dot(vector, orthogonal_vectors(:, j)) / dot(orthogonal_vectors(:, j), orthogonal_vectors(:, j));
            vector = vector - coefficient * orthogonal_vectors(:, j);
        end
        orthogonal_vectors(:, i) = vector;
    end
    % disp(size(orthogonal_vectors))
end


function [orthogonal_vectors] = randomOrtVecFromData(data,k)
    vector_dim = 784;
    vectors = data( :,randperm(size(data, 2), k));
    orthogonal_vectors = zeros(vector_dim, 10);

    for i = 1:10
        vector = vectors(:, i);
        for j = 1:i - 1
            coefficient = dot(vector, orthogonal_vectors(:, j)) / dot(orthogonal_vectors(:, j), orthogonal_vectors(:, j));
            vector = vector - coefficient * orthogonal_vectors(:, j);
        end
        orthogonal_vectors(:, i) = vector;
    end
    disp(size(orthogonal_vectors))
end

function [data, labels]=DataRead(numPic)
    load ./train_images.mat;    % Read image data
    load ./train_labels.mat;     %Read lable of images
    data = [];
    labels = [];
    ImgNum = numPic; 
    for i = 1:ImgNum
        GetOneImg = train_images(:,:,i);
        v_GetOneImg = GetOneImg(:);
        data = [data,v_GetOneImg];
        labels = [labels,train_labels(i)];
    end
end

function [clusters, be_labels, centroids, labels] = kmeans_algorithm(data, k, be_label, type)
    cnt = 0;
    clusters{1} = data;
    be_labels{1} = be_label;
    % 初始化聚类中心
    % centroids = data( :,randperm(size(data, 2), k));
    % centroids = randomOrtVecFromData(data,k);
    % centroids = randomOrtVec();
    if type == 1
        centroids = kmeans_init(data', k);
    elseif type == 2
        centroids = randomOrtVec();
    elseif type == 3
        centroids = data( :,randperm(size(data, 2), k));
    end
    % 初始化标签
    labels = zeros(size(data, 2), 1);
    % disp(labels');
    % 迭代更新聚类中心和标签直到收敛
    converged = false;
    while ~converged
        cnt = cnt + 1;
        % 计算每个样本与聚类中心的距离
        distances = pdist2(data', centroids');
        
        % 分配每个样本到最近的聚类中心
        [~, new_labels] = min(distances, [], 2);


        % disp(new_labels);
        % 如果标签没有变化，则收敛
        if isequal(labels, new_labels)
            converged = true;
        else
            labels = new_labels;
            
            % 更新聚类中心为每个簇的均值
            for i = 1:k
                % disp(mean(data( :,labels == i), 2));
                centroids( :,i) = mean(data( :,labels == i), 2);
                clusters{i} = data( :,labels == i);
                be_labels{i} = be_label(labels == i);

            end
        end
    end
    % disp(clusters{k});
    % disp(be_labels{k});
end

function initial_centers = kmeans_init(data, K)
    [N, ~] = size(data);
    initial_centers = zeros(K, size(data, 2));
    
    % 随机选择第一个聚类中心
    initial_centers(1, :) = data(randi(N), :);
    
    % 选择剩余的聚类中心
    for k = 2:K
        min_distances = inf(N, 1);
        for i = 1:k-1
            % 计算每个数据点到最近的聚类中心的距离
            distances = sum((data - initial_centers(i, :)).^2, 2);
            min_distances = min(min_distances, distances);
        end
        % 使用加权随机采样选择下一个聚类中心
        probabilities = min_distances / sum(min_distances);
        cumulative_probabilities = cumsum(probabilities);
        r = rand();
        for i = 1:N
            if r <= cumulative_probabilities(i)
                initial_centers(k, :) = data(i, :);
                break;
            end
        end
    end
    initial_centers = initial_centers';
end

 
function [cluster_to_split, be_labels,cluster_to_split_index] = chooseClusterToSplit(cluster_data, be_label)
    % 假设已经有一个簇的数据，存储在 cluster_data 变量中
    num_clusters = length(cluster_data);
    % disp(num_clusters);
    variances = zeros(num_clusters, 1);
    
    for i = 1:num_clusters
        cluster = cluster_data{i};
        sum1 = 0;
        % disp(cluster);
        len = size(cluster, 2);
        column = sum(cluster, 2);
        % variances(i) = norm(var(cluster,1,1));  % 计算每个簇内数据点的方差
        for j = 1 : len
            sum1 = sum1 + norm((cluster(j) - column));

        end
        variances(i) = sum1;
    end
    % disp(variances');
    % 选择具有最大方差的簇
    [~, cluster_to_split_index] = max(variances);
    cluster_to_split = cluster_data{cluster_to_split_index};
    % disp(cluster_to_split_index);
    be_labels = be_label{cluster_to_split_index};
end

function [clusters, be_label] = bikmeans_algorithm(data, k, label)
    clusters{1} = data;
    be_label{1} = label;
    centroid = zeros(size(data, 1), k);
    for i = 1:k-1
        
        [cluster_to_split, inputlabels, be] = chooseClusterToSplit(clusters, be_label);

        [subclusters, be_labels, centroids, ~] = kmeans_algorithm(cluster_to_split, 2 ,inputlabels);

        clusters{be} = subclusters{1};
        be_label{be} = be_labels{1};
        centroid(:,be)=centroids(:,1);
        clusters{i+1} = subclusters{2};
        be_label{i+1} = be_labels{2};
        centroid(:,i + 1)=centroids(:,2);
        % disp(clusters);
    end
    % disp(centroid);
    [clusters, be_label]=kmeans_algorithm1(data, k ,label, centroid);
end

function [centroids, be_labels] = kmeans_algorithm1(data, k, be_label,centroid)
    cnt = 0;
    clusters{1} = data;
    be_labels{1} = be_label;
    % 初始化聚类中心
    centroids = centroid;
    % centroids = randomOrtVecFromData(data,k);
    % centroids = randomOrtVec();
    % 初始化标签
    labels = zeros(size(data, 2), 1);
    % disp(labels');
    % 迭代更新聚类中心和标签直到收敛
    converged = false;
    while ~converged
        cnt = cnt + 1;
        % 计算每个样本与聚类中心的距离
        distances = pdist2(data', centroids');
        
        % 分配每个样本到最近的聚类中心
        [~, new_labels] = min(distances, [], 2);


        % disp(new_labels);
        % 如果标签没有变化，则收敛
        if isequal(labels, new_labels)
            converged = true;
        else
            labels = new_labels;
            
            % 更新聚类中心为每个簇的均值
            for i = 1:k
                % disp(mean(data( :,labels == i), 2));
                centroids( :,i) = mean(data( :,labels == i), 2);
                clusters{i} = data( :,labels == i);
                be_labels{i} = be_label(labels == i);

            end
        end
    end
    % disp(clusters{k});
    % disp(be_labels{k});
    % disp("迭代次数：" + cnt + "次。");
end