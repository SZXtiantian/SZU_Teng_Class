figFile1 = 'ran1000.fig';  % 替换成第一个.fig文件的路径
figFile2 = 'kmeans++1000.fig';  % 替换成第二个.fig文件的路径
figFile3 = 'kmeans++100.fig';  % 替换成第一个.fig文件的路径



fig1 = openfig(figFile1);
fig2 = openfig(figFile2);
fig3 = openfig(figFile3);
% 假设折线图对象在第一个.fig文件中
lines1 = findobj(fig1, 'Type', 'line');
% 假设折线图对象在第二个.fig文件中
lines2 = findobj(fig2, 'Type', 'line');
lines3 = findobj(fig3, 'Type', 'line');
% 遍历曲线对象并提取数据
for i = 1:length(lines1)
    % disp(lines1(i));
    xData = get(lines1(i), 'XData');
    yDataPre = get(lines1(i), 'YData');
    yDataKpp = get(lines2(i), 'YData');
    yDataRan = get(lines3(i), 'YData');
    % 这里可以使用 xData 和 yData 来处理数据
end

data_variance = var(yDataPre);
disp(['数据的方差为pre: ', num2str(data_variance)]);
data_variance1 = var(yDataKpp);
disp(['数据的方差为: ', num2str(data_variance1)]);
data_variance2 = var(yDataRan);
disp(['数据的方差为: ', num2str(data_variance2)]);