figFile1 = 'ran.fig';  % 替换成第一个.fig文件的路径
figFile2 = 'kmeans++100.fig';  % 替换成第二个.fig文件的路径
figFile3 = 'ranVec.fig';  % 替换成第一个.fig文件的路径



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
    disp(lines1(i));
    xData = get(lines1(i), 'XData');
    yDataPre = get(lines1(i), 'YData');
    yDataKpp = get(lines2(i), 'YData');
    yDataRan = get(lines3(i), 'YData');
    % 这里可以使用 xData 和 yData 来处理数据
end
mergedFig = figure;  % 创建一个新的图形窗口
ax = gca; % 获取当前的轴对象，或者创建一个新的轴对象
plot(ax, xData, yDataKpp,'o-','DisplayName', '相互最远向量'); % 绘制第一个折线图
hold on; % 保持图形窗口开启，以绘制第二个折线图
plot(ax, xData, yDataPre,'o-','DisplayName', '随机向量'); % 绘制第二个折线图
hold on; % 保持图形窗口开启，以绘制第二个折线图
plot(ax, xData, yDataRan,'o-','Color', [0.5, 0.2, 0.8],'DisplayName', '随机正交向量'); % 绘制第二个折线图
hold off; % 关闭保持模式，以防止进一步的绘图影响
legend('Location', 'northeast');