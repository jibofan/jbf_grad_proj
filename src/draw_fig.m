clear; clc; close all;

% ===== 参数设置 =====
mesh = 'ra';   % 网格名
e = 2;              % e 的值

% 横坐标
x = [25, 50, 100, 250, 500, 750, 1000, 1500, 2000];

% 三条折线的纵坐标
y1 = [];  % 首次k-way划分时间
y2 = [];  % 首次amd排序时间
y3 = [];  % 第二次重排序时间

% 数据目录：../tested/mesh_e
dataDir = sprintf('../tested/%s_%d', mesh, e);

% 顺次读取 runtime_N.txt 的最后一行后三个数
for i = 1:length(x)
    N = x(i);
    filename = fullfile(dataDir, sprintf('runtime_%d.txt', N));

    if ~isfile(filename)
        warning('文件不存在：%s', filename);
        y1(end+1) = NaN;
        y2(end+1) = NaN;
        y3(end+1) = NaN;
        continue;
    end

    % 读取所有行，取最后一行（忽略空行）
    lines = readlines(filename);
    lines = lines(strlength(strtrim(lines)) > 0);

    if isempty(lines)
        warning('文件为空：%s', filename);
        y1(end+1) = NaN;
        y2(end+1) = NaN;
        y3(end+1) = NaN;
        continue;
    end

    lastLine = char(lines(end));

    % 提取该行中的所有数字（支持整数/小数/科学计数法）
    nums = sscanf(lastLine, '%f');

    if numel(nums) < 3
        warning('最后一行数字不足3个：%s | 行内容: %s', filename, lastLine);
        y1(end+1) = NaN;
        y2(end+1) = NaN;
        y3(end+1) = NaN;
        continue;
    end

    % 取最后三个数，分别追加到 y1,y2,y3
    y1(end+1) = nums(end-2);
    y2(end+1) = nums(end-1);
    y3(end+1) = nums(end);
end

% 作图
fig = figure;
plot(x, y1, '-o', 'LineWidth', 1.5); hold on;
plot(x, y2, '-s', 'LineWidth', 1.5);
plot(x, y3, '-^', 'LineWidth', 1.5);
hold off;

% 标题、坐标轴、图例
title(sprintf('%s网格 e=%d 运行结果', mesh, e));
xlabel('叉数');
ylabel('时间/ms');
legend('首次k-way划分时间', '首次amd排序时间', '第二次重排序时间', 'Location', 'best');

% 美化
grid on;
xlim([min(x), max(x)]);
ylim('auto');

% ===== 保存图像到 mesh_e 文件夹，命名为 fig.png =====
outPath = fullfile(dataDir, 'fig.png');
exportgraphics(fig, outPath, 'Resolution', 300);   % 推荐
% 若你的MATLAB版本较老，可改用：saveas(fig, outPath);

fprintf('图像已保存到：%s\n', outPath);