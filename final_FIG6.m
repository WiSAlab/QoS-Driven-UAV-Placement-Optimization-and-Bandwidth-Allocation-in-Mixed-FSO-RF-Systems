% Dữ liệu
x = 1:6;  % số UAV

% Rate FSO 3 GHz
opt_3g = [413 1121 1664 2016 2320 2536];       % Xanh trơn
tla_3g = [646 1224 1919 2159 2435 2747]; % Đỏ trơn

% Rate FSO 2 GHz
opt_2g = [400 857 1163 1513 1642 1836];       % Xanh sọc
tla_2g = [558 1096 1484 1673 1731 1962];      % Đỏ sọc

data = [opt_3g(:), tla_3g(:), opt_2g(:), tla_2g(:)];

% Vẽ biểu đồ
figure;
b = bar(x, data, 'grouped');
hold on;

% Màu sắc
b(1).FaceColor = 'b';  % Xanh trơn
b(2).FaceColor = 'r';  % Đỏ trơn
b(3).FaceColor = [0.4 0.6 0.2];  % Xanh sọc
b(4).FaceColor = [1 0.6 0.2];  % Đỏ sọc

% Thêm số liệu trên đầu mỗi cột, căn giữa và không tràn
offset = 0.015;
for i = 1:4
    x_vals = b(i).XEndPoints;
    y_vals = b(i).YEndPoints;
    labels = string(round(b(i).YData, 1)) ; % Thêm % như trong hình
    text(x_vals - offset, y_vals + 25, labels, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 13, 'FontWeight', 'bold');
end

% Thêm marker trên thân 3 trong 4 cột với cùng số lượng (3 cái mỗi cột)
marker_list = {'o', '+', '*'};  % Marker khác nhau
marker_size = 6;
y_ratios = [0.3, 0.5, 0.7];    % Vị trí giống nhau cho mọi cột

for i = 1:length(x)
%     Cột 2: TLA (3GHz)
    x2 = b(2).XEndPoints(i); y2 = b(2).YEndPoints(i);
    for r = y_ratios
        plot(x2, y2 * r, marker_list{1}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

%     Cột 3: OPT-UAV (2GHz)
    x3 = b(3).XEndPoints(i); y3 = b(3).YEndPoints(i);
    for r = y_ratios
        plot(x3, y3 * r, marker_list{2}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

%     Cột 4: TLA (2GHz)
    x4 = b(4).XEndPoints(i); y4 = b(4).YEndPoints(i);
    for r = y_ratios
        plot(x4, y4 * r, marker_list{3}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
end

% Nhãn trục, giới hạn
xlabel('Number of UAVs in emergency zone', ...
    'FontSize', 16, 'FontWeight', 'bold');
ylabel('Total data rate (Mbps)', ...
    'FontSize', 16, 'FontWeight', 'bold');
ylim([0 3500]);

% Legend
legend({'OPT-UAV (3Gbps)', 'TLA (3Gbps)', 'OPT-UAV (2Gbps)', 'TLA (2Gbps)'}, ...
    'FontSize', 16, 'Location', 'northoutside', ...
    'Orientation', 'horizontal');

grid on;



