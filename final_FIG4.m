clc;
clear all;

% % Dữ liệu
rate_3GHz = [2016.40, 2159.00];
rate_2GHz = [1556.41, 1673.06];
percent_3GHz = [67, 62];
percent_2GHz = [64, 60];

labels = {'OPT-UAV', 'TLA'};

% Tăng khoảng cách giữa 2 nhóm
x = [1, 4];               % cách xa hơn nữa giữa 2 nhóm
bw = 0.18;                 % độ rộng mỗi cột
x_pos = [-1.5, -0.5, 0.5, 1.5] * bw * 3.4;  % giãn đều các cột trong nhóm

% figure('Color', 'w'); hold on;
figure('Color', 'w', 'Position', [100 100 580 580]); 
hold on;

% === Trục trái: Data rate ===
yyaxis left
b1 = bar(x + x_pos(2), rate_3GHz, bw, 'FaceColor', 'r', 'EdgeColor', 'r'); % Rate 3GHz
b2 = bar(x + x_pos(4), rate_2GHz, bw, 'FaceColor', [0.4 0.6 0.2], 'EdgeColor', [0.4 0.6 0.2]); % Rate 2GHz

ylabel('Total data rate (Mbps)', 'FontWeight', 'bold', ...
    'Color', 'k', 'FontSize', 16, 'FontName', 'Times New Roman');
ylim([1300 2300]);
set(gca, 'ycolor', 'k');

% === Trục phải: Percentage ===
yyaxis right
b3 = bar(x + x_pos(1), percent_3GHz, bw, 'FaceColor', 'b', 'EdgeColor', 'b'); % % 3GHz
b4 = bar(x + x_pos(3), percent_2GHz, bw, 'FaceColor', [1 0.6 0.2], 'EdgeColor', [1 0.6 0.2]); % % 2GHz

ylabel('Percentage of served MUs (%)', 'FontWeight', 'bold', ...
    'Color', 'k', 'FontSize', 16, 'FontName', 'Times New Roman');
ylim([40 100]);
set(gca, 'ycolor', 'k');

% Trục X
set(gca, 'xtick', x, 'xticklabel', labels, 'FontSize', 16, ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold');

% Hiển thị số liệu trên đầu cột
yyaxis left
for i = 1:2
    text(x(i)+x_pos(2), rate_3GHz(i)+50, sprintf('%.2f', rate_3GHz(i)), ...
        'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
    text(x(i)+x_pos(4), rate_2GHz(i)+50, sprintf('%.2f', rate_2GHz(i)), ...
        'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
end

yyaxis right
for i = 1:2
    text(x(i)+x_pos(1), percent_3GHz(i)+2, sprintf('%.f%%', percent_3GHz(i)), ...
        'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
    text(x(i)+x_pos(3), percent_2GHz(i)+2, sprintf('%.f%%', percent_2GHz(i)), ...
        'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
end

% Thêm marker vào 3 trong 4 cột ===
marker_list = {'o', '+', '*'};  % Mỗi loại marker cho mỗi cột
marker_size = 9;
y_ratios = [0.3, 0.7, 0.9];

for i = 1:length(x)
    % --- Cột Rate 3GHz (b1) - trục trái ---
    yyaxis left
    x1 = b1.XEndPoints(i); y1 = b1.YEndPoints(i);
    for r = y_ratios
        plot(x1, y1 * r, marker_list{1}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % --- Cột Rate 2GHz (b2) - trục trái ---
    x2 = b2.XEndPoints(i); y2 = b2.YEndPoints(i);
    for r = y_ratios
        plot(x2, y2 * r, marker_list{2}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % --- Cột % 2GHz (b4) - trục phải ---
    yyaxis right
    x4 = b4.XEndPoints(i); y4 = b4.YEndPoints(i);
    for r = y_ratios
        plot(x4, y4 * r, marker_list{3}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
end



% Legend đúng thứ tự dữ liệu
legend([b3, b1, b4, b2], ...
    {'Percentage of served MUs(3Gbps)', 'Total data rate(3Gbps)', 'Percentage of served MUs(2Gbps)', 'Total data rate(2Gbps)'}, ...
    'Location','northoutside', ...
    'Orientation','vertical', ...
    'FontSize',14, ...
    'FontName','Times New Roman', ...
    'FontWeight','bold', ...
    'EdgeColor','k', ...
    'Box','on');

% Thẩm mỹ
box on; grid on;
set(gca, 'LineWidth', 1.0);


