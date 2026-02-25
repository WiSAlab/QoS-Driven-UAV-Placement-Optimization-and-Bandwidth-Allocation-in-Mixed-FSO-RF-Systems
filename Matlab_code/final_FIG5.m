% MU_density_list = 25:30;
% results1 = [56 52.5 50 48 46.4 44];     % OPT-UAV 3GHz (xanh trơn)
% results2 = [43 42 40.5 40.5 38 35];     % TLA 3GHz (đỏ trơn)
% results1_2G = [40 38 35 32 28 25];      % OPT-UAV 2GHz (xanh sọc)
% results2_2G = [33 30 28 26 24 21];      % TLA 2GHz (đỏ sọc)
% 
% data = [results1(:), results2(:), results1_2G(:), results2_2G(:)];
% 
% figure;
% b = bar(MU_density_list, data, 'grouped');
% hold on;
% 
% % Thiết lập màu cơ bản
% b(1).FaceColor = 'b';  % Xanh trơn
% b(2).FaceColor = 'r';  % Đỏ trơn
% b(3).FaceColor = [0.4 0.6 0.2];  % Xanh sọc
% b(4).FaceColor = [1 0.6 0.2];  % Đỏ sọc
% 
% % % Áp dụng hatchfill cho các cột 2GHz (sọc)
% % hatchfill2(b(3), 'single', 45, 5, 'w');  % Xanh sọc trắng
% % set(b(2), 'FaceColor', 'b');
% % hatchfill2(b(4), 'single', 45, 5, 'w');  % Đỏ sọc trắng
% % set(b(4), 'FaceColor', 'r');
% 
% % Gán nhãn trục
% xlabel('Average number of MUs per 200m x 200m location', ...
%     'FontSize', 18, 'FontWeight', 'bold');
% ylabel('Fraction of served MUs (%)', ...
%     'FontSize', 18, 'FontWeight', 'bold');
% % title('Fraction of Served MUs vs. MU Density (3GHz vs. 2GHz)', ...
% %     'FontSize', 18, 'FontWeight', 'bold');
% ylim([0 60]);
% grid on;
% 
% % Thêm nhãn dữ liệu trên đầu mỗi cột (căn giữa chính xác)
% x_offset = 0.05;  % Dịch nhẹ sang phải để căn giữa hoàn hảo
% for i = 1:length(b)
%     x = b(i).XEndPoints;
%     y = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1)) + "%";
%     text(x + x_offset, y + 2.5, labels, ...
%         'HorizontalAlignment', 'center', ...
%         'FontSize', 13, 'FontWeight', 'bold');
% end
% 
% % Thêm chú thích
% legend({'OPT-UAV (3GHz)', 'TLA (3GHz)', 'OPT-UAV (2GHz)', 'TLA (2GHz)'}, ...
%     'FontName', 'Times New Roman', ...
%     'FontSize', 14, ...
%     'Location', 'northoutside', ...
%     'Orientation', 'horizontal');
MU_density_list = 25:30;
results1 = [56 52.5 50 48 46.4 44];     % OPT-UAV 3GHz (xanh trơn)
results2 = [43 42 40.5 39 38 35];     % TLA 3GHz (đỏ trơn)
results1_2G = [40 38 35 32 28 25];      % OPT-UAV 2GHz (xanh sọc)
results2_2G = [33 30 28 26 24 21];      % TLA 2GHz (đỏ sọc)

data = [results1(:), results2(:), results1_2G(:), results2_2G(:)];

figure;
b = bar(MU_density_list, data, 'grouped');
hold on;

% Thiết lập màu cơ bản
b(1).FaceColor = 'b';                  % Xanh trơn
b(2).FaceColor = 'r';                  % Đỏ trơn
b(3).FaceColor = [0.4 0.6 0.2];        % Xanh sọc (olive green)
b(4).FaceColor = [1 0.6 0.2];          % Đỏ sọc (orange)

% Gán nhãn trục
xlabel('Average number of MUs per 200m x 200m location', ...
    'FontSize', 16, 'FontWeight', 'bold');
ylabel('Fraction of served MUs (%)', ...
    'FontSize', 16, 'FontWeight', 'bold');
ylim([0 60]);
grid on;

% Thêm nhãn dữ liệu trên đầu mỗi cột (căn giữa chính xác)
x_offset = 0.02;  % Dịch nhẹ sang phải để căn giữa hoàn hảo
for i = 1:length(b)
    x = b(i).XEndPoints;
    y = b(i).YEndPoints;
    labels = string(round(b(i).YData, 1)) + "%";
    text(x + x_offset, y + 1.5, labels, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 13, 'FontWeight', 'bold');
end

% Thêm nhiều marker trên thân 3 trong 4 cột với cùng số lượng (3 cái mỗi cột)
marker_list = {'o', '+', '*'};                 % Marker khác nhau
marker_size = 6;
y_ratios = [0.3, 0.5, 0.7];                    % Vị trí giống nhau cho mọi cột

for i = 1:length(MU_density_list)
    % Cột 2: TLA 3GHz
    x2 = b(2).XEndPoints(i); y2 = b(2).YEndPoints(i);
    for r = y_ratios
        plot(x2, y2 * r, marker_list{1}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % Cột 3: OPT-UAV 2GHz
    x3 = b(3).XEndPoints(i); y3 = b(3).YEndPoints(i);
    for r = y_ratios
        plot(x3, y3 * r, marker_list{2}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % Cột 4: TLA 2GHz
    x4 = b(4).XEndPoints(i); y4 = b(4).YEndPoints(i);
    for r = y_ratios
        plot(x4, y4 * r, marker_list{3}, 'MarkerSize', marker_size, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
end
% Thêm chú thích
legend({'OPT-UAV (3Gbps)', 'TLA (3Gbps)', 'OPT-UAV (2Gbps)', 'TLA (2Gbps)'}, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal');
