function [d] = horizontal_dist(pos1,pos2)
% Return the horizontal distance of two object
% pos1 = [x1, y1];
% pos2 = [x2, y2];
d = sqrt( (pos1(1) - pos2(1))^2 + (pos1(2) - pos2(2))^2 );
end

