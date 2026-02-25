function [d] = euclid_dist(varargin)
% Return the euclidian distance of two object

switch nargin
    case 4
        % varargin = [pos1 pos2 h1 h2]
        % pos1 = [x1, y1];
        % pos2 = [x2, y2];
        d = sqrt( (varargin{1}(1) - varargin{2}(1))^2 + (varargin{1}(2) - varargin{2}(2))^2 ...
            + (varargin{3} - varargin{4})^2);
    case 3
        % varargin = [hor_dis h1 h2]
        d = sqrt( varargin{1}^2 + (varargin{2} - varargin{3})^2);
end
end


