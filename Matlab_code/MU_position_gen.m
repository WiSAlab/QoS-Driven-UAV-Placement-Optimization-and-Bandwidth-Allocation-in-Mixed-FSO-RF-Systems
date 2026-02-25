function [xy, numbPoints] = MU_position_gen(radius, numUser, center)
    % radius: in m

    xx0=center(1); yy0=center(2); %centre of disk
    
    %Point process parameters
    % areaTotal = radius^2; %area of disk
    % lambda = numUser/areaTotal; %intensity (ie mean density) of the Poisson process
    numbPoints = poissrnd(numUser);
    while numbPoints <= 0
        numbPoints = poissrnd(numUser); %Poisson number of points
    end
    
    theta = 2*pi*(rand(numbPoints,1)); %angular coordinates
    rho = radius*(rand(numbPoints,1)); %radial coordinates
    
    %Convert from polar to Cartesian coordinates
    [xx,yy] = pol2cart(theta,rho); %x/y coordinates of Poisson points
    %Shift centre of disk to (xx0,yy0)
    xx=xx+xx0;
    yy=yy+yy0;
    for i = 1:length(xx)
        xy(i,:) = [xx(i), yy(i)];
    end

end

