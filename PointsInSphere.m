function [ data ] = PointsInSphere( n, a, b, x, y, z )

    if nargin==3
       x = 0;
       y = 0;
       z = 0;
    end
    
    r1 = (rand(n,1)*(b^3-a^3)+a^3).^(1/3);
    cphi = -1 + 2*rand(n,1); 
    sphi = sqrt(1-cphi.^2); 
    th1 = 2*pi*rand(n,1);
    % Convert to cart.
    z = r1.*cphi + Z; 
    r1 = r1.*sphi; 
    x = r1.*sin(th1) + X; 
    y = r1.*cos(th1) + Y;



end

