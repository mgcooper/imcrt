%==========================================================================
function [xscat,yscat,zscat,rayl,tf]=rodintersect(x,y,z,ux,uy,uz,rayl,rod)
%==========================================================================
% 1. check if the line segment is entirely outside of the cylinder's domain
%       if yes, we are done, the line does not intersect the cylinder
% 2. if no, check if the ray pierces into the cylinder
% 3. if no, check if the y-z line projection intersects the y-z circle
%       if no, we are done, the line does not intersect the cylinder
% 4. if yes, minimize the distance from the line to the circle center
%       if this distance is > than the radius, the line segement does not
%       intersect the cylinder
%       otherwise, there is at least one intersection
% 5. choose the y-z solution from (2) that is nearest to the original point
%==========================================================================
% Inputs:   x: x-position of ray relative to the origin
%           y: y-position of ray relative to the origin
%           z: z-position of ray relative to the origin
%           dx: distance from prior x-position to x along the x-axis
%           dy: distance from prior y-position to y along the x-axis
%           dz: distance from prior z-position to z along the x-axis
%           x: x-position of cylinder base center relative to the origin
%           y: y-position of cylinder base center relative to the origin
%           z: z-position of cylinder base center relative to the origin
%           r: radius cylinder base center relative to the origin
%               cx/y/zlims: limits of the cylinder along each axis

% Outputs:  x/y/zscat: position of ray intersection with cylinder nearest
%                       the original x-position (i.e. the scattering point)
%           rayl: reflected (scattered) ray length
%           tf: true if the ray intersected the cylinder
%==========================================================================

% recover the start point x0,y0,z0 from the direction cosines and ray length
dx = ux*rayl; dy = uy*rayl; dz = uz*rayl;
x0 = x-dx; y0 = y-dy; z0 = z-dz;

% recover the detector rod geometry:
% for reference: rod = [xc0,yc0,zc0,rc,hc]; [xcenter ycenter zcenter radius height]
xc0 = rod(1);
yc0 = rod(2);
zc0 = rod(3);
rc  = rod(4);
hc  = rod(5);
cx  = [xc0,xc0+hc];     % x-limits of rod 
cy  = [yc0-rc,yc0+rc];  % y-limits of rod
cz  = [zc0-rc,zc0+rc];  % z-limits of rod

% 1. if any of these are true, then the ray does not intersect the cylinder
%==========================================================================
    if (    (x>=cx(2) && x0>=cx(2)) || (x<=cx(1) && x0<=cx(1))     || ...
            (y>=cy(2) && y0>=cy(2)) || (y<=cy(1) && y0<=cy(1))     || ...
            (z>=cz(2) && z0>=cz(2)) || (z<=cz(1) && z0<=cz(1))     )
        xscat = x; yscat = y; zscat = z; tf = false;
        return;
    end

% 2. check if the ray intersects the end of the rod. The start point can
% never be inside the cylinder, so it is sufficient to check if the ray
% starts outside and ends inside, and then to check if the end point is
% within the circle projection in the y-z. The case where the ray pokes
% through the end cap and into the cylinder surface is detected in step 4.
%==========================================================================
    if (x0<=cx(1) && x>cx(1)) || (x>cx(2) && x0<cx(2))
        if ((y-yc0)^2 + (z-zc0)^2) < rc^2
            t       = (xc0-x0)/dx;
            xscat   = xc0;
            yscat   = y0+dy*t;
            zscat   = z0+dz*t;
            rayl    = sqrt((x-xscat)^2+(y-yscat)^2+(z-zscat)^2);
            tf      = true;
            return;
        end
    end
    
% 3. find solutions to the line-circle equation in the y-z plane
%==========================================================================
    ystar   = y0-yc0;
    zstar   = z0-zc0;
    a       = dz*dz + dy*dy;
    b       = 2*zstar*dz + 2*ystar*dy;
    c       = zstar*zstar + ystar*ystar - rc*rc;
    
    if (b^2-4*a*c) < 0 % there are no intersections
        xscat = x; yscat = y; zscat = z; tf = false;
        return;
    else 
        s = (-b+[sqrt(b^2-4*a*c) -sqrt(b^2-4*a*c)])./(2*a);
    end

% 4. find intersections with the cylinder if they exist
%==========================================================================
% if s is outside the range [0,1] there are no intersections
    s = s(s>=0 & s<=1);
    if isempty(s)
        xscat = x; yscat = y; zscat = z; tf = false;
        return;
    else 
% check the x-solutions
        xsoln = (x0-xc0)+dx*s; 
% I commented this out, I think this is unneccesary
%         if (xsoln(1) <= cxlims(1) && xsoln(2) <= cxlims(1) || ...
%             xsoln(1) >= cxlims(2) && xsoln(2) >= cxlims(2) )
% % both x sol's are outside of the cylinder domain, there is no intersection
%             xscat   = x; yscat = y; zscat = z; tf = false;
%             return;
%         else
% choose the solution that is nearest to the original location                        
            ysoln   =   y0+dy*s;
            zsoln   =   z0+dz*s;
            d       =   sqrt((x0-xsoln).^2+(y0-ysoln).^2+(z0-zsoln).^2);
            xscat   =   xsoln(d==(min(d)));
            yscat   =   ysoln(d==(min(d)));
            zscat   =   zsoln(d==(min(d)));
            rayl    =   sqrt((x-xscat)^2+(y-yscat)^2+(z-zscat)^2);
            tf      =   true;
%         end
    end    
end


