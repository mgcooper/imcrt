function [ux,uy,uz] = chgdir2(ux,uy,uz,us,ps)
% if initial direction not straight up or straight down
if sqrt(1-uz*uz) < 1e-12
   ux = sqrt(1-us*us)*cos(ps);
   uy = sqrt(1-us*us)*sin(ps);
   uz = sign(uz)*us;
else
   % initial direction straight up or down (sin(theta)=0)
   ux = sqrt(1-us*us)/sqrt(1-uz*uz)*(ux*uz*cos(ps)-uy*sin(ps))+ux*us;
   uy = sqrt(1-us*us)/sqrt(1-uz*uz)*(uy*uz*cos(ps)+ux*sin(ps))+uy*us;
   uz = -sqrt(1-us*us)*sqrt(1-uz*uz)*cos(ps)+uz*us;
end