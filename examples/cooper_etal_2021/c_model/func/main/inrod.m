function tf = inrod(x,y,z,cgeom)
% cgeom  = [xc0,yc0,zc0,rc,hc]; [xcenter ycenter zcenter radius height]
tf=(x>cgeom(1)&&x<(cgeom(1)+cgeom(5))&&((y-cgeom(2))^2+(z-cgeom(3))^2)<(cgeom(4)^2));