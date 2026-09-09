function us = hgcos(g)
   %HGCOS Cosine of polar scattering angle from Henyey-Greenstein phase function
   %
   %  us = hgcos(g) returns the cosine of the polar scattering angle drawn from
   %  the Henyey-Greenstein phase function with asymmetry parameter g (Wang et
   %  al. 1995, Eq. 3.28, inverse transform of one uniform draw).
   %
   %  Pass g = 0 for isotropic scattering. The mean of us is g.
   %
   % See also: mcrt

   if g == 0
      us = 1-2*rand;
   else
      us = 1/(2*g)*((1+g^2)-((1-g^2)/((1+g)+(-2*g)*rand))^2);
   end
end
