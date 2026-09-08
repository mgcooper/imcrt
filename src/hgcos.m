function us = hgcos(g)
   % Cosine of the polar scattering angle drawn from the Henyey-Greenstein
   % phase function with asymmetry g (Wang et al. 1995, Eq. 3.28, inverse
   % transform of one uniform draw). g = 0 is isotropic, where the general
   % form divides by zero. The mean of us is g.
   if g == 0
      us = 1-2*rand;
   else
      us = 1/(2*g)*((1+g^2)-((1-g^2)/((1+g)+(-2*g)*rand))^2);
   end
end
