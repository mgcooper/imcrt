function M = RotationMatrix(a, b)
   %ROTATIONMATRIX Rotation matrix that aligns vector a with vector b.
   %
   %  M = RotationMatrix(a, b) returns the 3x3 rotation matrix M such that
   %  M*a/norm(a) equals b/norm(b). M is the minimal rotation, about the
   %  axis cross(a, b), in Rodrigues' axis-angle form. a and b are
   %  3-element vectors in any orientation; neither may be zero.
   %
   %  When a and b are collinear within 1e-12, M depends on their sense.
   %  The parallel case gives the identity. The antiparallel case gives a
   %  half-turn about an axis perpendicular to a. The azimuth of that axis
   %  is arbitrary, so a caller must compare only rotation-invariant
   %  quantities in that case.
   %
   %  Near antiparallel the axis direction is conditioned as eps divided by
   %  norm(cross(a, b)), but M*a still equals b to round-off.
   %
   %  This file is the independent oracle for the direction-cosine update
   %  in src/mcrt.m and src/chgdir.m. It is not on the production path.
   %
   %  See also cross, dot, chgdir.

   % The formula needs unit column vectors. A zero vector has no direction,
   % so reject it before dividing by its norm.
   assert(isequal(size(a(:)), [3 1]) && norm(a) > 0)
   assert(isequal(size(b(:)), [3 1]) && norm(b) > 0)
   a = a(:) / norm(a);
   b = b(:) / norm(b);

   % The rotation angle has sine norm(v) and cosine c. Below 1e-12 the cross
   % product is rounding noise, so its direction is unusable as an axis.
   v = cross(a, b);
   c = dot(a, b);
   s = norm(v);
   if s < 1e-12
      if c > 0
         M = eye(3);
      else
         M = halfturn(a);
      end
      return
   end

   % Rounding leaves cross(a, b) off perpendicular to a by about eps, and
   % that error is amplified by 1/s. Re-orthogonalize the unit axis against
   % a so that M*a lands on b to round-off even when s is tiny.
   k = v / s;
   k = k - a * (a' * k);
   k = k / norm(k);

   % Rodrigues' axis-angle form: I + sin(t) K + (1 - cos(t)) K^2, with K the
   % cross-product matrix of the unit axis. No factor divides by 1 + c.
   K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
   M = eye(3) + s*K + (1 - c)*(K*K);
end

function H = halfturn(a)
   % Half-turn about a unit axis p perpendicular to a: H = 2pp' - I. Build
   % p from the coordinate axis least aligned with a, so the cross product
   % is never degenerate. The azimuth of p is arbitrary.
   [~, i] = min(abs(a));
   p = cross(a, double((1:3)' == i));
   p = p / norm(p);
   H = 2*(p*p') - eye(3);
end
