function ref = vdhtable35()
   % Van de Hulst (1980), Multiple Light Scattering, Vol. 2, p. 435,
   % Table 35: plane-parallel slab of optical depth 2, single-scattering
   % albedo 0.9, Henyey-Greenstein g = 0.75, normal incidence, no surface
   % reflection. Hemispherical reflectance Rd is diffuse (the model has no
   % specular term); hemispherical transmittance Tt is total, diffuse plus
   % the direct exp(-2). The angular entries are the diffuse intensities
   % I(mu) at mu = cos(theta); the model's tallies per steradian equal
   % I(mu)*mu/pi, which fields R_sr and T_sr hold. The direct beam is not
   % part of the angular table.
   ref.ka = 10;
   ref.ks = 90;
   ref.g = 0.75;
   ref.Z = 0.02;
   ref.dz = 0.001;
   ref.Rd = 0.09739;
   ref.Tt = 0.66096;
   ref.Tdr = exp(-(ref.ka + ref.ks)*ref.Z);
   ref.mu = [0 0.1 0.3 0.5 0.7 0.9 1.0];
   ref.R = [0.10641 0.13304 0.13856 0.11956 0.09411 0.07132 0.06180];
   ref.T = [0.11811 0.16486 0.22473 0.27785 0.36889 0.72155 2.40270];
   ref.R_sr = ref.R.*ref.mu/pi;
   ref.T_sr = ref.T.*ref.mu/pi;
end
