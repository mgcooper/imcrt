# imcrt impact report

Generated 07-Sep-2026 15:06:23 at commit 90cb624 by tests/verify/impactreport.m with N = 200000 packets per run and M = 10 runs.

src/ is unchanged since tag fix-S-bin-measures, so the post-fix column is the shipped kernel.

## 1. Impact of each fix

Versions in fix order: pre-fix commit b9d51d5, then tags fix-B-aliasing, fix-A-direct-tally, fix-C-fluence, fix-N-roulette, fix-S-bin-measures. Every version runs the same seed, so each table is one paired realization. A fix that changes trajectories or random-number use still leaves sampling error in that realization. A change is attributed to a fix when it exceeds 1e-6 relative. That separates a real change from summation-order rounding. r50 is quantized to one radial bin. The paper_radial case has the paper albedo and asymmetry (0.999, 0.9) at optical depth 10. Its lengths are scaled to the kernel's 2 cm radial grid. The paper's rod geometry is not part of src/mcrt.m.

### vdh_reflectance (ka 10, ks 90, g 0.75, Z 0.02, dz 0.001, seed 42)

| quantity | pre-fix | post-fix | change | fixes | note |
|---|---|---|---|---|---|
| Rdf | 0.0970871 | 0.097114 | +0.03% | A | hemispherical diffuse reflectance, a sum of packet weights |
| Tdf | 0.527148 | 0.527192 | +0.01% | A | hemispherical diffuse transmittance, a sum of packet weights |
| Tdr | 0.13408 | 0.134035 | -0.03% | A | direct transmittance; A moves scattered grazing exits out of it |
| Rdr | 2.68779e-05 | 0 | -100.00% | A | direct reflectance; zero for the vertical source after A |
| Adf | 0.155062 | 0.155062 | -0.00% | none | total diffuse absorbed fraction, a sum of packet weights |
| Rdf_a(1) | 0.0150691 | 0.0178577 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Tdf_a(1) | 0.640231 | 0.758708 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Rdf_r(1) | 273.22 | 193.682 | -29.11% | B, S | first radial bin per area; B moves packets, S corrects the area |
| Tdf_r(1) | 5610.84 | 5325.41 | -5.09% | B, S | first radial bin per area; B moves packets, S corrects the area |
| r50 R | 0.0195043 | 0.0215039 | +10.25% | B | radius enclosing half the reflected diffuse weight; B |
| r50 T | 0.00950877 | 0.00950877 | +0.00% | none | radius enclosing half the transmitted diffuse weight; B |
| phi_rz(1,1) | 15543.1 | 321653 | +1969.42% | B, C, S | on-axis surface fluence density; C confines the pencil, S sets dV |
| phi_z(1) | 1.25083 | 1.25083 | +0.00% | none | surface plane-integrated fluence; unchanged by every fix |

### absorbing (ka 90, ks 10, g 0.5, Z 0.05, dz 0.005, seed 43)

| quantity | pre-fix | post-fix | change | fixes | note |
|---|---|---|---|---|---|
| Rdf | 0.00522792 | 0.00523302 | +0.10% | A | hemispherical diffuse reflectance, a sum of packet weights |
| Tdf | 0.00160243 | 0.00160243 | -0.00% | none | hemispherical diffuse transmittance, a sum of packet weights |
| Tdr | 0.00692 | 0.00692 | +0.00% | none | direct transmittance; A moves scattered grazing exits out of it |
| Rdr | 5.106e-06 | 0 | -100.00% | A | direct reflectance; zero for the vertical source after A |
| Adf | 0.0924726 | 0.0924726 | -0.00% | none | total diffuse absorbed fraction, a sum of packet weights |
| Rdf_a(1) | 0.000673046 | 0.000797595 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Tdf_a(1) | 0.00119062 | 0.00141095 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Rdf_r(1) | 35.325 | 39.9597 | +13.12% | B, S | first radial bin per area; B moves packets, S corrects the area |
| Tdf_r(1) | 6.25468 | 6.98483 | +11.67% | B, S | first radial bin per area; B moves packets, S corrects the area |
| r50 R | 0.00950877 | 0.0105079 | +10.51% | B | radius enclosing half the reflected diffuse weight; B |
| r50 T | 0.0135062 | 0.0135062 | +0.00% | none | radius enclosing half the transmitted diffuse weight; B |
| phi_rz(1,1) | 2943.84 | 253240 | +8502.39% | B, C, S | on-axis surface fluence density; C confines the pencil, S sets dV |
| phi_z(1) | 0.819286 | 0.819286 | -0.00% | none | surface plane-integrated fluence; unchanged by every fix |

### mini_fluence (ka 1, ks 99, g 0.9, Z 0.1, dz 0.005, seed 44)

| quantity | pre-fix | post-fix | change | fixes | note |
|---|---|---|---|---|---|
| Rdf | 0.247452 | 0.247512 | +0.02% | A | hemispherical diffuse reflectance, a sum of packet weights |
| Tdf | 0.589189 | 0.589234 | +0.01% | A | hemispherical diffuse transmittance, a sum of packet weights |
| Tdr | 7.49604e-05 | 3e-05 | -59.98% | A | direct transmittance; A moves scattered grazing exits out of it |
| Rdr | 5.97668e-05 | 0 | -100.00% | A | direct reflectance; zero for the vertical source after A |
| Adf | 0.153225 | 0.153225 | +0.00% | none | total diffuse absorbed fraction, a sum of packet weights |
| Rdf_a(1) | 0.0465004 | 0.0551054 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Tdf_a(1) | 0.301283 | 0.357036 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Rdf_r(1) | 151.525 | 56.837 | -62.49% | B, S | first radial bin per area; B moves packets, S corrects the area |
| Tdf_r(1) | 682.94 | 110.014 | -83.89% | B, S | first radial bin per area; B moves packets, S corrects the area |
| r50 R | 0.0705012 | 0.083501 | +18.44% | B | radius enclosing half the reflected diffuse weight; B |
| r50 T | 0.0455018 | 0.0515016 | +13.19% | B | radius enclosing half the transmitted diffuse weight; B |
| phi_rz(1,1) | 48516.8 | 305321 | +529.31% | B, C, S | on-axis surface fluence density; C confines the pencil, S sets dV |
| phi_z(1) | 1.56451 | 1.56451 | -0.00% | none | surface plane-integrated fluence; unchanged by every fix |

### paper_radial (ka 0.01, ks 9.99, g 0.9, Z 1, dz 0.05, seed 45)

| quantity | pre-fix | post-fix | change | fixes | note |
|---|---|---|---|---|---|
| Rdf | 0.300546 | 0.30059 | +0.01% | A | hemispherical diffuse reflectance, a sum of packet weights |
| Tdf | 0.680887 | 0.680946 | +0.01% | A | hemispherical diffuse transmittance, a sum of packet weights |
| Tdr | 0.000133909 | 7.5e-05 | -43.99% | A | direct transmittance; A moves scattered grazing exits out of it |
| Rdr | 4.43516e-05 | 0 | -100.00% | A | direct reflectance; zero for the vertical source after A |
| Adf | 0.0173887 | 0.0173887 | -0.00% | none | total diffuse absorbed fraction, a sum of packet weights |
| Rdf_a(1) | 0.0685096 | 0.0811875 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Tdf_a(1) | 0.314354 | 0.372526 | +18.51% | S | first angular bin per steradian; S corrects its solid angle |
| Rdf_r(1) | 5.35006 | 1.58996 | -70.28% | B, S | first radial bin per area; B moves packets, S corrects the area |
| Tdf_r(1) | 13.2831 | 0 | -100.00% | B | first radial bin per area; B moves packets, S corrects the area |
| r50 R | 0.5515 | 0.7695 | +39.53% | B | radius enclosing half the reflected diffuse weight; B |
| r50 T | 0.3985 | 0.5065 | +27.10% | B | radius enclosing half the transmitted diffuse weight; B |
| phi_rz(1,1) | 15972.2 | 269996 | +1590.41% | B, C, S | on-axis surface fluence density; C confines the pencil, S sets dV |
| phi_z(1) | 1.71025 | 1.71025 | -0.00% | none | surface plane-integrated fluence; unchanged by every fix |

## 2. Paper-kernel retrospective

The retrospective runs the van de Hulst case M = 10 times at N = 200000, seeds 101 to 110. A scratch copy of the frozen 2021 verification script, extracted with its helpers from commit b9d51d5, and the fixed kernel each run every seed. The copy carries the production kernel's direct clause (ns==0 || ia==0) in place of the script's grazing clause. It therefore scores exits the way the published runs did. It keeps the out-of-place chgdir (no defect B) and the shifted bin measures (defect S). It also keeps the old roulette (defect N) at the script's wmin = 1e-4, which the fixed kernel shares. The production kernel used wmin = 1e-5. The text under the table gives the measured effect of that threshold and of the production kernel's random-number use. se is the run spread over sqrt(M). The two kernels share seeds and the same random-number order, so their runs are paired. diff/se is the mean paired difference (fixed minus 2021) in units of its own standard error; identical runs give none.

| quantity | 2021 mean | 2021 se | fixed mean | fixed se | paired diff | diff/se |
|---|---|---|---|---|---|---|
| Rdf | 0.0974559 | 0.00017 | 0.0974559 | 0.00017 | 0 | none |
| Tdf | 0.525521 | 0.00029 | 0.525521 | 0.00029 | 0 | none |
| Tdr | 0.135455 | 0.00029 | 0.135455 | 0.00029 | 0 | none |
| Rdr | 0 | 0 | 0 | 0 | 0 | none |
| Adf | 0.155114 | 0.00012 | 0.155114 | 0.00012 | 0 | none |
| Rdf_a(1) | 0.0180871 | 0.00068 | 0.0214342 | 0.00081 | 0.00335 | 26.43 |
| Tdf_a(1) | 0.632731 | 0.0063 | 0.74982 | 0.0075 | 0.117 | 100.62 |
| Rdf_r(1) | 166.388 | 4.6 | 197.201 | 5.5 | 30.8 | 35.84 |
| Tdf_r(1) | 4398.49 | 19 | 5213.02 | 23 | 815 | 227.92 |
| r50 R | 0.0219038 | 0.00016 | 0.0219038 | 0.00016 | 0 | none |
| r50 T | 0.00950877 | 0 | 0.00950877 | 0 | 0 | none |
| phi_rz(1,1) | 15203.8 | 1.7e+02 | 321578 | 6.4e+02 | 3.06e+05 | 545.99 |
| phi_z(1) | 1.25164 | 0.0024 | 1.25164 | 0.0024 | 0 | none |

Reference values: Rd 0.09739, Tt 0.66096, Tdr 0.13534.

The copy's complete tallies at wmin = 1e-5, the production threshold, and at wmin = 1e-4 are identical run for run. A packet that reaches wmin = 1e-4 draws a roulette number at that threshold only, so identical tallies mean that no packet reached it. Roulette is unbiased, so the threshold changes the variance only. The hemispherical shifts in standard errors of the paired difference:

- Rdf: none
- Tdf: none
- Tdr: none
- Rdr: none
- Adf: none

A second copy also takes the production kernel's random-number use: the path length from log(1-rand) and the Henyey-Greenstein terms hg4 = 1-g, hg5 = 2*g. That relabels the random stream and leaves every sampled distribution unchanged. This copy runs seeds 111 to 120, independent of the first copy's. Its hemispherical means at wmin = 1e-5, in standard errors of the unpaired difference from the first copy at wmin = 1e-5, shift by:

- Rdf: -0.39
- Tdf: 0.50
- Tdr: -0.82
- Rdr: none
- Adf: 0.94

### How the published numbers relate to the corrected model

The 2021 kernel used the out-of-place chgdir, so defect B never touched it. Its radial spread (r50) matches the fixed kernel. Its production loop counted only unscattered packets as direct: ia is 0 only for an exit exactly along the normal. Defect A therefore applies only to its verification script, and section 1 measures that clause as the A column. Its grid builder is the shipped one, so defect S applies. In any 2021 output that reports them, the first angular and radial bins per steradian or per area are 15.6% below the corrected value. The corrected value is 18.5% higher. The second bins are 3.3% above it. Defect S changes the per-bin densities only, not the hemispherical sums. Its fluence density adds the per-depth direct absorption to every radial column (defect C). The fix confines that pencil to bin 1 as a volume density. Defect N drops a boosted packet that still sits below wmin. That needs an albedo below 1/wrr = 0.1, so it never fires at the paper albedos. The hemispherical reflectance, transmittance, and diffuse absorption of the fixed kernel and the paired copy are identical run for run, as the paired differences show. The production-stream copy differs from that copy by the listed standard-error ratios; a relabeled random stream predicts ratios of order one.
