Q4 result: material-ALE xi=r/R(t), attachment4 properties, radius attachment2 piecewise linear, Kirchhoff moisture flux, BDF sparse.
N160 drying 182969.317394s = 50.8248104h
N320 182967.668476s = 50.8243524h
N640 rtol1e-9 182967.268554s = 50.8242413h; use 50.8242 h.
Files: private/code/Q4/q4_shrink.py; private/results/Q4/q4_ALE_N640_tol1e-09.json; Q4_convergence.csv; Q4_ALE_summary.md.
Table rows in JSON, columns physical 0,0.5,1,1.5,2 cm + surface. NaN means physical position outside shrunken radius. At endpoint R=1.2cm; surface C=0.0524963, center C=.15. Integer end 182968s terminal max C=.149999622.
ALE equations: rho cp T_t = R^-2 xi^-1 (xi k T_xi)_xi; C_t=R^-2 xi^-1(xi D C_xi)_xi. Robin boundary physical flux scales 1/R in xi form; implemented Rt* h(Tinf-Ts)/Rt^2 equivalent.
