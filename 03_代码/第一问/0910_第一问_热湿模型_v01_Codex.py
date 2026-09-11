"""A Q1: conservative axisymmetric FVM, nonlinear trapezoidal time stepping.

2026-09-10, Codex. Run with Python having numpy/scipy; input is audited CSV.
The independent 1D baseline uses exactly the same radial discretization.
No phase-change heat is added to the approved effective-property model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import LinearOperator, cg, splu

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO.parent
INPUT = WORKSPACE / 'CUMCM2026Problems/A题/附件/附件1.xlsx'
BOUNDARY_CSV = REPO / '02_数据/清洗数据/0910_第一问_烘房边界_v01_Codex.csv'
OUT = REPO / '04_结果/第一问'
R, HALF_LENGTH = .02, .125
RHO, CP, K, H, HM = 820., 2600., .36, 25., 8e-7
ALPHA = K / (RHO * CP)
TABLE_TIMES = np.array([100, 300, 600, 900, 1200, 1500, 1800])
OUTPUT_R = np.arange(21) * .001


def environment():
    a = np.loadtxt(BOUNDARY_CSV, delimiter=',', skiprows=1, encoding='utf-8-sig')
    assert a.shape == (241, 3) and np.isfinite(a).all()
    assert np.array_equal(a[:, 0], np.arange(0, 14401, 60))
    return a


def diffusivity(c):
    if np.any(c <= 0):
        raise ValueError('Nonpositive moisture in nonlinear solve; reduce timestep.')
    return 7e-9 * np.exp(-.89 / c)


class Grid:
    def __init__(self, refinement=1, dimension=2, axial_refinement=None):
        # Every required 1 mm output radius is a node. Last 1 mm is finer.
        n = 2 * refinement
        self.r = np.r_[np.linspace(0, .019, 19*n+1),
                       np.linspace(.019, R, 8*n+1)[1:]]
        if dimension == 2:
            q = np.linspace(0, 1, 24*(axial_refinement or refinement)+1)
            self.z = HALF_LENGTH * (1 - (1-q)**2)
            zf = np.r_[0, (self.z[:-1]+self.z[1:])/2, HALF_LENGTH]
            dz = np.diff(zf)
        else:
            self.z = np.array([0.])
            dz = np.array([1.])  # per-unit axial length
        rf = np.r_[0, (self.r[:-1]+self.r[1:])/2, R]
        ar = np.diff(rf**2)/2  # volume / (2*pi)
        self.shape = (len(self.z), len(self.r))
        self.volume = (dz[:, None]*ar[None, :]).ravel()
        self.sqrtv = np.sqrt(self.volume)
        self.size = len(self.volume)
        ids = np.arange(self.size).reshape(self.shape)
        p = [ids[:, :-1].ravel()]
        qidx = [ids[:, 1:].ravel()]
        geom = [(dz[:, None]*rf[None, 1:-1] /
                 np.diff(self.r)[None, :]).ravel()]
        if dimension == 2:
            p.append(ids[:-1, :].ravel())
            qidx.append(ids[1:, :].ravel())
            geom.append((ar[None, :] / np.diff(self.z)[:, None]).ravel())
        self.p, self.q = np.concatenate(p), np.concatenate(qidx)
        self.geom = np.concatenate(geom)
        self.side = np.zeros(self.shape)
        self.side[:, -1] = R*dz
        self.end = np.zeros(self.shape)
        if dimension == 2:
            self.end[-1, :] = ar
        self.side, self.end = self.side.ravel(), self.end.ravel()
        self.boundary = self.side+self.end
        self.rows = np.r_[self.p, self.q, np.arange(self.size)]
        self.cols = np.r_[self.q, self.p, np.arange(self.size)]

    def operator(self, a, beta):
        if np.ndim(a):
            af = 2*a[self.p]*a[self.q]/(a[self.p]+a[self.q])
        else:
            af = a
        g = self.geom*af
        diagonal = -(np.bincount(self.p, weights=g, minlength=self.size)
                     + np.bincount(self.q, weights=g, minlength=self.size)
                     + beta*self.boundary)/self.volume
        offdiag = g/(self.sqrtv[self.p]*self.sqrtv[self.q])
        mat = coo_matrix((np.r_[offdiag, offdiag, diagonal],
                          (self.rows, self.cols)),
                         shape=(self.size, self.size)).tocsr()
        forcing = beta*self.boundary/self.sqrtv
        return mat, forcing


def integrate(refinement, dimension, dt, label, axial_refinement=None, startup=True):
    start = time.perf_counter()
    env = environment()
    grid = Grid(refinement, dimension, axial_refinement)
    temperature = np.full(grid.size, 28.)
    moisture = np.full(grid.size, 2.55)
    kt, bt = grid.operator(ALPHA, H/(RHO*CP))
    factor_cache = {}
    moisture_factor_cache = {}
    kreference, _ = grid.operator(float(diffusivity(np.array(2.55))), HM)
    temperature_output, moisture_output = [], []
    snapshots_t, snapshots_c = [], []
    balances = []
    cumulative_t, cumulative_c = 0., 0.
    max_balance_t = max_balance_c = 0.
    max_iter = total_iter = nsteps = cg_iterations = 0
    max_nonlinear_residual = 0.
    min_c, max_c, min_t, max_t = 2.55, 2.55, 28., 28.
    totalv = grid.volume.sum()
    output_ids = np.array([np.argmin(abs(grid.r-r)) for r in OUTPUT_R])
    assert np.max(abs(grid.r[output_ids]-OUTPUT_R)) < 1e-14

    def ambient(t):
        return tuple(np.interp(t, env[:, 0], env[:, j]) for j in (1, 2))

    def flux(u, boundary_value, beta):
        return beta*np.dot(grid.boundary, boundary_value-u)

    def advance(t0, step, theta):
        nonlocal temperature, moisture, cumulative_t, cumulative_c
        nonlocal max_iter, total_iter, nsteps, cg_iterations
        nonlocal max_balance_t, max_balance_c, max_nonlinear_residual
        nonlocal min_c, max_c, min_t, max_t
        old_t, old_c = temperature, moisture
        et0, ec0 = ambient(t0)
        et1, ec1 = ambient(t0+step)
        x0 = grid.sqrtv*old_t
        cache_key = (round(step, 12), theta)
        if cache_key not in factor_cache:
            factor_cache[cache_key] = splu((eye(grid.size, format='csc')-
                                           theta*step*kt).tocsc())
        rhs = x0+(1-theta)*step*(kt@x0+bt*et0)+theta*step*bt*et1
        temperature = factor_cache[cache_key].solve(rhs)/grid.sqrtv

        kc0, bc = grid.operator(diffusivity(old_c), HM)
        xc0 = grid.sqrtv*old_c
        rhs_c = xc0+(1-theta)*step*(kc0@xc0+bc*ec0)+theta*step*bc*ec1
        guess = old_c.copy()
        for iteration in range(1, 31):
            km, _ = grid.operator(diffusivity(guess), HM)
            am = eye(grid.size, format='csr')-theta*step*km
            if cache_key not in moisture_factor_cache:
                moisture_factor_cache[cache_key] = splu(
                    (eye(grid.size, format='csc')-theta*step*kreference).tocsc())
            preconditioner = LinearOperator(am.shape,
                matvec=moisture_factor_cache[cache_key].solve)
            counter = [0]
            def callback(_):
                counter[0] += 1
            xn, info = cg(am, rhs_c, x0=grid.sqrtv*guess,
                          rtol=2e-13, atol=1e-16,
                          M=preconditioner, maxiter=1500, callback=callback)
            cg_iterations += counter[0]
            if info:
                raise RuntimeError(f'CG failed: {info}')
            candidate = xn/grid.sqrtv
            change = np.max(abs(candidate-guess))
            guess = candidate
            if change < 2e-10:
                kcheck, _ = grid.operator(diffusivity(guess), HM)
                residual = np.max(abs((xn-theta*step*(kcheck@xn)-rhs_c)
                                      /grid.sqrtv))
                if residual < 2e-9:
                    break
        else:
            raise RuntimeError('Picard iteration did not converge')
        moisture = guess
        max_nonlinear_residual = max(max_nonlinear_residual, residual)
        max_iter = max(max_iter, iteration)
        total_iter += iteration
        nsteps += 1
        cumulative_t += step*((1-theta)*flux(old_t, et0, H/(RHO*CP))+
                              theta*flux(temperature, et1, H/(RHO*CP)))
        cumulative_c += step*((1-theta)*flux(old_c, ec0, HM)+
                              theta*flux(moisture, ec1, HM))
        balance_t = np.dot(grid.volume, temperature-28.)-cumulative_t
        balance_c = np.dot(grid.volume, moisture-2.55)-cumulative_c
        max_balance_t = max(max_balance_t, abs(balance_t)/totalv)
        max_balance_c = max(max_balance_c, abs(balance_c)/totalv)
        min_t, max_t = min(min_t, temperature.min()), max(max_t, temperature.max())
        min_c, max_c = min(min_c, moisture.min()), max(max_c, moisture.max())
        if moisture.min() <= 0 or moisture.max() > 2.55000001:
            raise RuntimeError('Moisture maximum principle check failed')
        if temperature.min() < 27.999999 or temperature.max() > max(et0, et1, old_t.max())+1e-7:
            raise RuntimeError('Thermal maximum principle check failed')

    t = 0.
    # Resolve the incompatible initial Robin data before using the main step.
    first_dt = dt/5000 if startup else dt
    advance(0., first_dt/2, 1.)
    advance(first_dt/2, first_dt/2, 1.)
    t = first_dt
    for target in range(1, 1801):
        while t < target-1e-10:
            step = min(dt, target-t)
            if startup:
                for end, divisor in [(.001, 5000), (.01, 500), (.1, 50),
                                     (1., 25), (10., 10), (60., 4)]:
                    if t < end-1e-10:
                        step = min(step, dt/divisor, end-t)
                        break
            step = round(step, 12)
            advance(t, step, .5)
            t += step
        temperature_output.append(temperature.reshape(grid.shape)[0, output_ids].copy())
        moisture_output.append(moisture.reshape(grid.shape)[0, output_ids].copy())
        balances.append([target, np.dot(grid.volume, temperature)/totalv,
                         np.dot(grid.volume, moisture)/totalv,
                         cumulative_t/totalv, cumulative_c/totalv])
        if target in TABLE_TIMES:
            snapshots_t.append(temperature.reshape(grid.shape).copy())
            snapshots_c.append(moisture.reshape(grid.shape).copy())
        if target % 300 == 0:
            print(f'{label}: t={target}, nodes={grid.size}, '
                  f'elapsed={time.perf_counter()-start:.1f}s, max_iter={max_iter}', flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT/f'{label}.npz', r=grid.r, z=grid.z,
                        temperature=np.asarray(temperature_output),
                        moisture=np.asarray(moisture_output),
                        snapshots_t=np.asarray(snapshots_t),
                        snapshots_c=np.asarray(snapshots_c),
                        balances=np.asarray(balances), times=np.arange(1, 1801))
    summary = {
        'label': label, 'dimension': dimension, 'refinement': refinement,
        'dt_s': dt, 'grid_shape_z_r': grid.shape,
        'graded_startup': startup, 'first_dt_s': first_dt,
        'axial_refinement': axial_refinement or refinement,
        'radial_min_spacing_m': float(np.diff(grid.r).min()),
        'axial_min_spacing_m': float(np.diff(grid.z).min()) if dimension == 2 else None,
        'steps': nsteps, 'max_picard_iterations': max_iter,
        'mean_picard_iterations': total_iter/nsteps,
        'cg_iterations': cg_iterations,
        'max_nonlinear_residual_C': float(max_nonlinear_residual),
        'max_heat_balance_residual_mean_Celsius': float(max_balance_t),
        'max_water_balance_residual_mean_kg_kg': float(max_balance_c),
        'temperature_bounds': [float(min_t), float(max_t)],
        'moisture_bounds': [float(min_c), float(max_c)],
        'final_center_T_C': [float(temperature[0]), float(moisture[0])],
        'final_midplane_surface_T_C': [float(temperature[len(grid.r)-1]),
                                        float(moisture[len(grid.r)-1])],
        'elapsed_seconds': time.perf_counter()-start,
        'input_sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__},
        'status': 'success',
    }
    (OUT/f'{label}_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=True), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--refinement', type=int, default=16)
    parser.add_argument('--dimension', type=int, choices=(1, 2), default=2)
    parser.add_argument('--axial-refinement', type=int, default=1)
    parser.add_argument('--legacy-startup', action='store_true')
    parser.add_argument('--dt', type=float, default=.5)
    parser.add_argument('--label', default='q1_main_r16')
    args = parser.parse_args()
    if args.dt <= 0 or args.dt > 1 or abs(round(1/args.dt)-1/args.dt) > 1e-10:
        parser.error('dt must divide one second and lie in (0,1]')
    integrate(args.refinement, args.dimension, args.dt, args.label,
              args.axial_refinement, not args.legacy_startup)
