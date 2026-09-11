"""Independent Q1 checks: exact 1D heat series and FVM identities.

Run with: python q1_independent_reference.py --with-rerun --save
The optional repeat runs the official 1D baseline in temporary storage.
The Bessel-series heat solution is derived independently of its discretization.
"""
from __future__ import annotations

import hashlib
import importlib.util
import io
import json
from pathlib import Path
from contextlib import redirect_stdout
from tempfile import TemporaryDirectory

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1

SOLVER = Path(__file__).resolve().parents[1] / "0910_第一问_热湿模型_v01_Codex.py"
MAIN_LABEL = "q1_main_r16"
BASELINE_LABEL = "q1_baseline_r16"
RUN_LABELS = (MAIN_LABEL, BASELINE_LABEL, "q1_reference_r16", "q1_spacecheck_r32")
spec = importlib.util.spec_from_file_location("reviewed_q1", SOLVER)
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)


def exact_radial_heat(modes=256):
    """Return T(t,r), t=1..1800, for piecewise-linear measured ambient T.

    For lambda*J1(lambda)=Bi*J0(lambda), the unit initial profile has
    coefficient 2 J1/[lambda (J0^2+J1^2)].  q_n obeys q_n'=Tinf'-beta*q_n.
    The recurrence integrates every linear ambient segment analytically.
    """
    env = solver.environment()
    bi = solver.H * solver.R / solver.K
    root_function = lambda x: x * j1(x) - bi * j0(x)
    brackets = np.arange(1e-8, (modes + 2) * np.pi, 0.1)
    roots = np.array([
        brentq(root_function, left, right)
        for left, right in zip(brackets[:-1], brackets[1:])
        if root_function(left) * root_function(right) < 0
    ])[:modes]
    assert len(roots) == modes
    coeff = 2*j1(roots)/(roots*(j0(roots)**2+j1(roots)**2))
    rates = solver.ALPHA*roots**2/solver.R**2
    basis = j0(np.outer(solver.OUTPUT_R/solver.R, roots))
    q = np.full(modes, env[0, 1]-28.)
    decay = np.exp(-rates)
    response = -np.expm1(-rates)/rates
    previous = env[0, 1]
    temperatures = []
    for second in range(1, 1801):
        ambient = np.interp(second, env[:, 0], env[:, 1])
        q = q*decay+(ambient-previous)*response
        temperatures.append(ambient-basis@(coeff*q))
        previous = ambient
    return np.array(temperatures), roots, coeff


def check_geometry_and_operator():
    grid = solver.Grid(16, 2, axial_refinement=1)
    coefficient = 2.
    matrix, _ = grid.operator(coefficient, 0.)
    physical = matrix.multiply(grid.sqrtv[None, :]).multiply(
        (1/grid.sqrtv)[:, None]).tocsr()
    radial_actual = float(physical[0, 1])
    radial_expected = 4*coefficient/np.diff(grid.r)[0]**2
    axial_actual = float(physical[0, len(grid.r)])
    axial_expected = 2*coefficient/np.diff(grid.z)[0]**2
    # An arbitrary nonuniform, positive profile exercises variable D faces.
    moisture = np.linspace(.2, 2.5, grid.size)
    matrix, forcing = grid.operator(solver.diffusivity(moisture), solver.HM)
    rhs = matrix@(grid.sqrtv*moisture)+forcing*.02
    balance = np.dot(grid.sqrtv, rhs)-solver.HM*np.dot(grid.boundary, .02-moisture)
    result = {
        "nodes": grid.size,
        "radial_center_actual": radial_actual,
        "radial_center_expected": float(radial_expected),
        "axial_center_actual": axial_actual,
        "axial_center_expected": float(axial_expected),
        "matrix_symmetry_max": float(abs(matrix-matrix.T).max()),
        "arbitrary_state_mass_balance_abs": float(abs(balance)),
        "volume_scaled": float(grid.volume.sum()),
        "volume_scaled_expected": solver.R**2/2*solver.HALF_LENGTH,
        "external_area_scaled": float(grid.boundary.sum()),
        "external_area_scaled_expected": solver.R*solver.HALF_LENGTH+solver.R**2/2,
    }
    assert np.isclose(radial_actual, radial_expected, rtol=1e-14)
    assert np.isclose(axial_actual, axial_expected, rtol=1e-14)
    assert result["matrix_symmetry_max"] < 1e-14
    assert result["arbitrary_state_mass_balance_abs"] < 1e-19
    assert np.isclose(result["volume_scaled"], result["volume_scaled_expected"], rtol=1e-14)
    assert np.isclose(result["external_area_scaled"], result["external_area_scaled_expected"], rtol=1e-14)
    return result


def exact_axisymmetric_heat(radial_nodes, axial_nodes, radial_modes=256, axial_modes=1024):
    """Finite-cylinder product eigenfunctions, evaluated at seven table times."""
    env = solver.environment()
    _, roots, radial_coeff = exact_radial_heat(radial_modes)
    axial_biot = solver.H*solver.HALF_LENGTH/solver.K
    axial_roots = np.array([
        brentq(lambda x: x*np.tan(x)-axial_biot,
               n*np.pi+1e-9, (n+.5)*np.pi-1e-9)
        for n in range(axial_modes)
    ])
    axial_coeff = 4*np.sin(axial_roots)/(2*axial_roots+np.sin(2*axial_roots))
    rates = solver.ALPHA*((axial_roots[:, None]/solver.HALF_LENGTH)**2
                         +(roots[None, :]/solver.R)**2)
    axial_basis = np.cos(np.outer(axial_nodes/solver.HALF_LENGTH, axial_roots))
    radial_basis = j0(np.outer(radial_nodes/solver.R, roots))
    state = np.full((axial_modes, radial_modes), env[0, 1]-28.)
    knots = sorted(set(env[env[:, 0] <= 1800, 0].tolist()+solver.TABLE_TIMES.tolist()))
    result = []
    last, previous = 0., env[0, 1]
    for current in knots[1:]:
        ambient = np.interp(current, env[:, 0], env[:, 1])
        step = current-last
        state = (state*np.exp(-step*rates)
                 +(ambient-previous)*(-np.expm1(-step*rates))/(step*rates))
        if current in solver.TABLE_TIMES:
            result.append(ambient-axial_basis@(
                axial_coeff[:, None]*state*radial_coeff[None, :])@radial_basis.T)
        last, previous = current, ambient
    return np.array(result)


def inspect_saved_runs(reference):
    results = {}
    for label in RUN_LABELS:
        path = solver.OUT/(label+".npz")
        summary_path = path.with_name(path.stem+"_summary.json")
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        with np.load(path) as data:
            temperatures, moisture = data["temperature"], data["moisture"]
            assert temperatures.shape == moisture.shape == (1800, 21)
            assert np.array_equal(data["times"], np.arange(1, 1801))
            assert np.isfinite(temperatures).all() and np.isfinite(moisture).all()
            assert moisture.min() > 0 and moisture.max() <= 2.5500001
            error = temperatures-reference
            b = data["balances"]
            heat_closure = np.max(abs((b[:, 1]-28)-b[:, 3]))
            water_closure = np.max(abs((b[:, 2]-2.55)-b[:, 4]))
            final = np.array([temperatures[-1, 0], moisture[-1, 0]])
            assert np.allclose(final, summary["final_center_T_C"], rtol=0, atol=1e-12)
            snapshot_final = np.array([data["snapshots_t"][-1, 0, -1],
                                       data["snapshots_c"][-1, 0, -1]])
            assert np.allclose(snapshot_final, summary["final_midplane_surface_T_C"], rtol=0, atol=1e-12)
            results[path.stem] = {
                "dimension": summary["dimension"],
                "refinement": summary["refinement"],
                "dt_s": summary["dt_s"],
                "heat_error_analytic_max_Celsius": float(abs(error).max()),
                "heat_error_analytic_t100_1800_max_Celsius": float(abs(error[99:]).max()),
                "heat_error_analytic_final_max_Celsius": float(abs(error[-1]).max()),
                "saved_heat_balance_max_mean_Celsius": float(heat_closure),
                "saved_water_balance_max_mean_kg_kg": float(water_closure),
                "nonlinear_residual_max_kg_kg": summary["max_nonlinear_residual_C"],
                "final_center_T_C": final.tolist(),
                "final_midplane_surface_T_C": snapshot_final.tolist(),
                "npz_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
    return results


def inspect_axisymmetric_heat():
    # One complete 2D run suffices for this independent boundary/axis test.
    path = solver.OUT/(MAIN_LABEL+".npz")
    if not path.exists():
        return {"status": "not_run", "reason": f"{MAIN_LABEL}.npz does not exist"}
    with np.load(path) as data:
        high = exact_axisymmetric_heat(data["r"], data["z"], axial_modes=1024)
        low = exact_axisymmetric_heat(data["r"], data["z"], axial_modes=512)
        truncation = float(abs(high-low).max())
        assert truncation < 5e-7
        error = data["snapshots_t"]-high
        return {
            "run": path.stem,
            "radial_modes": 256,
            "axial_modes": 1024,
            "table_times_s": solver.TABLE_TIMES.tolist(),
            "axial_mode_doubling_max_Celsius": truncation,
            "full_2d_table_snapshots_error_max_Celsius": float(abs(error).max()),
            "final_snapshot_error_max_Celsius": float(abs(error[-1]).max()),
            "final_corner_reference_Celsius": float(high[-1, -1, -1]),
            "final_corner_computed_Celsius": float(data["snapshots_t"][-1, -1, -1]),
        }


def check_baseline_reproducibility():
    original_output = solver.OUT
    with np.load(original_output/(BASELINE_LABEL+".npz")) as baseline:
        with TemporaryDirectory(prefix="q1_independent_repeat_") as temporary:
            try:
                solver.OUT = Path(temporary)
                with redirect_stdout(io.StringIO()):
                    repeat_summary = solver.integrate(16, 1, .5, "repeat", startup=True)
                with np.load(solver.OUT/"repeat.npz") as repeated:
                    differences = {
                        field: float(abs(baseline[field]-repeated[field]).max())
                        for field in ("temperature", "moisture", "snapshots_t", "snapshots_c", "balances")
                    }
                    assert max(differences.values()) < 1e-12
            finally:
                solver.OUT = original_output
    return {"baseline": BASELINE_LABEL, "max_abs_differences": differences,
            "elapsed_seconds": repeat_summary["elapsed_seconds"],
            "temporary_outputs_removed": True}


def run_checks(with_rerun=False):
    reference, roots, coeff = exact_radial_heat(512)
    reference_half, _, _ = exact_radial_heat(256)
    truncation = float(abs(reference-reference_half).max())
    assert truncation < 2e-8
    result = {
        "canonical_main": MAIN_LABEL,
        "canonical_baseline": BASELINE_LABEL,
        "solver_sha256": hashlib.sha256(SOLVER.read_bytes()).hexdigest(),
        "geometry_and_operator": check_geometry_and_operator(),
        "analytic_reference": {
            "modes": 512,
            "comparison_256_512_modes_max_Celsius": truncation,
            "Biot": solver.H*solver.R/solver.K,
            "first_three_roots": roots[:3].tolist(),
            "first_three_coefficients": coeff[:3].tolist(),
            "final_center_and_surface_Celsius": reference[-1, [0, -1]].tolist(),
        },
        "saved_runs": inspect_saved_runs(reference),
        "finite_cylinder_analytic_reference": inspect_axisymmetric_heat(),
    }
    if with_rerun:
        result["baseline_reproducibility"] = check_baseline_reproducibility()
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--with-rerun", action="store_true")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    results = run_checks(args.with_rerun)
    payload = json.dumps(results, ensure_ascii=False, indent=2)
    if args.save:
        Path(__file__).with_name("q1_independent_checks.json").write_text(
            payload+"\n", encoding="utf-8")
    print(payload)
