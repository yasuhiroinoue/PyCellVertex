#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import math
import random
import sys
from pathlib import Path

# Setup import path for local src directory
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from celldynamics.constants import L_THRESHOLD, STEP_RECONNECT_DEFAULT
from celldynamics.division import cell_division, is_convex, update_centers
from celldynamics.init_plain import InitParams, init_plain
from celldynamics.intersection import cell_intersection
from celldynamics.rearrange import cell_rearrange2
from celldynamics.sim_step import DELTA_TIME, motion_vertex_second_step, update_pulse
from celldynamics.vec import Vec3
from celldynamics.vtk_output import dump_vtk_snapshot

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Celldynamics-Python Main Simulation Entrypoint",
        fromfile_prefix_chars="@",
    )
    p.add_argument(
        "--out-dir", type=Path, default=Path("output"), help="Output directory for VTK files"
    )
    p.add_argument("--steps", type=int, default=1000, help="Total number of simulation steps")

    # Init settings
    p.add_argument("--num-x", type=int, default=5, help="Number of cells in X direction")
    p.add_argument("--num-y", type=int, default=5, help="Number of cells in Y direction")
    p.add_argument("--k-area", type=float, default=1.0, help="Area elasticity constant")
    p.add_argument("--area-eq", type=float, default=2.60, help="Equilibrium area")
    p.add_argument("--k2-length", type=float, default=1.0, help="Line tension constant K2")
    p.add_argument("--k1-length", type=float, default=-6.0, help="Line tension constant K1")
    p.add_argument("--k1-pcp", type=float, default=0.0, help="PCP line tension constant")
    p.add_argument("--length-eq", type=float, default=0.0, help="Equilibrium length")
    p.add_argument("--pulse-t", type=float, default=55.0, help="Pulse period T")
    p.add_argument("--phase-x", type=float, default=0.7854, help="Phase shift in X")
    p.add_argument("--phase-y", type=float, default=1.5708, help="Phase shift in Y")

    # Simulation settings
    p.add_argument("--power-pcp", type=float, default=2.0, help="Power parameter for PCP")
    p.add_argument("--enable-division", action=argparse.BooleanOptionalAction, default=True, help="Enable cell division")
    p.add_argument(
        "--division-time", type=float, default=1000.0, help="Physical time before division eligibility"
    )
    p.add_argument(
        "--division-stagger-frac",
        type=float,
        default=0.25,
        help="Initial/daughter cell_time randomization range as a fraction of division-time",
    )
    p.add_argument(
        "--enable-rearrange", action=argparse.BooleanOptionalAction, default=True, help="Enable T1 transitions (rearrangement)"
    )
    p.add_argument("--enable-intersection", action=argparse.BooleanOptionalAction, default=True, help="Enable cell intersections")
    p.add_argument(
        "--step-reconnect",
        type=int,
        default=STEP_RECONNECT_DEFAULT,
        help="Interval steps for reconnect/rearrange",
    )
    p.add_argument(
        "--vacant-area-th",
        type=float,
        default=0.01,
        help="Area threshold for triangular vacant region fix",
    )

    # Output settings
    p.add_argument(
        "--dump-vtk", action="store_true", help="Emit C++-format 2dv_line/2dv_face VTK files"
    )
    p.add_argument("--vtk-step", type=int, default=100, help="Output VTK every N steps")
    p.add_argument(
        "--vtk-step-scale",
        type=int,
        default=10000,
        help="Filename step scaling (matches C++ default: step*10000)",
    )
    p.add_argument(
        "--vtk-alive-only", action="store_true", help="Filter out dead lines/vertices from VTK output"
    )
    p.add_argument(
        "--vtk-mode", type=str, choices=["raw", "alive", "both"], default="raw", help="VTK output mode"
    )

    p.add_argument("--config", type=Path, help="Optional JSON config file with parameters")

    # Parse config first if present
    initial_args = p.parse_known_args()[0]
    if initial_args.config and initial_args.config.exists():
        import json

        with open(initial_args.config) as f:
            config_data = json.load(f)
        p.set_defaults(**config_data)

    args = p.parse_args()

    if args.division_stagger_frac < 0.0:
        p.error("--division-stagger-frac must be >= 0")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Starting simulation. Output directory: {args.out_dir}")

    # Initialize state
    state = init_plain(
        InitParams(
            num_x=args.num_x,
            num_y=args.num_y,
            k_area=args.k_area,
            area_eq=args.area_eq,
            k2_length=args.k2_length,
            k1_length=args.k1_length,
            k1_pcp_length=args.k1_pcp,
            length_eq=args.length_eq,
            pulse_t=args.pulse_t,
            phase_x=args.phase_x,
            phase_y=args.phase_y,
        )
    )
    logging.info(
        f"Initialized grid: {args.num_x}x{args.num_y} with {len(state.p_c)} cells, {len(state.p_v)} vertices, {len(state.p_l)} lines."
    )

    if args.enable_division:
        for c in state.p_c:
            c.cell_time = random.uniform(0.0, args.division_time * args.division_stagger_frac)

    event_log: list[str] = []

    if args.dump_vtk:
        dump_vtk_snapshot(state, step=0, out_dir=args.out_dir, step_scale=args.vtk_step_scale, alive_only=args.vtk_alive_only, vtk_mode=args.vtk_mode)
        logging.info("Wrote initial VTK (step 0).")

    for step in range(1, args.steps + 1):
        update_centers(state)

        if (
            args.enable_intersection
            and args.step_reconnect > 0
            and (step % args.step_reconnect) == (args.step_reconnect // 2)
        ):
            cell_intersection(state, step=step, event_log=event_log)
            update_centers(state)

        if args.enable_rearrange and args.step_reconnect > 0 and (step % args.step_reconnect) == 0:
            cell_rearrange2(state, threshold=L_THRESHOLD, step=step, event_log=event_log, vacant_area_th=args.vacant_area_th)
            update_centers(state)

        if args.enable_division and state.p_c:
            for cidx in range(len(state.p_c) - 1, -1, -1):
                if state.p_c[cidx].cell_time >= args.division_time and is_convex(
                    state, cidx
                ):
                    theta = random.uniform(0, 2 * math.pi)
                    axis = Vec3(math.cos(theta), math.sin(theta), 0.0)
                    try:
                        cell_division(state, cidx, axis)
                        event_log.append(f"{step},division,{cidx}")
                        state.p_c[cidx].cell_time = random.uniform(0.0, args.division_time * args.division_stagger_frac)
                        state.p_c[-1].cell_time = random.uniform(0.0, args.division_time * args.division_stagger_frac)
                        update_centers(state)
                    except RuntimeError:
                        pass

        update_pulse(state)

        # Force calculation and vertex position update
        motion_vertex_second_step(state, power_pcp=args.power_pcp, delta_time=DELTA_TIME)

        if args.dump_vtk and (step % args.vtk_step == 0 or step == args.steps):
            dump_vtk_snapshot(
                state, step=step, out_dir=args.out_dir, step_scale=args.vtk_step_scale, alive_only=args.vtk_alive_only, vtk_mode=args.vtk_mode
            )
            logging.info(f"Wrote VTK at step {step} (time {step * DELTA_TIME:.2f})")

    # Output event log if needed
    if event_log:
        log_path = args.out_dir / "event_log.txt"
        log_path.write_text("\n".join(event_log) + "\n", encoding="utf-8")
        logging.info(f"Wrote event log to {log_path}")

    logging.info(f"Simulation completed successfully ({args.steps} steps).")


if __name__ == "__main__":
    main()
