"""Command-line interface for extension-bsf."""

from __future__ import annotations

import argparse
import json
import sys

from .bilayer import BilayerError, generate_scan, write_bilayer


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="extension-bsf")
    sub = parser.add_subparsers(dest="command", required=True)

    bilayer = sub.add_parser("bilayer", help="build one validated bilayer")
    bilayer.add_argument("input")
    bilayer.add_argument("output")
    bilayer.add_argument("--interlayer-distance", type=float, required=True)
    bilayer.add_argument("--vacuum", type=float, required=True)
    bilayer.add_argument("--shift", nargs=2, type=float, default=(0.0, 0.0))
    bilayer.add_argument("--minimum-distance", type=float, default=1.0)
    bilayer.add_argument("--max-cross-layer-pairs", type=int, default=4_000_000)
    bilayer.add_argument("--force", action="store_true")

    scan = sub.add_parser("scan", help="build a half-open periodic stacking grid")
    scan.add_argument("input")
    scan.add_argument("output_directory")
    scan.add_argument("--grid", nargs=2, type=int, required=True)
    scan.add_argument("--interlayer-distance", type=float, required=True)
    scan.add_argument("--vacuum", type=float, required=True)
    scan.add_argument("--minimum-distance", type=float, default=1.0)
    scan.add_argument("--max-scan-points", type=int, default=10_000)
    scan.add_argument("--max-total-atoms", type=int, default=2_000_000)
    scan.add_argument("--max-cross-layer-pairs", type=int, default=4_000_000)
    scan.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "bilayer":
            result = write_bilayer(
                args.input,
                args.output,
                interlayer_distance=args.interlayer_distance,
                vacuum=args.vacuum,
                shift=tuple(args.shift),
                minimum_distance=args.minimum_distance,
                max_cross_layer_pairs=args.max_cross_layer_pairs,
                force=args.force,
            )
        else:
            result = generate_scan(
                args.input,
                args.output_directory,
                grid=tuple(args.grid),
                interlayer_distance=args.interlayer_distance,
                vacuum=args.vacuum,
                minimum_distance=args.minimum_distance,
                max_scan_points=args.max_scan_points,
                max_total_atoms=args.max_total_atoms,
                max_cross_layer_pairs=args.max_cross_layer_pairs,
                force=args.force,
            )
    except (BilayerError, OSError) as exc:
        print(f"extension-bsf: error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
