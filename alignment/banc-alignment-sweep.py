"""
Optic Lobe Alignment: Parameter Sweep

Two-stage parameter optimization:
  Stage 1: Fast one-at-a-time screening with cosine metric (2s/iter)
  Stage 2: Validate top configs with ensemble metric (20s/iter)

Imports run_alignment from banc-alignment-run.py.

Usage:
  python alignment/banc-alignment-sweep.py --side right
  python alignment/banc-alignment-sweep.py --side right --stage 2
"""

import argparse
import functools
import sys
import time
import os
import pandas as pd

# Add parent dir to path so we can import the align module
sys.path.insert(0, os.path.dirname(__file__))
from importlib import import_module

print = functools.partial(print, flush=True)


def run_config(data, config, label, screen_iters=15, data_dir="data/optic_lobe"):
    """Run alignment with given config, return metrics dict."""
    banc_meta, target_meta, banc_el, target_el, nblast, seeds = data

    t0 = time.time()
    _, metrics = align_mod.run_alignment(
        banc_meta, target_meta, banc_el, target_el, nblast, seeds,
        data_dir=data_dir, **config)
    elapsed = time.time() - t0

    metrics["label"] = label
    metrics["elapsed"] = elapsed
    for k, v in config.items():
        metrics[f"param_{k}"] = v
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Parameter sweep for optic lobe alignment")
    parser.add_argument("--side", default="right", choices=["right", "left"])
    parser.add_argument("--data-dir", default="data/optic_lobe")
    parser.add_argument("--stage", type=int, default=0, choices=[0, 1, 2, 3],
                        help="0=all stages, 1=fast screening, 2=ensemble validation, "
                             "3=individual profile sweep")
    parser.add_argument("--screen-iters", type=int, default=15,
                        help="Iterations for screening runs (default: 15)")
    parser.add_argument("--validate-iters", type=int, default=25,
                        help="Iterations for validation runs (default: 25)")
    args = parser.parse_args()

    # Import alignment module (hyphenated name needs importlib)
    global align_mod
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "align", os.path.join(os.path.dirname(__file__), "banc-alignment-run.py"))
    align_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(align_mod)

    # Load data once
    print("=== Loading data ===")
    data = align_mod.load_data(args.side, args.data_dir)

    all_results = []

    # Default config
    defaults = dict(
        max_iter=args.screen_iters, convergence_pct=0.5,
        tau_start=2.0, tau_end=0.5,
        alpha_start=0.3, alpha_end=0.8,
        nblast_threshold=0.25,
        capacity_scale=1.5, sinkhorn_iters=5,
        metric="cosine", prior_weight=0.0,
        hop2_weight=0.5, ensemble_top_k=30, ensemble_blend=0.5,
        nt_weight=0.0, ind_weight=0.0, ind_metric="auto",
    )

    # =========================================================
    # Stage 1: Fast screening with cosine (one-at-a-time)
    # =========================================================
    if args.stage in (0, 1):
        print("\n" + "=" * 60)
        print("STAGE 1: Fast screening (cosine, %d iters)" % args.screen_iters)
        print("=" * 60)

        sweeps = {
            "capacity_scale": [1.0, 1.2, 1.5, 2.0],
            "hop2_weight": [0.0, 0.25, 0.5, 1.0],
            "nblast_threshold": [0.15, 0.25, 0.35],
            "tau_start": [1.0, 2.0, 4.0],
            "alpha_start": [0.1, 0.3, 0.5],
            "alpha_end": [0.6, 0.8, 0.95],
            "ind_weight": [0.0, 0.1, 0.2, 0.3, 0.5],
        }

        # Baseline
        print("\n--- Baseline (defaults) ---")
        m = run_config(data, defaults, "baseline", args.screen_iters, args.data_dir)
        print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
              f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
        all_results.append(m)

        # Sweep each parameter
        best_per_param = {}
        for param, values in sweeps.items():
            print(f"\n--- Sweeping {param}: {values} ---")
            param_results = []
            for val in values:
                config = defaults.copy()
                config[param] = val
                label = f"{param}={val}"
                m = run_config(data, config, label, args.screen_iters, args.data_dir)
                print(f"  {label:30s}  holdout={m['best_holdout']:.1f}%, "
                      f"Mi1={m['mi1_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, "
                      f"{m['elapsed']:.0f}s")
                all_results.append(m)
                param_results.append((m["best_holdout"], val, m))

            # Best value for this param
            param_results.sort(key=lambda x: -x[0])
            best_val = param_results[0][1]
            best_per_param[param] = best_val
            print(f"  Best: {param}={best_val} ({param_results[0][0]:.1f}%)")

        # Combined best from stage 1
        print(f"\n--- Stage 1 combined best: {best_per_param} ---")
        combined = defaults.copy()
        combined.update(best_per_param)
        m = run_config(data, combined, "stage1_combined", args.screen_iters, args.data_dir)
        print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
              f"rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
        all_results.append(m)

    # =========================================================
    # Stage 2: Validate with ensemble metric
    # =========================================================
    if args.stage in (0, 2):
        print("\n" + "=" * 60)
        print("STAGE 2: Ensemble validation (%d iters)" % args.validate_iters)
        print("=" * 60)

        # If stage 1 ran, use its best params; otherwise use defaults
        if args.stage == 0:
            base = defaults.copy()
            base.update(best_per_param)
        else:
            base = defaults.copy()

        base["metric"] = "ensemble"
        base["max_iter"] = args.validate_iters

        configs = [
            ("ensemble_defaults", {**defaults, "metric": "ensemble", "max_iter": args.validate_iters}),
            ("ensemble_best_s1", base),
        ]

        # Also try a few ensemble-specific variations
        for ek in [10, 30, 50]:
            c = base.copy()
            c["ensemble_top_k"] = ek
            configs.append((f"ensemble_topk={ek}", c))

        for eb in [0.3, 0.5, 0.7]:
            c = base.copy()
            c["ensemble_blend"] = eb
            configs.append((f"ensemble_blend={eb}", c))

        for label, config in configs:
            print(f"\n--- {label} ---")
            m = run_config(data, config, label, args.validate_iters, args.data_dir)
            print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
                  f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
            all_results.append(m)

    # =========================================================
    # Stage 3: Individual profile sweep (ensemble + ind_weight)
    # =========================================================
    if args.stage in (0, 3):
        print("\n" + "=" * 60)
        print("STAGE 3: Individual profile sweep (ensemble, %d iters)" % args.validate_iters)
        print("=" * 60)

        # Use best known config as base
        best_known = dict(
            max_iter=args.validate_iters, convergence_pct=0.5,
            tau_start=4.0, tau_end=0.5,
            alpha_start=0.05, alpha_end=0.95,
            nblast_threshold=0.15,
            capacity_scale=1.5, sinkhorn_iters=5,
            metric="ensemble", prior_weight=0.0,
            hop2_weight=1.0, ensemble_top_k=30, ensemble_blend=0.3,
            nt_weight=0.0, ind_weight=0.0, ind_metric="auto",
        )

        # Baseline: no individual scoring
        print("\n--- Baseline (ind_weight=0, ensemble) ---")
        m = run_config(data, best_known, "ind_baseline", args.validate_iters, args.data_dir)
        print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
              f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
        all_results.append(m)

        # Sweep ind_weight
        for iw in [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]:
            c = best_known.copy()
            c["ind_weight"] = iw
            label = f"ind_weight={iw}"
            print(f"\n--- {label} ---")
            m = run_config(data, c, label, args.validate_iters, args.data_dir)
            print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
                  f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
            all_results.append(m)

        # Sweep ind_weight with cosine-only metric (to see if Jaccard matters for ind)
        for iw in [0.1, 0.2, 0.3]:
            c = best_known.copy()
            c["ind_weight"] = iw
            c["metric"] = "cosine"
            label = f"ind_weight={iw}_cosine"
            print(f"\n--- {label} ---")
            m = run_config(data, c, label, args.validate_iters, args.data_dir)
            print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
                  f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
            all_results.append(m)

        # Sweep ind_metric=jaccard (pure Jaccard for individual scoring, ensemble for centroids)
        for iw in [0.1, 0.2, 0.3]:
            c = best_known.copy()
            c["ind_weight"] = iw
            c["ind_metric"] = "jaccard"
            label = f"ind_weight={iw}_ind_jaccard"
            print(f"\n--- {label} ---")
            m = run_config(data, c, label, args.validate_iters, args.data_dir)
            print(f"  holdout={m['best_holdout']:.1f}%, Mi1={m['mi1_at_best']:.1f}%, "
                  f"NT={m['nt_at_best']:.1f}%, rho={m['type_rho_at_best']:.3f}, {m['elapsed']:.0f}s")
            all_results.append(m)

    # Save results
    results_df = pd.DataFrame(all_results)
    out_file = f"{args.data_dir}/sweep_results_{args.side}.csv"
    results_df.to_csv(out_file, index=False)
    print(f"\nSaved {len(results_df)} results to {out_file}")

    # Print summary table
    print("\n=== SUMMARY (sorted by holdout) ===")
    cols = ["label", "best_holdout", "mi1_at_best", "nt_at_best", "type_rho_at_best", "elapsed"]
    summary = results_df[cols].sort_values("best_holdout", ascending=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
