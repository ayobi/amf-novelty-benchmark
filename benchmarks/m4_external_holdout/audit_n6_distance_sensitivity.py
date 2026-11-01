#!/usr/bin/env python3

"""Compare directional and symmetric ITS distances for N6."""

import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

M4 = Path(__file__).resolve().parent
PROJECTS = M4.parents[2]

sys.path.insert(0, str(M4))

import external_holdout as v1
from its_rule import ItsRule, dist

DECISIONS = (
    PROJECTS / "amfsplit"
    / "parity/naming/n6_decision_parity/python_decisions.json"
)

MODES = ("query_first", "sym_min", "sym_mean", "sym_max")


def main():
    with DECISIONS.open() as fh:
        fixture = json.load(fh)

    data = v1.Data(v1.FROZEN / "inputs/benchmark")
    rule = ItsRule(data, fixture["quantile"])

    sequences = data.fa["full_its"]
    ids = rule.its_ok

    all_by_culture = defaultdict(list)
    main_by_culture = defaultdict(list)

    for accession in ids:
        culture = data.cop[accession]["culture"]
        all_by_culture[culture].append(accession)

        if accession in rule.mains:
            main_by_culture[culture].append(accession)

    # Compute both directions once per accession pair.
    pair_cache = {}

    def directed(x, y):
        key = (x, y)
        if key not in pair_cache:
            pair_cache[key] = dist(sequences[x], sequences[y])
        return pair_cache[key]

    def pair_distance(mode, x, y):
        forward = directed(x, y)

        if mode == "query_first":
            return forward

        reverse = directed(y, x)

        if mode == "sym_min":
            return min(forward, reverse)
        if mode == "sym_mean":
            return (forward + reverse) / 2
        if mode == "sym_max":
            return max(forward, reverse)

        raise ValueError(mode)

    # Recalculate thresholds using the exact N6 training-only
    # sampling scheme, but with the selected distance definition.
    def calibrated_threshold(training, mode):
        values = []

        for species, cultures in data.by_species.items():
            tc = cultures & training

            if len(tc) < 2:
                continue

            for culture in tc:
                for query in main_by_culture[culture]:
                    distances = []

                    for other in tc - {culture}:
                        for reference in main_by_culture[other]:
                            distances.append(
                                pair_distance(mode, query, reference)
                            )

                    if distances:
                        values.append(min(distances))

        if len(values) < 2:
            return None, len(values)

        threshold = statistics.quantiles(
            values, n=100
        )[fixture["quantile"] - 1]

        return threshold, len(values)

    # Baseline/original cache verification.
    modes = {
        "original": {
            "fires": set(),
            "distance_changed": 0,
            "nearest_changed": 0,
            "thresholds": [],
        }
    }

    for mode in MODES:
        modes[mode] = {
            "fires": set(),
            "distance_changed": 0,
            "nearest_changed": 0,
            "thresholds": [],
        }

    print("RECALIBRATING FOLD-SPECIFIC THRESHOLDS", flush=True)

    fold_thresholds = {}

    for fold_index, fold in enumerate(fixture["folds"], 1):
        target = fold["target"]
        training = frozenset(fold["training_cultures"])

        original_t, original_n = rule.threshold(training)

        fold_thresholds[(target, "original")] = (
            original_t, original_n
        )

        modes["original"]["thresholds"].append(original_t)

        for mode in MODES:
            threshold, n = calibrated_threshold(training, mode)

            if threshold is None:
                raise RuntimeError(
                    f"Missing threshold: {target}, {mode}"
                )

            fold_thresholds[(target, mode)] = (threshold, n)
            modes[mode]["thresholds"].append(threshold)

        print(
            f"[{fold_index:02}/{len(fixture['folds'])}] "
            f"{target}: original={original_t:.4f} "
            + " ".join(
                f"{m}={fold_thresholds[(target, m)][0]:.4f}"
                for m in MODES
            ),
            flush=True,
        )

    print("\nEVALUATING QUERY DECISIONS", flush=True)

    total = 0
    original_fires = set()

    for fold in fixture["folds"]:
        target = fold["target"]
        training = frozenset(fold["training_cultures"])

        original_t, _ = fold_thresholds[(target, "original")]

        alternative_distances = {
            mode: {}
            for mode in MODES
        }

        for row in fold["rows"]:
            accession = row["accession"]
            own = data.cop[accession]["culture"]
            status = row["baseline"]["status"]
            eligible = status.startswith("known")

            expected = row["expected"]
            original_fire = int(expected["its_rule_fired"])

            # Original N6 distance and threshold.
            old_candidates = {
                culture: distance
                for culture, distance
                in rule.dmin["all"].get(accession, {}).items()
                if culture in training and culture != own
            }

            if not old_candidates:
                raise RuntimeError(
                    f"No original reference: {accession}"
                )

            old_distance = min(old_candidates.values())

            if round(old_distance, 3) != float(
                expected["its_nearest_training_pct"]
            ):
                raise RuntimeError(
                    f"Original distance mismatch: {accession}"
                )

            if round(original_t, 3) != float(
                expected["its_threshold_pct"]
            ):
                raise RuntimeError(
                    f"Original threshold mismatch: {target}"
                )

            old_nearest = {
                culture for culture, distance
                in old_candidates.items()
                if abs(distance - old_distance) < 1e-12
            }

            total += 1

            if original_fire:
                original_fires.add((target, accession))
                modes["original"]["fires"].add(
                    (target, accession)
                )

            for mode in MODES:
                directed_by_culture = {}

                for culture in training:
                    if culture == own:
                        continue

                    references = all_by_culture.get(culture, [])

                    if references:
                        directed_by_culture[culture] = min(
                            pair_distance(mode, accession, reference)
                            for reference in references
                        )

                new_distance = min(directed_by_culture.values())

                new_nearest = {
                    culture for culture, distance
                    in directed_by_culture.items()
                    if abs(distance - new_distance) < 1e-12
                }

                threshold, _ = fold_thresholds[(target, mode)]
                fires = eligible and new_distance > threshold

                key = (target, accession)

                if fires:
                    modes[mode]["fires"].add(key)

                if abs(old_distance - new_distance) > 1e-9:
                    modes[mode]["distance_changed"] += 1

                if old_nearest != new_nearest:
                    modes[mode]["nearest_changed"] += 1

    print("\nFINAL COMPARISON")
    print("=" * 100)

    print(
        f"{'Method':15} {'Fires':>8} {'Removed':>9} "
        f"{'New':>8} {'Dist. changed':>15} "
        f"{'Nearest changed':>17} {'Median threshold':>18}"
    )
    print("-" * 100)

    for mode, result in modes.items():
        fires = result["fires"]

        removed = original_fires - fires
        newly_fired = fires - original_fires

        med_threshold = statistics.median(
            result["thresholds"]
        )

        print(
            f"{mode:15} {len(fires):8} "
            f"{len(removed):9} {len(newly_fired):8} "
            f"{result['distance_changed']:15} "
            f"{result['nearest_changed']:17} "
            f"{med_threshold:18.4f}"
        )

    print(f"\nTotal query-fold instances: {total}")
    print(f"Original firing count:      {len(original_fires)}")
    print(f"Unique pair distances cached: {len(pair_cache)}")

    print("\nDECISION CHANGES BY METHOD")
    for mode in MODES:
        fires = modes[mode]["fires"]
        removed = sorted(original_fires - fires)
        newly_fired = sorted(fires - original_fires)

        print(f"\n{mode}:")
        print("  Removed:")
        for target, accession in removed:
            print(f"    {target} | {accession}")

        print("  Newly fired:")
        for target, accession in newly_fired:
            print(f"    {target} | {accession}")


if __name__ == "__main__":
    main()
