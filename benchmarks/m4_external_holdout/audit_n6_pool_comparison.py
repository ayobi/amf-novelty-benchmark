#!/usr/bin/env python3

"""Compare N6 all-copy and main-only reference pools."""

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

M4 = Path(__file__).resolve().parent
PROJECTS = M4.parents[2]

sys.path.insert(0, str(M4))

import external_holdout as v1
from its_rule import ItsRule

DECISIONS = (
    PROJECTS / "amfsplit"
    / "parity/naming/n6_decision_parity/python_decisions.json"
)

OUTPUT = PROJECTS / "amf-novelty-benchmark/audit_outputs/n6_all_vs_main_pool.tsv"


def nearest(cached, training, own):
    candidates = {
        culture: distance
        for culture, distance in cached.items()
        if culture in training and culture != own
    }

    if not candidates:
        return None, []

    minimum = min(candidates.values())

    cultures = sorted(
        culture for culture, distance in candidates.items()
        if abs(distance - minimum) <= 1e-12
    )

    return minimum, cultures


def main():
    with DECISIONS.open() as fh:
        fixture = json.load(fh)

    data = v1.Data(v1.FROZEN / "inputs/benchmark")
    rule = ItsRule(data, fixture["quantile"])

    stats = Counter()
    all_fires = set()
    main_fires = set()
    changed_rows = []
    fold_stats = defaultdict(Counter)

    for fold in fixture["folds"]:
        target = fold["target"]
        training = frozenset(fold["training_cultures"])

        threshold, threshold_n = rule.threshold(training)

        for row in fold["rows"]:
            accession = row["accession"]
            own = data.cop[accession]["culture"]
            known = row["baseline"]["status"].startswith("known")
            expected = row["expected"]

            all_distance, all_nearest = nearest(
                rule.dmin["all"].get(accession, {}),
                training,
                own,
            )

            main_distance, main_nearest = nearest(
                rule.dmin["main"].get(accession, {}),
                training,
                own,
            )

            if all_distance is None or main_distance is None:
                raise RuntimeError(
                    f"Missing reference distance: {accession}"
                )

            # Validate the original N6 implementation against its fixture.
            if round(all_distance, 3) != float(
                expected["its_nearest_training_pct"]
            ):
                raise RuntimeError(
                    f"Original distance mismatch: {accession}"
                )

            if round(threshold, 3) != float(
                expected["its_threshold_pct"]
            ):
                raise RuntimeError(
                    f"Threshold mismatch: {target}"
                )

            if threshold_n != int(expected["its_threshold_n"]):
                raise RuntimeError(
                    f"Threshold sample count mismatch: {target}"
                )

            old_fire = known and all_distance > threshold
            new_fire = known and main_distance > threshold

            key = (target, accession)

            if old_fire != bool(expected["its_rule_fired"]):
                raise RuntimeError(
                    f"Original firing mismatch: {accession}"
                )

            if old_fire:
                all_fires.add(key)

            if new_fire:
                main_fires.add(key)

            stats["total"] += 1
            stats["all_fires"] += int(old_fire)
            stats["main_fires"] += int(new_fire)
            stats["distance_changed"] += int(
                abs(all_distance - main_distance) > 1e-9
            )
            stats["nearest_culture_changed"] += int(
                set(all_nearest) != set(main_nearest)
            )
            stats["all_ties"] += int(len(all_nearest) > 1)
            stats["main_ties"] += int(len(main_nearest) > 1)

            fold_stats[target]["all_fires"] += int(old_fire)
            fold_stats[target]["main_fires"] += int(new_fire)

            if old_fire and not new_fire:
                stats["removed"] += 1

            if new_fire and not old_fire:
                stats["new"] += 1

            if (
                abs(all_distance - main_distance) > 1e-9
                or set(all_nearest) != set(main_nearest)
                or old_fire != new_fire
            ):
                changed_rows.append({
                    "fold_target": target,
                    "accession": accession,
                    "source_culture": own,
                    "baseline_status": row["baseline"]["status"],
                    "baseline_placement": row["baseline"]["placement"],
                    "all_nearest_cultures": ";".join(all_nearest),
                    "main_nearest_cultures": ";".join(main_nearest),
                    "all_distance_pct": f"{all_distance:.9f}",
                    "main_distance_pct": f"{main_distance:.9f}",
                    "threshold_pct": f"{threshold:.9f}",
                    "all_fired": int(old_fire),
                    "main_fired": int(new_fire),
                })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    if changed_rows:
        with OUTPUT.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=list(changed_rows[0]),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(changed_rows)

    print("ALL vs MAIN REFERENCE POOL")
    print("=" * 65)

    for key, value in stats.items():
        print(f"{key:30} {value}")

    print("\nFIRING DECISION OVERLAP")
    print("=" * 65)
    print("Retained from all pool:", len(all_fires & main_fires))
    print("Removed in main pool:  ", len(all_fires - main_fires))
    print("New in main pool:      ", len(main_fires - all_fires))

    print("\nBY SPECIES FOLD")
    print("=" * 65)

    print(f"{'Fold target':35} {'All':>8} {'Main':>8}")
    for target, s in sorted(fold_stats.items()):
        print(
            f"{target:35} "
            f"{s['all_fires']:8} "
            f"{s['main_fires']:8}"
        )

    print("\nChanged records written to:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
