#!/usr/bin/env python3

"""Reconstruct nearest ITS training references for every N6 firing."""

import argparse
import csv
import json
import sys
from collections import Counter

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M4 = Path(__file__).resolve().parent

# Import the frozen M2b data layer and the original N6 implementation.
sys.path.insert(0, str(M4))

import external_holdout as v1
from its_rule import ItsRule, dist as its_dist


def split_placement(value):
    return {x.strip() for x in (value or "").split(";") if x.strip()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--decisions",
        type=Path,
        default=(
            ROOT.parent
            / "amfsplit"
            / "parity/naming/n6_decision_parity/python_decisions.json"
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "audit_outputs/n6_nearest_references.tsv",
    )
    args = parser.parse_args()

    print("Loading frozen reference data...", flush=True)
    data = v1.Data(v1.FROZEN / "inputs/benchmark")

    with args.decisions.open() as fh:
        decisions = json.load(fh)

    if decisions["pool"] != "all":
        raise ValueError("Expected the frozen N6 all-copy pool")

    quantile = decisions["quantile"]

    print(
        f"Reference copies: {len(data.cop)} | "
        f"Cultures: {len(data.cultures)} | "
        f"ITS sequences: loading N6 cache...",
        flush=True,
    )

    rule = ItsRule(data, quantile)

    # Preserve the accession ordering used by the original N6 cache.
    pair_order = {
        accession: i
        for i, accession in enumerate(rule.its_ok)
    }

    print("Reconstructing firing decisions...", flush=True)

    output = []
    checks = Counter()
    placement_comparison = Counter()
    unique_nearest_cultures = set()

    for fold in decisions["folds"]:
        target = fold["target"]
        training = frozenset(fold["training_cultures"])

        threshold, threshold_n = rule.threshold(training)

        for row in fold["rows"]:
            expected = row["expected"]

            if int(expected["its_rule_fired"]) != 1:
                continue

            accession = row["accession"]
            own_culture = data.cop[accession]["culture"]
            query_sequence = data.fa["full_its"][accession]

            # Preserve the exact N6 training-culture selection.
            candidates = [
                (culture, distance)
                for culture, distance
                in rule.dmin["all"].get(accession, {}).items()
                if culture in training and culture != own_culture
            ]

            if not candidates:
                raise RuntimeError(
                    f"No eligible training reference for {accession}"
                )

            nearest_distance = min(distance for _, distance in candidates)

            nearest_cultures = sorted(
                culture
                for culture, distance in candidates
                if abs(distance - nearest_distance) <= 1e-12
            )

            # Verify distance, threshold and sample count against
            # the frozen Python N6 decision fixture.
            expected_distance = float(
                expected["its_nearest_training_pct"]
            )
            expected_threshold = float(expected["its_threshold_pct"])

            distance_ok = round(nearest_distance, 3) == expected_distance
            threshold_ok = (
                round(threshold, 3) == expected_threshold
                and threshold_n == int(expected["its_threshold_n"])
            )

            checks["records"] += 1
            checks["distance_matches"] += int(distance_ok)
            checks["threshold_matches"] += int(threshold_ok)

            if not distance_ok:
                raise RuntimeError(
                    f"Distance mismatch for {accession}: "
                    f"{nearest_distance} vs {expected_distance}"
                )

            if not threshold_ok:
                raise RuntimeError(
                    f"Threshold mismatch for {target}: "
                    f"{threshold} vs {expected_threshold}"
                )

            baseline = row["baseline"]
            placement = baseline["placement"]
            placement_species = split_placement(placement)

            source_species = data.species.get(own_culture, "")

            nearest_reference_records = []

            for culture in nearest_cultures:
                copy_ids = [
                    y for y in data.by_culture[culture]
                    if y in rule.its_ok
                ]

                # Recover the actual accession(s) supplying the
                # minimum distance within each nearest culture.
                copy_distances = []

                for y in copy_ids:
                    reference_sequence = data.fa["full_its"][y]

                    # Reproduce the exact accession ordering used by ItsRule.__init__.
                    if pair_order[accession] < pair_order[y]:
                        d = its_dist(query_sequence, reference_sequence)
                    else:
                        d = its_dist(reference_sequence, query_sequence)

                    copy_distances.append((d, y))

                culture_minimum = min(d for d, _ in copy_distances)

                nearest_accessions = sorted(
                    accession_id
                    for distance, accession_id in copy_distances
                    if abs(distance - culture_minimum) <= 1e-12
                )

                if abs(culture_minimum - nearest_distance) > 1e-9:
                    raise RuntimeError(
                        f"Accession-level distance mismatch: "
                        f"{accession} / {culture}"
                    )

                reference_species = data.species.get(culture, "")
                unique_nearest_cultures.add(culture)

                nearest_reference_records.append({
                    "culture": culture,
                    "species": reference_species,
                    "accessions": nearest_accessions,
                    "copy_count": len(copy_ids),
                })

            nearest_species = {
                x["species"]
                for x in nearest_reference_records
                if x["species"]
            }

            placement_matches_reference = bool(
                placement_species & nearest_species
            )
            reference_matches_source = bool(
                source_species and source_species in nearest_species
            )
            placement_matches_source = bool(
                source_species and source_species in placement_species
            )

            placement_comparison[
                "baseline matches nearest reference"
                if placement_matches_reference else
                "baseline differs from nearest reference"
            ] += 1

            placement_comparison[
                "nearest reference matches source species"
                if reference_matches_source else
                "nearest reference differs from source species"
            ] += 1

            output.append({
                "fold_target": target,
                "accession": accession,
                "source_culture": own_culture,
                "source_species": source_species,
                "baseline_status": baseline["status"],
                "baseline_placement": placement,
                "nearest_cultures": ";".join(nearest_cultures),
                "nearest_reference_species": ";".join(sorted(nearest_species)),
                "nearest_reference_accessions": ";".join(
                    f"{ref['culture']}:[{','.join(ref['accessions'])}]"
                    for ref in nearest_reference_records
                ),
                "nearest_culture_copy_counts": ";".join(
                    f"{ref['culture']}:{ref['copy_count']}"
                    for ref in nearest_reference_records
                ),
                "nearest_distance_pct": f"{nearest_distance:.9f}",
                "threshold_pct": f"{threshold:.9f}",
                "excess_pct": f"{nearest_distance - threshold:.9f}",
                "baseline_matches_nearest_species": int(
                    placement_matches_reference
                ),
                "nearest_reference_matches_source_species": int(
                    reference_matches_source
                ),
                "baseline_matches_source_species": int(
                    placement_matches_source
                ),
            })

    if checks["records"] != decisions["total_fired"]:
        raise RuntimeError(
            f"Expected {decisions['total_fired']} firings, "
            f"reconstructed {checks['records']}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)

    fields = list(output[0].keys())

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output)

    print("\nRECONSTRUCTION VALIDATION")
    print("=" * 55)
    print(f"Firing records:             {checks['records']}")
    print(f"Distance matches:           {checks['distance_matches']}")
    print(f"Threshold matches:          {checks['threshold_matches']}")
    print(f"Unique nearest cultures:    {len(unique_nearest_cultures)}")

    print("\nREFERENCE IDENTITY COMPARISON")
    print("=" * 55)

    for label, count in placement_comparison.most_common():
        print(f"{label:48} {count:4}")

    print(f"\nAudit written to: {args.out}")


if __name__ == "__main__":
    main()
