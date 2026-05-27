#!/usr/bin/env python
"""Read Abricate flu hits and decide the one call the workflow should use.

This file reads the competing A and B hits from Abricate, works out which type
is best supported, then turns that into one simple result for the rest of the
workflow. That same result is used both for picking the second-pass reference
panel and for building the final typing report.
"""

import re

import pandas as pd
from pandas.errors import EmptyDataError


TYPE_A = "Type_A"
TYPE_B = "Type_B"
LINEAGES = {"Victoria", "Yamagata"}


def load_typing_table(path):
    """Load and score an Abricate typing table."""
    try:
        abricate = pd.read_csv(path, delimiter="\t", keep_default_na=False)
    except EmptyDataError:
        return pd.DataFrame()

    if abricate.empty:
        return abricate

    for column in ("%COVERAGE", "%IDENTITY"):
        abricate[column] = pd.to_numeric(abricate[column], errors="coerce").fillna(0)

    # Use a simple additive score so we can compare competing hits consistently
    # across type-, subtype-, and lineage-level matches.
    abricate["_score"] = abricate["%COVERAGE"] + abricate["%IDENTITY"]

    return abricate.sort_values(
        by=["_score", "%COVERAGE", "%IDENTITY"], ascending=False
    ).reset_index(drop=True)


def best_row(df, gene=None, allowed=None, pattern=None):
    """Return the strongest Abricate hit for a category."""
    if df.empty:
        return None

    subset = df
    if gene is not None:
        subset = subset[subset["GENE"] == gene]
    if allowed is not None:
        subset = subset[subset["RESISTANCE"].isin(allowed)]
    if pattern is not None:
        subset = subset[subset["RESISTANCE"].str.match(pattern)]

    if subset.empty:
        return None
    return subset.iloc[0]


def row_score(row):
    """Return the numeric ranking score for a row."""
    if row is None:
        return 0.0
    return float(row["_score"])


def canonical_subtype(value, prefix):
    """Reduce subtype variants like N3v1 to their canonical form."""
    if not value:
        return None
    match = re.match(rf"^({prefix}\d+)", value)
    if match:
        return match.group(1)
    return value


def choose_type(df):
    """Resolve the sample type from competing A/B signals.

    Abricate cross-hits are common for influenza. We therefore rank type
    evidence using the number of independent signals first, then total score.
    This prevents a lone A matrix hit from overriding stronger B matrix+lineage
    support, while still allowing A samples with clear HA/NA support to win.
    """
    a_type = best_row(df, gene="M1", allowed=[TYPE_A])
    a_ha = best_row(df, gene="HA", pattern=r"^H\d+")
    a_na = best_row(df, gene="NA", pattern=r"^N\d+")
    b_type = best_row(df, gene="M1", allowed=[TYPE_B])
    b_lineage = best_row(df, gene="HA", allowed=sorted(LINEAGES))

    a_hits = [row for row in (a_type, a_ha, a_na) if row is not None]
    b_hits = [row for row in (b_type, b_lineage) if row is not None]

    a_score = sum(row_score(row) for row in a_hits)
    b_score = sum(row_score(row) for row in b_hits)

    if b_lineage is not None and len(b_hits) > len(a_hits):
        return TYPE_B
    if len(a_hits) > len(b_hits):
        return TYPE_A
    if b_lineage is not None and b_score > a_score:
        return TYPE_B
    if a_score > b_score:
        return TYPE_A
    if b_lineage is not None:
        return TYPE_B
    if a_hits:
        return TYPE_A
    if b_hits:
        return TYPE_B
    return "undetermined"


def resolve_typing_calls(df):
    """Resolve a single type/subtype call from an Abricate typing table."""
    if df.empty:
        return {
            "type": "undetermined",
            "HA": "undetermined",
            "NA": "undetermined",
            "lineage": "undetermined",
        }

    flu_type = choose_type(df)
    if flu_type == TYPE_A:
        ha_row = best_row(df, gene="HA", pattern=r"^H\d+")
        na_row = best_row(df, gene="NA", pattern=r"^N\d+")
        return {
            "type": TYPE_A,
            "HA": canonical_subtype(
                ha_row["RESISTANCE"] if ha_row is not None else None, "H"
            ) or "undetermined",
            "NA": canonical_subtype(
                na_row["RESISTANCE"] if na_row is not None else None, "N"
            ) or "undetermined",
            "lineage": "undetermined",
        }

    if flu_type == TYPE_B:
        lineage_row = best_row(df, gene="HA", allowed=sorted(LINEAGES))
        lineage = lineage_row["RESISTANCE"] if lineage_row is not None else None
        return {
            "type": TYPE_B,
            "HA": lineage or "undetermined",
            "NA": "undetermined",
            "lineage": lineage or "undetermined",
        }

    return {
        "type": "undetermined",
        "HA": "undetermined",
        "NA": "undetermined",
        "lineage": "undetermined",
    }


def as_processed_json(calls):
    """Render resolved typing calls using the legacy JSON structure."""
    return {
        "HA": [calls["HA"]],
        "NA": [calls["NA"]],
        "type": [calls["type"]],
    }
