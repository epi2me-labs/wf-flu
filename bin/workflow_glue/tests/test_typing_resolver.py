"""Tests for influenza typing resolution."""
import pandas as pd
from workflow_glue.select_reference import choose_reference_headers
from workflow_glue.typing_resolver import as_processed_json, resolve_typing_calls


def make_typing_df(rows):
    """Build a scored Abricate-like table for unit tests."""
    df = pd.DataFrame(rows)
    df["_score"] = df["%COVERAGE"] + df["%IDENTITY"]
    return df.sort_values(
        by=["_score", "%COVERAGE", "%IDENTITY"], ascending=False
    ).reset_index(drop=True)


def test_resolve_type_b_when_lineage_outweighs_lone_a_matrix_hit():
    """B lineage and B matrix support should beat a lone A matrix cross-hit."""
    df = make_typing_df([
        {
            "GENE": "M1",
            "RESISTANCE": "Type_A",
            "%COVERAGE": 100.00,
            "%IDENTITY": 97.86,
        },
        {
            "GENE": "HA",
            "RESISTANCE": "Victoria",
            "%COVERAGE": 95.83,
            "%IDENTITY": 93.11,
        },
        {"GENE": "M1", "RESISTANCE": "Type_B", "%COVERAGE": 99.82, "%IDENTITY": 95.97},
    ])

    calls = resolve_typing_calls(df)

    assert calls == {
        "type": "Type_B",
        "HA": "Victoria",
        "NA": "undetermined",
        "lineage": "Victoria",
    }
    assert as_processed_json(calls) == {
        "HA": ["Victoria"],
        "NA": ["undetermined"],
        "type": ["Type_B"],
    }


def test_resolve_type_a_ignores_b_cross_hits_when_a_has_ha_and_na_support():
    """A sample with clear H/N support should stay A despite B cross-hits."""
    df = make_typing_df([
        {"GENE": "HA", "RESISTANCE": "H3", "%COVERAGE": 99.94, "%IDENTITY": 97.71},
        {"GENE": "M1", "RESISTANCE": "Type_A", "%COVERAGE": 99.90, "%IDENTITY": 98.58},
        {"GENE": "NA", "RESISTANCE": "N2", "%COVERAGE": 99.86, "%IDENTITY": 98.51},
        {
            "GENE": "HA",
            "RESISTANCE": "Yamagata",
            "%COVERAGE": 95.06,
            "%IDENTITY": 93.89,
        },
        {"GENE": "M1", "RESISTANCE": "Type_B", "%COVERAGE": 99.82, "%IDENTITY": 95.72},
    ])

    calls = resolve_typing_calls(df)

    assert calls == {
        "type": "Type_A",
        "HA": "H3",
        "NA": "N2",
        "lineage": "undetermined",
    }


def test_resolve_type_a_canonicalizes_variant_na_calls():
    """Variant subtype labels such as N3v2 should be normalized for reporting."""
    df = make_typing_df([
        {"GENE": "HA", "RESISTANCE": "H7", "%COVERAGE": 98.0, "%IDENTITY": 97.0},
        {"GENE": "M1", "RESISTANCE": "Type_A", "%COVERAGE": 99.0, "%IDENTITY": 98.0},
        {"GENE": "NA", "RESISTANCE": "N3v2", "%COVERAGE": 96.35, "%IDENTITY": 91.43},
    ])

    calls = resolve_typing_calls(df)

    assert calls["type"] == "Type_A"
    assert calls["HA"] == "H7"
    assert calls["NA"] == "N3"


def test_choose_reference_headers_limits_a_panel_to_selected_ha_and_na():
    """Reference selection should keep A internals plus the chosen HA/NA."""
    headers = [
        "A_MP",
        "A_NP",
        "A_HA_H1",
        "A_HA_H3",
        "A_NA_N1",
        "A_NA_N2",
        "B_HA",
        "B_MP",
    ]

    selected = choose_reference_headers(
        {"type": "Type_A", "HA": "H1", "NA": "N1", "lineage": "undetermined"},
        headers,
    )

    assert selected == ["A_MP", "A_NP", "A_HA_H1", "A_NA_N1"]


def test_choose_reference_headers_returns_all_b_segments_for_type_b():
    """Reference selection should use the B segment panel for influenza B."""
    headers = [
        "A_MP",
        "A_HA_H1",
        "A_NA_N1",
        "B_HA",
        "B_MP",
        "B_NA",
    ]

    selected = choose_reference_headers(
        {
            "type": "Type_B",
            "HA": "Yamagata",
            "NA": "undetermined",
            "lineage": "Yamagata",
        },
        headers,
    )

    assert selected == ["B_HA", "B_MP", "B_NA"]
