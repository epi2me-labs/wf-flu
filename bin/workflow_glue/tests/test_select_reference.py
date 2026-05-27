"""Tests for reference-panel selection helpers."""
import pysam
from workflow_glue.select_reference import choose_reference_headers, write_reference
from workflow_glue.typing_resolver import load_typing_table


def test_load_typing_table_handles_empty_file(tmp_path):
    """An empty typing file should resolve to an empty dataframe, not a crash."""
    typing_file = tmp_path / "empty.tsv"
    typing_file.write_text("")

    df = load_typing_table(typing_file)

    assert df.empty


def test_choose_reference_headers_keeps_full_panel_when_type_is_undetermined():
    """Undetermined first-pass typing should fall back to the full panel."""
    headers = ["A_MP", "A_HA_H1", "A_NA_N1", "B_HA", "B_NA"]

    selected = choose_reference_headers(
        {
            "type": "undetermined",
            "HA": "undetermined",
            "NA": "undetermined",
            "lineage": "undetermined",
        },
        headers,
    )

    assert selected == headers


def test_choose_reference_headers_falls_back_to_all_ha_and_na_when_missing():
    """Missing subtype references should widen selection instead of dropping them."""
    headers = ["A_MP", "A_HA_H1", "A_HA_H3", "A_NA_N1", "A_NA_N2"]

    selected = choose_reference_headers(
        {"type": "Type_A", "HA": "H5", "NA": "N7", "lineage": "undetermined"},
        headers,
    )

    assert selected == ["A_MP", "A_HA_H1", "A_HA_H3", "A_NA_N1", "A_NA_N2"]


def test_write_reference_emits_only_selected_records(tmp_path):
    """Selected reference FASTA should contain only the requested headers."""
    reference = tmp_path / "reference.fasta"
    reference.write_text(
        ">A_MP\nACGTACGT\n"
        ">A_HA_H1\nTTTTGGGG\n"
        ">B_HA\nCCCCAAAA\n"
    )
    pysam.faidx(str(reference))

    output = tmp_path / "selected.fasta"
    write_reference(reference, ["A_MP", "B_HA"], output)

    assert output.read_text() == ">A_MP\nACGTACGT\n>B_HA\nCCCCAAAA\n"
