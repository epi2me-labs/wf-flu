"""Tests for report helpers."""

from workflow_glue import report


def test_selected_panel_label_prefers_specific_panel_call():
    """Panel labels should use the most specific resolved call available."""
    assert report.selected_panel_label(
        {"type": "Type_B", "HA": "Victoria", "NA": "undetermined"}
    ) == "B / Victoria"
    assert report.selected_panel_label(
        {"type": "Type_A", "HA": "H3", "NA": "N1"}
    ) == "A / H3N1"
    assert report.selected_panel_label(
        {"type": "Type_A", "HA": "undetermined", "NA": "undetermined"}
    ) == "A"


def test_mapped_read_counts_reads_count_files(tmp_path):
    """Mapped-read helper should assemble a report dataframe from count files."""
    sample_files = {
        "sample-a": {"mapped_reads": str(tmp_path / "sample-a.txt")},
        "sample-b": {"mapped_reads": str(tmp_path / "sample-b.txt")},
    }
    (tmp_path / "sample-a.txt").write_text("12\n")
    (tmp_path / "sample-b.txt").write_text("7\n")
    sample_details = [
        {"sample": "sample-a", "type": "Type_B"},
        {"sample": "sample-b", "type": "Type_A"},
    ]

    df = report.mapped_read_counts(sample_details, sample_files)

    assert df.to_dict("records") == [
        {"Sample": "sample-a", "Mapped reads": 12, "Selected panel": "B"},
        {"Sample": "sample-b", "Mapped reads": 7, "Selected panel": "A"},
    ]
