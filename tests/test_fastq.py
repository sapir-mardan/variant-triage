import pytest
import pandas as pd
from variant_triage.fastq import count_records, validate_record, qc_summary, read_stats

def test_count_records(tmp_path):
    f = tmp_path / "tiny.fastq"
    f.write_text("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n")
    assert count_records(f) == 2

def test_valid_record():
    #the test's name should describe the scenario being tested
    lines = ["@ssr-name\n", "AGTA\n", "+\n", "@II-\n"]
    assert validate_record(lines) == True

@pytest.fixture
def valid_record():
    return(["@ssr-name\n", "AGTA\n", "+\n", "@II-\n"])

def test_valid_record_with_fixture(valid_record): #in pytest every parameter is a call for a fixture hence one should be defined before.
    assert validate_record(valid_record) == True

@pytest.mark.parametrize("lines, reason", [
    (["Ssr-name\n", "AGTA\n", "+\n", "@II-\n"], 'bad header'),
    (["@ssr-name\n", "AGTA\n", "@\n", "@II-\n"], "bad plus sign"),
    (["@ssr-name\n", "AGT\n", "+\n", "@II-\n"], "dna and score not same lenght")
])

def test_non_valid_records(lines, reason):
    assert validate_record(lines) == False

def test_count_records_raises_exception_on_files_with_truncated_records(tmp_path): # must be exactly "tmp_path"
    bad_file = tmp_path / "trancuated.fastq"
    bad_file.write_text("@read1\nACGT\n+\n") #creating a file with 3 lined record

    with pytest.raises(ValueError):
        count_records(bad_file)


@pytest.mark.parametrize("seq, expected_gc", [
    ("GGCC", 100.0),
    ("ggcc", 100.0),   # or just write "ggcc"
    ("GgCc", 100.0),
    ("AATT", 0.0),
    ("", None)
    ])

def test_gc_content_case_insensitive(tmp_path, seq, expected_gc):
    f = tmp_path / "test.fastq"
    f.write_text(f"@r1\n{seq}\n+\n{'I' * len(seq)}\n")
    result = qc_summary(f)
    assert result["gc_content"] == expected_gc

def test_read_stats_runs_with_empty_file(tmp_path):
    f = tmp_path / "empty.fastq"
    f.write_text("")
    df = read_stats(f)
    assert len(df) == 0