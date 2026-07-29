from variant_triage.fastq import count_records

def test_count_records(tmp_path):
    f = tmp_path / "tiny.fastq"
    f.write_text("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n")
    assert count_records(f) == 2