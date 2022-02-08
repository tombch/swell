from swell import swell
from swell import old_swell

# python3 -m pytest -v

bed_path = '../samples/nCoV-2019.scheme.bed'
fasta_path = '../samples/s1.fasta'
depth_path = '../samples/s1.bam.depth'
bam_path = '../samples/s1.bam'
ref = 'MN908947.3'
thresholds = [1, 5, 10, 20, 50, 100, 200]


def test_load_scheme():
    assert swell.load_scheme(bed_path) == old_swell.load_scheme(bed_path)


def test_summarise_swell_from_fasta():
    old_header, old_fields = old_swell.swell_from_fasta(fasta_path)
    new_header, new_fields = swell.summarise_swell_from_fasta(fasta_path)
    new_fields = new_fields[0]
    new_indexes = [1, 7] # Positions of new data returned by swell (header and pc_ambiguous columns)
    matching_new_header = []
    matching_new_fields = []
    for i, (header, field) in enumerate(zip(new_header, new_fields)):
        if not (i in new_indexes):
            matching_new_header.append(header)
            matching_new_fields.append(field)
    assert matching_new_header == old_header and matching_new_fields == old_fields


def test_swell_from_fasta():
    old_header, old_fields = old_swell.swell_from_fasta(fasta_path)
    new_header, new_fields = swell.swell_from_fasta(fasta_path)
    new_fields = new_fields[0]
    new_indexes = [1, 7] # Positions of new data returned by swell (header and pc_ambiguous columns)
    matching_new_header = []
    matching_new_fields = []
    for i, (header, field) in enumerate(zip(new_header, new_fields)):
        if not (i in new_indexes):
            matching_new_header.append(header)
            matching_new_fields.append(field)
    assert matching_new_header == old_header and matching_new_fields == old_fields


def test_swell_from_depth():
    new_tiles = swell.load_scheme(bed_path)
    new_header, new_fields = swell.swell_from_depth(depth_path, new_tiles, ref, thresholds)
    new_fields = new_fields[0]
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields


def test_swell_from_bam():
    new_tiles = swell.load_scheme(bed_path)
    new_header, new_fields = swell.swell_from_bam(bam_path, new_tiles, ref, thresholds)
    new_fields = new_fields[0]
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields