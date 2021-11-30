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


def test_group_swell_from_fasta():
    old_header, old_fields = old_swell.swell_from_fasta(fasta_path)
    new_header, new_fields = swell.group_swell_from_fasta(fasta_path)
    new_indexes = [1, 6] # positions of new data returned by swell
    for new_i in new_indexes:
        new_header = new_header[:new_i] + new_header[new_i+1:]
        new_fields = new_fields[:new_i] + new_fields[new_i+1:]
    assert new_header == old_header and new_fields == old_fields


def test_swell_from_depth():
    new_tiles = swell.load_scheme(bed_path)
    new_header, new_fields = swell.swell_from_depth(depth_path, new_tiles, ref, thresholds) 
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields


def test_swell_from_bam():
    new_tiles = swell.load_scheme(bed_path)
    new_header, new_fields = swell.swell_from_bam(bam_path, new_tiles, ref, thresholds) 
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields