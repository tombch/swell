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


def test_swell_from_fasta():
    old_header, old_fields = old_swell.swell_from_fasta(fasta_path)
    new_header, new_fields = swell.swell_from_fasta(fasta_path)
    ambiguous_data_index = 6 # position of new data (pc_ambiguous) returned by swell
    new_header = new_header[:ambiguous_data_index] + new_header[ambiguous_data_index+1:] # remove pc_ambiguous header
    new_fields = new_fields[:ambiguous_data_index] + new_fields[ambiguous_data_index+1:] # remove pc_ambiguous field
    assert new_header == old_header and new_fields == old_fields


def test_swell_from_depth():
    new_tiles = swell.load_scheme(bed_path)
    new_tile_vector, new_header, new_fields = swell.swell_from_depth(depth_path, new_tiles, ref, thresholds) 
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields


def test_swell_from_bam():
    new_tiles = swell.load_scheme(bed_path)
    new_tile_vector, new_header, new_fields = swell.swell_from_bam(bam_path, new_tiles, ref, thresholds) 
    old_tiles = old_swell.load_scheme(bed_path)
    old_header, old_fields = old_swell.swell_from_depth(depth_path, old_tiles, ref, thresholds)
    assert new_header == old_header and new_fields == old_fields
