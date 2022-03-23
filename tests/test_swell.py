import glob
from swell import swell
from tests import swell_original

# Put paths here
bed_new = ""
bed_old = ""
samples_dir = ""

fastas = glob.glob(f'{samples_dir}/*.fasta')
bams = glob.glob(f'{samples_dir}/*.bam')
depths = [x + '.depth' for x in bams]
genomes = ['NC_045512', 'NC045512', 'MN908947.3']
thresholds = [1, 5, 10, 20, 50, 100, 200]


def test_load_scheme():
    assert swell.load_scheme(bed_new) == swell_original.load_scheme(bed_old)


def test_swell_from_fasta():
    for fasta in fastas:
        old_header, old_fields = swell_original.swell_from_fasta(fasta)
        new_header, new_fields = swell.swell_from_fasta(fasta)
        new_fields = new_fields[0]
        new_headers = ['header', 'pc_ambiguous']
        matching_new_header = []
        matching_new_fields = []
        for header, field in zip(new_header, new_fields):
            if not (header in new_headers):
                matching_new_header.append(header)
                matching_new_fields.append(field)
        assert matching_new_header == old_header and matching_new_fields == old_fields


def test_swell_from_fasta_summarise():
    for fasta in fastas:
        old_header, old_fields = swell_original.swell_from_fasta(fasta)
        new_header, new_fields = swell.summarise_swell_from_fasta(fasta)
        new_fields = new_fields[0]
        new_headers = ['header', 'pc_ambiguous']
        matching_new_header = []
        matching_new_fields = []
        for header, field in zip(new_header, new_fields):
            if not (header in new_headers):
                matching_new_header.append(header)
                matching_new_fields.append(field)
        assert matching_new_header == old_header and matching_new_fields == old_fields


def test_swell_from_bam():
    new_tiles = swell.load_scheme(bed_new)
    old_tiles = swell_original.load_scheme(bed_old)
    for bam, depth in zip(bams, depths):
        new_header, new_fields = swell.swell_from_bam(bam, new_tiles, genomes, thresholds)
        new_fields = new_fields[0]
        old_header, old_fields = swell_original.swell_from_depth(depth, old_tiles, genomes, thresholds)
        assert new_header == old_header and new_fields == old_fields


def test_swell_from_bam_no_bed():
    for bam, depth in zip(bams, depths):
        new_header, new_fields = swell.swell_from_bam(bam, {}, genomes, thresholds)
        new_fields = new_fields[0]
        old_header, old_fields = swell_original.swell_from_depth(depth, {}, genomes, thresholds)
        assert new_header == old_header and new_fields == old_fields


def test_swell_from_depth():
    new_tiles = swell.load_scheme(bed_new)
    old_tiles = swell_original.load_scheme(bed_old)
    for depth in depths:
        new_header, new_fields = swell.swell_from_depth(depth, new_tiles, genomes, thresholds)
        new_fields = new_fields[0]
        old_header, old_fields = swell_original.swell_from_depth(depth, old_tiles, genomes, thresholds)
        assert new_header == old_header and new_fields == old_fields


def test_swell_from_depth_no_bed():
    for depth in depths:
        new_header, new_fields = swell.swell_from_depth(depth, {}, genomes, thresholds)
        new_fields = new_fields[0]
        old_header, old_fields = swell_original.swell_from_depth(depth, {}, genomes, thresholds)
        assert new_header == old_header and new_fields == old_fields
