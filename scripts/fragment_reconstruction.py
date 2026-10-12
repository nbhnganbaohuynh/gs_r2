import sys
import pysam
import itertools


bam_path, refseq_path, fragment_path = sys.argv[1:4]

# fragment rules
MIN_MAPQ = 30
MIN_LENGTH = 20
MAX_LENGTH = 600

# reverse complement
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


# validate motif
MOTIF_LENGTH = 4

def valid_motif(motif):
    return len(motif) == MOTIF_LENGTH and all(base in "ACGT" for base in motif)


with pysam.AlignmentFile(bam_path, "rb") as bam_file, \
     pysam.FastaFile(refseq_path) as reference_file, \
     open(fragment_path, "w") as fragment_file:
    fragment_file.write("chrom\tlength\tleft_motif\tright_motif\n")
    
    for _, group in itertools.groupby(bam_file, key=lambda r: r.query_name):
        # filter reads
        reads = [
            r for r in group
            if r.is_paired
            and r.is_proper_pair
            and not r.is_unmapped
            and not r.mate_is_unmapped
            and not r.is_secondary
            and not r.is_supplementary
            and not r.is_duplicate
            and not r.is_qcfail
            and r.mapping_quality >= MIN_MAPQ
        ]
        
        r1_list = [r for r in reads if r.is_read1]
        r2_list = [r for r in reads if r.is_read2]
        
        # check read pairs
        if len(r1_list) != 1 or len(r2_list) != 1:
            continue
        
        r1, r2 = r1_list[0], r2_list[0]
        
        if r1.reference_id != r2.reference_id:
            continue
        if r1.reference_end is None or r2.reference_end is None:
            continue
        
        # chromosome name
        chrom = r1.reference_name
        
        # fragment length
        start = min(r1.reference_start, r2.reference_start)
        end = max(r1.reference_end, r2.reference_end)
        length = end - start
        
        if not MIN_LENGTH <= length <= MAX_LENGTH:
            continue
        
        # end motif
        left_motif = reference_file.fetch(chrom, start, start + MOTIF_LENGTH).upper()
        right_motif = reverse_complement(reference_file.fetch(chrom, end - MOTIF_LENGTH, end).upper())
        
        if not valid_motif(left_motif):
            left_motif = "NA"
        
        if not valid_motif(right_motif):
            right_motif = "NA"
        
        # write fragment information
        fragment_file.write(f"{chrom}\t{length}\t{left_motif}\t{right_motif}\n")
