'''
Created on Aug 18, 2013

@author: cambray
'''

### PATHS
output_folder = "output_20mer"  # No trailing slash
mapping_file  = "mapping.csv"
barcodes_file = "barcodes.csv"
fastq_R1      = "20140421/OLS1_S1_L001_R1"
fastq_R2      = "20140421/OLS1_S1_L001_R2"
file_ids      = [1]
db_fasta      = "combined.fasta"
water_path    = "/usr/local/bin/water"
output_mode   = "w"


### POSITIONS
bc_offset     = 0     # account for discrepancy between bc id and number of N
start_bc      = 0     # which sequence position to start looking for bc
stop_bc       = 16    # which sequence position to stop looking for (e.g. maximum spacer + length bc)
ROI_start_R1  = 23    # distance from barcode to start of region of interest
ROI_length    = 20    # length of the region of interest
ROI_start_R2  = 30    # length of NNNNNNGTTGCCTTCATACTGTCGAGCCAT
ROI_offset    =  0    # extra nt to include to account for potential indels (adds this amount on both sides)


### OTHERS
aln_threshold = 0.8
make_index    = True
verbose       = True

def find_match (seq, table, possible_starts):
    for k in possible_starts:
        if seq[k:k+20] in table:
            return seq[k:k+20]
    return None