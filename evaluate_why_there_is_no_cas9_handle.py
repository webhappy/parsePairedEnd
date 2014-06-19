from subprocess import call
from Bio import SeqIO
import sys
import csv
from Bio import pairwise2
from parameters import *
from call_bowtie_to_align import *

#MAPPINGS_TO_EVALUATE = ['pT461']  # ['pT461','pT462']

CAS9_HANDLE = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAAC'
R1_EXPECTED_PROMOTER = 'gcaactctctactgtttctccat'.upper()
R2_EXPECTED_PROMOTER = 'GTTGCCTTCATACTGTCGAGCCAT'
POST_sgRNA_SEQ_1 = 'AAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGAAGCTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGACTCCAGCTTGGCTGTTTTGGCGGATGAGAGAAGATTTTCAGCCTGATACAGATTAAATCAGAACGCAGAAGCGGTCTGATAAAACAGAATTTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGACTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACTAGATTCAGATCCTCTTCTGAGATGAGTTTTTGTTCGGGCCCAAGCTTCAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTT'
POST_sgRNA_SEQ_2 = 'AAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGAAGCTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATCTAGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTGGCAGTTCCCTACTCTCGCATGGGGAGTCCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACTGCCGCCAGGCAAATTCTGTTTTATCAGACCGCTTCTGCGTTCTGATTTAATCTGTATCAGGCTGAAAATCTTCTCTCATCCGCCAAAACAGCCAAGCTGGAGTCCGTTTAAACTCAATGATGATGATGATGATGGTCGACGGCGCTATTCAGATCCTCTTCTGAGATGAGTTTTTGTTCGGGCCCAAGCTTCAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTT'

#### FUNCTIONS
def find_degenerate_barcode(seq, orientation="r1", offset=-1):
    """
    find degenerate match where expected
    """

    # where expected
    start = offset
    stop  = MAX_BARCODE_LENGTH + offset
    cand = str(seq[start:stop])
    if cand in degenerate_barcodes[orientation] :
        bc_id = degenerate_barcodes[orientation][cand]
        bc_len = MAX_BARCODE_LENGTH - BARCODE_PAD_LENGTH[bc_id]
        return bc_id, start+bc_len, (bc_len-1)/float(bc_len)

    return -1, -1, -1

def water_bc(seq, orientation, id, threshold=0.8):
    """
    Align barcodes to sequence and return higest above threshold
    """
    match = {}
    # write sequence to temp file
    h=open("%s/temp%i" % (output_folder, id), "w")
    h.write(str(seq))
    h.close()
    
    # call alignment and read output file
    water_params = '-sformat2 fasta -aformat fasta -awidth3 200 -auto -snucleotide1 -snucleotide2'
    call('%s "%s/temp%i" "%s/bc_%s.fas" -gapopen 4 -gapextend 10 "%s/temp%i.aln" %s'
         % (water_path, output_folder, id, output_folder, orientation, output_folder, id, water_params), shell=True)
    h=open("%s/temp%i.aln" % (output_folder, id), "r")
    lines = h.readlines()
    h.close()
    
    #parse alignments and get scores 
    for i in range(0,len(lines)-1,4):
        if bc_len-2 <= len(lines[i+1]) <= bc_len+2:
            score = 0.
            for j in range(len(lines[i+1])):
                if lines[i+1][j] == lines[i+3][j]:
                    score+=1
            id = lines[i+2].strip()[1:]
            score /= max(bc_len, len(lines[i+1]))
            match[score] = (int(lines[i+2][1]), lines[i+1], lines[i+3], score)
    if match:
        hi_score = max(match.keys())  
        motif = match[hi_score][1].strip().replace("-","")
        pos   = seq.find(motif)+len(motif)
        if hi_score >= threshold:
            return match[hi_score][0], pos, hi_score
    return -1, -1, -1

def printPairedStatus():
    print """
        # %i
        Promoter_matches_r1 = %.2f%%
        Promoter_matches_r2 = %.2f%%
        Handle_matches = %.2f%%
        %% per bc  = %s
        """ % (total_count, r1_promoter_matches_count/float(total_count)*100, r2_promoter_matches_count/float(total_count)*100,
               handle_matches_count/float(total_count)*100, " | ".join(
        ["%s: %.2f" % (key, readperbc[key] / float(total_count) * 100) for key in readperbc if readperbc[key] != 0]))

def printUnpairedStatus():
    print """
        # %i
        Promoter_matches = %.2f%%
        Handle_matches = %.2f%%
        %% per bc  = %s
        """ % (total_count, promoter_matches_count/float(total_count)*100, handle_matches_count/float(total_count)*100, " | ".join(
        ["%s: %.2f" % (key, readperbc[key] / float(total_count) * 100) for key in readperbc if readperbc[key] != 0]))

MAX_BARCODE_LENGTH = 0   # Initialize while making the hash table
BARCODE_PAD_LENGTH = {}

def recurseAllCombinations (orig, cur, storeTable, valToAssign, deviations_allowed = 1):
    if len(cur) == len(orig):  # Base this to end recursion
        if cur in storeTable and valToAssign != storeTable[cur]:
            print cur,'is already in the hash table, there is a bug here'
            raise
        else:
            storeTable[cur] = valToAssign
            return

    LETTERS = ['A', 'T', 'G', 'C', 'N']

    cur_idx = len(cur)
    if orig[cur_idx] == 'N':
        for k in xrange(0,len(LETTERS)):
            recurseAllCombinations(orig, cur+LETTERS[k], storeTable, valToAssign, deviations_allowed)
    else:
        recurseAllCombinations(orig, cur+orig[cur_idx], storeTable, valToAssign, deviations_allowed)

        if deviations_allowed > 0:
            for k in xrange(0,len(LETTERS)):
                if LETTERS[k] == orig[cur_idx]:
                    continue
                recurseAllCombinations(orig, cur+LETTERS[k], storeTable, valToAssign, deviations_allowed-1)

def generateBarcodeHashTable (candidates):
    """

    :param candidates: Hashmap of barcode string to ID
    """
    global MAX_BARCODE_LENGTH
    ret = {}

    max_length = 0
    for k in candidates.iterkeys():
        if len(k) > max_length:
            max_length = len(k)

    MAX_BARCODE_LENGTH = max_length;
    modified_candidates = {}
    for k in candidates.iterkeys():
        BARCODE_PAD_LENGTH[candidates[k]] = max_length-len(k)
        to_append = 'N'*(max_length-len(k))
        modified_candidates[k+to_append] = candidates[k]

    for seq in modified_candidates.iterkeys():
        ret[seq] = modified_candidates[seq]
        recurseAllCombinations(seq,'',ret,modified_candidates[seq])
    return ret


def compute_average_Q_scores():
    flagQ = False
    average_q_R1 = sum(record_R1.letter_annotations["phred_quality"]) / float(len(record_R1))
    if average_q_R1 < 30:
        flagQ = True
    if paired:
        average_q_R2 = sum(record_R2.letter_annotations["phred_quality"]) / float(len(record_R2))
        if average_q_R2 < 30:
            flagQ *= True
    if not paired:
        average_q_R2 = 0
    return flagQ, average_q_R1, average_q_R2

def find_match (seq, table, possible_starts):
    for k in possible_starts:
        if seq[k:k+20] in table:
            return seq[k:k+20]
    return None

if __name__ == "__main__":
    
    ##### INITIALIZE OUTPUTS AND DATA STRUCTURE
    
    # Read STDIN (superseed parameter)
    if len(sys.argv)>1:
        file_ids = [int(id) for id in sys.argv[1].split("-")]
    
    # Determine if paired-end (ie is a R2 file provided?)
    paired = True
    if fastq_R2 == "":
        paired = False

    # get mapping
    mapping = {}  # Tuple of (id_R1, id_R2) maps to Description
    ids     = {"r1":[],"r2":[]}
    h = open(mapping_file)  # mapping.csv contains id_R1, id_R2, Description
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = int(data[0])
        data[1] = int(data[1])
        data[2] = str(data[2])
       # if data[2] in MAPPINGS_TO_EVALUATE:
        mapping[(data[0],data[1])] = data[2]
        if not data[0] in ids["r1"]:
            ids["r1"].append(data[0])
        if paired and not data[1] in ids["r2"]:
            ids["r2"].append(data[1])
    h.close()

    # get necessary barcodes according to mapping
    global bc_len
    barcodes = {"r1":{},"r2":{}}
    h = open(barcodes_file)  # barcodes.csv contains orientation,sequence,id
    h.readline()
    for l in h:
        data = l.strip().split(",")
        data[0] = str(data[0]).lower()
        data[1] = str(data[1]).upper()
        data[2] = int(data[2])
        if data[2] in ids[data[0]]:
            barcodes[data[0]][data[1]]=data[2]
            bc_len = len(data[1]) # all barcodes assumed to have same length!
                                  # Notice this assumption in above line!!!
    h.close()
     
    # write barcodes to fasta files (used for alignment)
    for direction in ("r1","r2"):
        h = open("%s/bc_%s.fas" % (output_folder, direction) ,"w")
        for seq in barcodes[direction]:
            h.write(">%s\n%s\n" % (barcodes[direction][seq],seq))
        h.close()
    
    
    # Produce degenerate barcodes to speed up matching via hash table
    degenerate_barcodes = {"r1": generateBarcodeHashTable(barcodes['r1']), "r2": {}}

    # I/O
    fastq_R1s   = {}
    h_dmpx_R1s  = {}
    if paired:
        fastq_R2s   = {}
        h_dmpx_R2s  = {}
    h_logs       = {}

    allOligos = {}
    db_fasta = open("OLS Master key - unique oligos 32992 with 5bp context.fas")  # Assumes structure is > line followed by non-> line
    for l in db_fasta:
        if ">" in l:
            name = l[1:].strip()
        if not ">" in l:
            allOligos[l.strip()[6:26] ] = name
            #print name, refs[name]
    db_fasta.close()

    print "The following file(s) will be analyzed:"

    id = 1
    # inputs
    fastq_R1s[id] = "%s_%03d.fastq" % (fastq_R1, id)
    print fastq_R1s[id]
    if paired:
        fastq_R2s[id] = "%s_%03d.fastq" % (fastq_R2, id)
        print fastq_R2s[id]

    # fastq outputs
    # This is all GC's original method of dumping to new fastq files
    h_dmpx_R1s[id] = {}
    if paired:
        h_dmpx_R2s[id] = {}


    # outFiles = {}  # Map from description to csv file handler
    # for m, v in mapping.iteritems():
    #     outFiles[v] = csv.writer(open("output/%s_reads.csv" % v,'wb'))

    #########   MAIN ANALYSIS
    for id in [1]:
        readperbc = {}
        for map in mapping:
            readperbc[mapping[map]] = 0

        total_count   = 0
        handle_matches_count = 0
        r1_promoter_matches_count = 0
        r2_promoter_matches_count = 0
        readperbc[mapping[map]] = 0
        records_R1 = SeqIO.parse(fastq_R1s[id], "fastq")
        if paired:
            records_R2 = SeqIO.parse(fastq_R2s[id], "fastq")

        csv_output = csv.writer(open('summary.csv','w'))
        csv_output.writerow(['Expt',
                             'R1_Seq','R1_aveQ','R1_prom_score','R1_sgRNA_match','R1_handle_score','R1_full_handle_score','R1_extended_score',
                             'R2_Seq','R2_aveQ','R2_prom_score','R2_sgRNA_match','R2_handle_score','R2_full_handle_score','R2_extended_score'
                             ])

        R1_fastq_out_if_no_sgRNA_but_handle_matches = open('no_sgRNA_but_has_handle_R1.fastq', 'w')
        if paired:
            R2_fastq_out_if_no_sgRNA_but_handle_matches = open('no_sgRNA_but_has_handle_R2.fastq','w')

        for record_R1 in records_R1:
            total_count += 1
            if paired:
                record_R2 = records_R2.next()
            if verbose and total_count % 1000 == 0:
                if paired:
                    printPairedStatus()
                else:
                    printUnpairedStatus()

            # Assess quality of the reads
            # Only keep paired reads with one mate average score >= Q30
            flagQ, average_Q_R1, average_Q_R2 = compute_average_Q_scores()
            if flagQ:
                continue

            bc_R1, pos_R1, score_R1 = find_degenerate_barcode(record_R1.seq[start_bc:stop_bc], "r1", bc_offset)
            # if still fail look for pairwise alignment (extend the region searched to allow for 2 insertion)
            if bc_R1 == -1:
                bc_R1, pos_R1, score_R1 = water_bc(record_R1.seq[start_bc:stop_bc+2], "r1", id, aln_threshold)

            pos_R2 = 0  # since there is no barcode for read_2
            # Check that we can map this read's multiplex barcode
            map_key = (bc_R1, 0)
            if bc_R1 == -1:
                #raise  Exception('should already be de-multiplexed')
                continue
            elif not map_key in mapping:
                continue

            read1_after_promoter = str(record_R1[pos_R1 + ROI_start_R1:].seq)
            read2_after_promoter = str(record_R2[pos_R2 + ROI_start_R2:].seq)

            r1_match = find_match(str(record_R1.seq).lower(), allOligos, xrange(pos_R1+ROI_start_R1-2,pos_R1+ROI_start_R1+3))
            if not r1_match:
                SeqIO.write(record_R1[pos_R1+ROI_start_R1-ROI_offset: pos_R1+ROI_start_R1+ROI_length+ROI_offset], open('temp.fastq','wb'),'fastq')
                r1_match = compute_bowtie_results('temp.fastq')

            if paired:
                r2_match = find_match(str(record_R2.seq).lower(), allOligos, xrange(pos_R2+ROI_start_R2-2,pos_R2+ROI_start_R2+3))
                r2_bowtie = ''
                if not r2_match:
                    SeqIO.write(record_R2[pos_R2+ROI_start_R2-ROI_offset: pos_R2+ROI_start_R2+ROI_length+ROI_offset], open('temp.fastq','wb'),'fastq')
                    r2_match = compute_bowtie_results('temp.fastq')

            r1_handle = str(record_R1[pos_R1+ROI_start_R1+20:].seq).upper()
            r1_promoter = str(record_R1[pos_R1 : pos_R1+len(R1_EXPECTED_PROMOTER) ].seq)
            r1_promoter_match = pairwise2.align.localxx(R1_EXPECTED_PROMOTER, r1_promoter)[0]

            if r1_promoter_match[2] > 20:
                r1_promoter_matches = True
                r1_promoter_matches_count += 1
            else:
                r1_promoter_matches = False

            if paired:
                r2_promoter = str(record_R2[6: 6+len(R2_EXPECTED_PROMOTER) ].seq)
                r2_handle = str(record_R2[pos_R2+ROI_start_R2+20:].seq).upper()
                r2_promoter_match = pairwise2.align.localxx(R2_EXPECTED_PROMOTER, r2_promoter)[0]
                if r2_promoter_match[2] > 20:
                    r2_promoter_matches = True
                    r2_promoter_matches_count += 1
                else:
                    r2_promoter_matches = False

                r2_cas9_match = pairwise2.align.globalxx(CAS9_HANDLE[0:len(r2_handle)], r2_handle)[0]
                if r2_cas9_match[2] >= 16:
                    r2_handle_matches = True
                else:
                    r2_handle_matches = False

            r1_cas9_match = pairwise2.align.localxx(CAS9_HANDLE[0:len(r1_handle)+1], r1_handle)[0]
            r1_handle_matches = False
            r1_full_handle_score = pairwise2.align.localxs('ATTTTAACTTGCTATTTCTAGCTCTAAAAC', read1_after_promoter,-2,-.1)[0][2]
            r1_score_after_prom = pairwise2.align.localxs(POST_sgRNA_SEQ_1, read1_after_promoter,-2, -.1 )[0][2]

            if paired:
                r2_cas9_match = pairwise2.align.localxx(CAS9_HANDLE[0:len(r2_handle)+1], r2_handle)[0]
                r2_handle_matches = False
                r2_full_handle_score = pairwise2.align.localxs('ATTTTAACTTGCTATTTCTAGCTCTAAAAC', read2_after_promoter,-2,-.1)[0][2]
                r2_score_after_prom = pairwise2.align.localxs(POST_sgRNA_SEQ_2, read2_after_promoter,-2, -.1 )[0][2]


            if not r1_match and r1_handle_matches:
                SeqIO.write(record_R1[pos_R1+ROI_start_R1-ROI_offset:pos_R1+ROI_start_R1+ROI_length+ROI_offset],R1_fastq_out_if_no_sgRNA_but_handle_matches,'fastq')
            if paired and not r2_match and r2_handle_matches:
                SeqIO.write(record_R1[pos_R2+ROI_start_R2-ROI_offset:pos_R2+ROI_start_R2+ROI_length+ROI_offset],R2_fastq_out_if_no_sgRNA_but_handle_matches,'fastq')
            csv_output.writerow([ mapping[map_key],
                                 read1_after_promoter, average_Q_R1, r1_promoter_match[2], r1_match, r1_cas9_match[2], r1_full_handle_score, r1_score_after_prom,
                                 read2_after_promoter, average_Q_R2, r2_promoter_match[2], r2_match, r2_cas9_match[2], r2_full_handle_score, r2_score_after_prom,
                                 ])
            #print r1_match, r2_match
              #outFiles[mapping[map_key]].writerow([record_R1[pos_R1+ROI_start_R1:pos_R1+ROI_start_R1+20], record_R2[30:50]])
            readperbc[mapping[map_key]] += 1

        ### Final print to screen of status
        if paired:
            printPairedStatus()
        else:
            printUnpairedStatus()


