from parameters import *
from Bio import SeqIO
import csv
from subprocess import call

mapping = {}  # Tuple of (id_R1, id_R2) maps to Description
ids     = {"r1":[],"r2":[]}
h = open(mapping_file)  # mapping.csv contains id_R1, id_R2, Description
h.readline()
for l in h:
    data = l.strip().split(",")
    data[0] = int(data[0])
    data[1] = int(data[1])
    data[2] = str(data[2])
    mapping[(data[0],data[1])] = data[2]
    if not data[0] in ids["r1"]:
        ids["r1"].append(data[0])
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

fastq_R1s   = {}
h_dmpx_R1s  = {}
fastq_R2s   = {}
h_dmpx_R2s  = {}
h_logs       = {}

allOligos = {}  # 20-mer to oligo name
# files = ["./createFastaFiles/List of FBA sgRNAs.csv",'./createFastaFiles/List of sRNA sgRNAs.csv']
# for f in files:
#     inFile = csv.reader(open(f,'rU'))
#     for (oligoName, seq) in inFile:
#         allOligos[seq] = oligoName

db_fasta = open("OLS Master key - unique oligos 32992 with 5bp context.fas")  # Assumes structure is > line followed by non-> line
for l in db_fasta:
    if ">" in l:
        name = l[1:].strip()
    if not ">" in l:
        allOligos[l.strip()[6:26] ] = name
        #print name, refs[name]
db_fasta.close()

#  --- Done reading the list of possible oligos to check against

for id in file_ids:
    # inputs
    fastq_R1s[id] = "%s_%03d.fastq" % (fastq_R1, id)
    fastq_R2s[id] = "%s_%03d.fastq" % (fastq_R2, id)

    # fastq outputs
    # This is all GC's original method of dumping to new fastq files
    h_dmpx_R1s[id] = {}
    h_dmpx_R2s[id] = {}

    for map in mapping:
        h_dmpx_R1s[id][map] = open("%s/R1_%s_%03d_dmpx.fastq" % (output_folder, mapping[map], id), 'r')
        h_dmpx_R2s[id][map] = open("%s/R2_%s_%03d_dmpx.fastq" % (output_folder, mapping[map], id), 'r')

# Now read in the files that were previously made
for id in file_ids:
    for map in mapping:
        outFile = csv.writer(open("results_"+mapping[map]+'.csv','wb'))
        outFile.writerow(['R1_match','R1_score','R2_match','R2_score'])

        num_R1_maps = 0
        num_R2_maps = 0
        num_both_map = 0
        total_count = 0

        r1 = h_dmpx_R1s[id][map]
        r2 = h_dmpx_R2s[id][map]

        records_R1 = SeqIO.parse(r1, "fastq")
        records_R2 = SeqIO.parse(r2, "fastq")

        for record_R1 in records_R1:
            total_count += 1
            record_R2 = records_R2.next()

            #seq for R1
            r1_match = find_match(str(record_R1.seq).lower(), allOligos, xrange(2,9))
            r2_match = find_match(str(record_R2.seq).lower(), allOligos, xrange(2,9))
            average_q_R1 = sum(record_R1.letter_annotations["phred_quality"])/float(len(record_R1))
            average_q_R2 = sum(record_R2.letter_annotations["phred_quality"])/float(len(record_R2))
            outFile.writerow([r1_match,average_q_R1,r2_match,average_q_R2])
        #    r1_match = find_match(str(record_R1.seq).lower(), allOligos, xrange(5,6))
        #    r2_match = find_match(str(record_R2.seq).lower(), allOligos, xrange(5,6))

            if r1_match:
                num_R1_maps += 1
            else:
                print average_q_R1
            if r2_match:
                num_R2_maps += 1
            else:
                print average_q_R1
            if r1_match and r2_match:
                num_both_map += 1

            if total_count % 100000 == 0:
                print "Current count:",total_count,', # that R1 maps:',num_R1_maps,', # that R2 maps:',num_R2_maps,', and # that both map:',num_both_map

        print "Total count:",total_count,', # that R1 maps:',num_R1_maps,', # that R2 maps:',num_R2_maps,', and # that both map:',num_both_map
