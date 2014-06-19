import subprocess
from parameters import *
import os
import parse_both_ends

# subprocess.call('python demultiplex_fastq.py', shell=True)
print "Done demultiplexing, now aligning with bowtie2"
mapping = {}  # Tuple of (id_R1, id_R2) maps to Description
h = open(mapping_file)  # mapping.csv contains id_R1, id_R2, Description
h.readline()
for l in h:
    data = l.strip().split(",")
    data[0] = int(data[0])
    data[1] = int(data[1])
    data[2] = str(data[2])
    mapping[(data[0],data[1])] = data[2]
h.close()

id = 1
for m in mapping:
    input_1 ="%s/R1_%s_%03d_dmpx.fastq" % (output_folder, mapping[m], id)
    input_2 = "%s/R2_%s_%03d_dmpx.fastq" % (output_folder, mapping[m], id)
    out_1 = "%s/R1_%s_%03d.sam" % (output_folder, mapping[m], id)
    out_2 = "%s/R2_%s_%03d.sam" % (output_folder, mapping[m], id)
    print "Calling bowtie from", input_1,"and",input_2
    subprocess.call('bowtie2-2.2.2/bowtie2 --norc --local -p 4 --no-hd --reorder -x bowtie2-2.2.2/all_oligos %s -S %s'%
                    (input_1, out_1), shell=True)
    subprocess.call('bowtie2-2.2.2/bowtie2 --norc --local -p 4 --no-hd --reorder -x bowtie2-2.2.2/all_oligos %s -S %s'%
                    (input_2, out_2), shell=True)
    parse_both_ends.process_both_ends(out_1,out_2)