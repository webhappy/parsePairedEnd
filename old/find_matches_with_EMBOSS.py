from subprocess import call
from Bio import SeqIO, AlignIO
from Bio.Alphabet import generic_dna

def find_matches_supermatcher (seq, reference_file='combined.fasta'):
    tmp = open('tmp.txt','w')
    tmp.write(seq+'\n')
    tmp.close()
    call('supermatcher tmp.txt %s -gapopen 10 -gapextend 10 -aformat3 simple -outfile out_detailed.supermatcher -minscore 40' % reference_file, shell=True)
    call('supermatcher tmp.txt %s -gapopen 10 -gapextend 10 -aformat3 score -outfile out.supermatcher -minscore 40' % reference_file, shell=True)

test_seq = 'TNCATTGAGACCCATGTGAAGGTCAGGATT'.lower()

#find_matches_supermatcher(test_seq)

t = AlignIO.parse('out.supermatcher','emboss')
print t
for j in t:
    print j