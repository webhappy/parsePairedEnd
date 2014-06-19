db_fasta = open("../OLS Master key - unique oligos 32992 with 5bp context.fas")  # Assumes structure is > line followed by non-> line
out_file = open('all_oligos.fas','wb')
for l in db_fasta:
    if ">" in l:
        name = l[1:].strip()
    if not ">" in l:
        sgRNA_seq = l.strip()[6:26]
        out_file.write('>'+sgRNA_seq+'\n'+sgRNA_seq+'\n')
        #print name, refs[name]
db_fasta.close()