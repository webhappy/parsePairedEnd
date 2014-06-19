import numpy as np
import mysql.connector
import pandas as pd
import time



def process_both_ends (file_for, file_rev):
    cnx = mysql.connector.connect(user='davidc', password='mysql_password', host='127.0.0.1',  database='CRISPR')
    sgRNAs = pd.io.sql.read_sql('select id, gene_name, seq FROM sgRNAs',cnx, index_col='seq')
    sgRNAs['id'] = sgRNAs['id'] - 1  # so id's go from 0 to (#-1)
    num_sgRNAs = sgRNAs.shape[0]
    genes = sorted(sgRNAs.groupby('gene_name').groups.keys())
    NUM_GENES = len(genes)
    id_to_sgRNA_seq = pd.DataFrame(sgRNAs.index, sgRNAs['id'])
    #sgRNA_to_gene = sgRNAs['gene_name'].to_dict()  # using pandas takes almost 3 seconds longer than converting to dict!!!

    colnames = ['read_name','flags','sgRNA','offset','qual','CIGAR',
                        'na','na2','na3','read_seq','read_qual',
                        'as','xs','xm','x0','xg','nm','md','yt']
    f_for = open(file_for, 'r')
    f_rev = open(file_rev, 'r')

    total_count = 0
    for_count = 0
    rev_count = 0
    both_count = 0
    both_map_to_genes = 0

    #counts = np.zeros(( 32992 * 32992, 2))
    genes_to_idx = {}
    for i in xrange(len(genes)):
        genes_to_idx[genes[i]] = i
    counts_genes = np.zeros((NUM_GENES,NUM_GENES))

    for line_for in f_for:
        total_count += 1
        line_rev = f_rev.next()  # iterate reverse end so it corresponds to forward end
        cols_for = line_for.split('\t')
        cols_rev = line_rev.split('\t')
        if cols_for[1] == '0':
            for_mapped = True
            for_gene = sgRNAs.at[cols_for[2], 'gene_name']
            #for_gene = sgRNA_to_gene[cols_for[2]]
        else:
            for_mapped = False

        if cols_rev[1] == '0':
            rev_mapped = True
            rev_gene = sgRNAs.at[cols_rev[2], 'gene_name']
            #rev_gene = sgRNA_to_gene[cols_rev[2]]
        else:
            rev_mapped = False

        if for_mapped:
            for_count += 1
        if rev_mapped:
            rev_count += 1
        if for_mapped and rev_mapped:
       #     if for_gene and rev_gene:
      #          both_map_to_genes += 1
     #           counts_genes[genes_to_idx[for_gene], genes_to_idx[rev_gene]] += 1  # add the count at the matrix position corresponding to these genes
            both_count += 1

    print '''
    Total count: %i
    For_mapped: %i
    Rev_mapped: %i
    Both map to genes: %i
    Both_mapped: %i''' % (total_count, for_count, rev_count, both_map_to_genes, both_count)

    # sgRNA_results_file = open('sgRNA_counts.csv','w')
    # sgRNA_results_file.write('sgRNA1, sgRNA2, counts_orientation_1, counts_orientation_2, counts_sum')
    # cur_idx = 0
    # for j in xrange(num_sgRNAs):
    #     for k in xrange(num_sgRNAs):
    #         assert get_sgRNA_index(id_to_sgRNA_seq[j], id_to_sgRNA_seq[k])[1]
    #         sgRNA_results_file.write(counts[cur_idx][0] + ',' + counts[cur_idx][1]+','+(counts[cur_idx][0]+counts[cur_idx][1]) )
    #         cur_idx += 1

#start = time.time()
#elapsed = (time.time() - start)
#print "Took",elapsed,'time'