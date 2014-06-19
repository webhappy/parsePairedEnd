__author__ = 'davidc'
import csv

#inFile = csv.reader(open('List of FBA sgRNAs.csv','rU'))
#outFile = open('FBA.fasta','wu')

files = ['List of FBA sgRNAs.csv','List of sRNA sgRNAs.csv']
outFile = open('../combined.fasta','wu')

for f in files:
    inFile = csv.reader(open(f,'rU'))

    for (oligoName, seq) in inFile:
        outFile.write('>'+oligoName+'\n')
        outFile.write(seq+'\n')