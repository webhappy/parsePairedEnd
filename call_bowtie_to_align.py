import subprocess

#bowtie settings unused
# '-N','1','--score-min','G,25,8',

def compute_bowtie_results (file_name='temp.fastq'):
    p = subprocess.Popen(['bowtie2-2.2.2/bowtie2','--local', '-x', 'bowtie2-2.2.2/all_oligos', file_name],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    last_line = ''
    for line in p.stdout:
        if line.strip() != '':
            last_line = line.strip()

    results = last_line.split('\t')

    #print last_line
    if results[2] == '*':
        return None
    return results[2]