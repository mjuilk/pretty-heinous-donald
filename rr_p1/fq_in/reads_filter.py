#########################
# script for filtering  #
# from a fastq file     #
# into two haplotypes   #
#########################

import gzip, os, sys
import pandas as pd

#path for fastq
fastqin = sys.argv[1]
print(f"Input fastq file : {fastqin}")

#list of reads
filename = sys.argv[2]
df = pd.read_csv(filename)
haplols1 = [x for x in list(df.iloc[:,0]) if str(x) != 'nan']
haplols2 = [x for x in list(df.iloc[:,1]) if str(x) != 'nan']

#filtering reads    
with gzip.open(fastqin, 'r') as all_reads:
    print("decoding reads...")
    all_reads_lines = [ln.decode("utf-8") for ln in all_reads]
    outfile1 = open('haplo1.fastq', 'w')
    outfile2 = open('haplo2.fastq', 'w')
    counter = 0
    print("looping through reads now...")
    for ln in all_reads_lines:
        if (ln.startswith('@')) and (ln.rstrip('\n')[1:] in haplols1):
            print(ln)
            idx = all_reads_lines.index(ln)
            joined = ''.join(all_reads_lines[idx : idx + 4])
            outfile1.write(joined)
        if (ln.startswith('@')) and (ln.rstrip('\n')[1:] in haplols2):
            print(ln)
            idx = all_reads_lines.index(ln)
            joined = ''.join(all_reads_lines[idx : idx + 4])
            outfile2.write(joined)
        else:
            pass

        if str(counter)[-1] == "0":
            print(f"{round(counter / len(all_reads_lines) * 100, 2)} of input FASTQ file read")
        counter += 1
    outfile1.close()
    outfile2.close()
