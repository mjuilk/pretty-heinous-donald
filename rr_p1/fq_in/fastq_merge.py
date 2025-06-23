#!/usr/bin/python

import gzip, os, sys
path = sys.argv[1]
print(path)
outfile = open(path.split('/')[-2] + '.fastq', 'w')
fns = [fn for fn in os.listdir(path + '/') if fn.endswith('fastq.gz')]
print(len(fns))
counter = 0
nr = len(fns)

for fn in fns:
    fullpath = path + '/' + fn
    file = gzip.open(fullpath)
    print(fullpath)
    print(round((counter / nr) * 100, 2))
    counter += 1
    for record in file:
        outfile.write(record.decode("utf-8"))
