#!/bin/python3

from subprocess import getoutput

barcodes = []

for num in range(1,13):
    if num < 10:
        barcodes.append('barcode0' + str(num))
    else:
        barcodes.append('barcode' + str(num))

barcodes.append('barcode12a')
barcodes.append('unclassified')

filtering = ['fail', 'pass']

total = 0

for filt in filtering:
    for barcode in barcodes:
        nbline = getoutput('zcat ' + filt + '/' + barcode + '/*.fastq.gz | wc -l').split(' ')[0]
        nbreads = int(nbline) / 4
        total = total + nbreads
        print(filt + '/' + barcode + ' = ' + str(nbreads))
print('Total = ' + str(total))
