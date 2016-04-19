#! usr/bin/env python

from collections import Counter, defaultdict
import pybedtools

lamina_bed = '../data-sets/bed/lamina.bed'
lamina = pybedtools.BedTool(lamina_bed)

#On what chrom is the region with the largest start value?

largest_start = 0

for record in lamina:
    if largest_start < record.start:
        largest_start = record.start
        start_chrom = record.chrom
print 'answer-1:', start_chrom

#What is the region with the largest end value on chrY?

largest_stop = 0

for record in lamina:
    if record.chrom == 'chrY':
        if largest_stop < record.end:
            largest_stop = record.end
            stop_chrom = record.chrom
            stop_start = record.start
            stop_value = record.fields[3]

print 'answer-2:', stop_chrom, stop_start, largest_stop, \
    stop_value, largest_stop - stop_start

##PROBLEM 2 (FASTQ FILES)

filename = "../data-sets/fastq/SP1.fq"

#parse_fastq will take a fastq file on imput and yield the name, seq, and
#quality for each record

def parse_fastq(filename, num_of_records):
    line_num = 0
    for line in open(filename):
        if line_num > (num_of_records *4): break
        line_type = line_num % 4
        if line_type == 0:
            name = line.strip()
        elif line_type == 1: 
            seq = line.strip()
        elif line_type == 3:
            quals = line.strip()
            
            yield name, seq, quals
      
        line_num +=1

#sum_quals take quality scores on imput and will output numerical scores
#for each character in the quality score

def sum_quals(qual): 
    val = sum([ord(i) for i in qual])
    return val


#reverse_complement will take a seq on imput and otput the reverse
#complement sequence

def reverse_complement(seq):
    comp = []
    for char in seq: 
        if char == 'A':
           comp.append('T')
        elif char == 'G':
            comp.append('C')
        elif char == 'T':
            comp.append('A')
        elif char == 'C': 
            comp.append('G')
    return ''.join(reversed(comp))

#Which of the first 10 seq records has the largest number fo cs

CCount = 0
for name, seq, quals in parse_fastq(filename, 10):
    if CCount < Counter(seq)['C']:
        CCount = Counter(seq)['C']
        Cname = name
        Cseq = seq
print 'answer-3:', Cname

#For each record, convert each character quality into a number and sum the
#numbers. Report the largest quality score

Qscore = 0
for name, seq, quals in parse_fastq(filename, 'inf'):
    if Qscore < sum_quals(quals):
        Qscore = sum_quals(quals)

print 'answer-4:', Qscore

#Report the reverse complement for each of the first 10 records 

rev_comp = []
for name, seq, quals in parse_fastq(filename, 10):
    rev_comp.append(reverse_complement(seq))

print 'answer-5:', rev_comp
