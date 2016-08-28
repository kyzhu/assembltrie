#Author: kzhu
#Test whether decompress is correct

#import random

def readFasta(fp):
    seq = []
    for line in fp:
        line = line.rstrip()
        seq.append(''.join(line))
    return seq

def isEqual(seq1, seq2):
    seq1.sort()
    seq2.sort()
    for i in range (0, len(seq1)):
	if seq1[i] != seq2[i]: return False
    return True

def isEqualN(seq1, seq2):
    for i in range (0, len(seq1)):
	seq1[i].replace('N', 'A')
    seq1.sort()
    seq2.sort()
    for i in range (0, len(seq1)):
	#if 'N' in seq1[i]: continue
	if seq1[i] != seq2[i]: return False
    return True

fp = open("chr2.fa", 'r')
seq1 = readFasta(fp)
fp.close()
fp = open("decompressed.out", 'r')
seq2 = readFasta(fp)
fp.close()
print(isEqual(seq1, seq2))
