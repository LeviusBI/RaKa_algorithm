# Rabin-Karp algorithm in python

import sys
from Bio import SeqIO


def Rabin_Karp(pattern, text, q, d):
    m = len(pattern)
    n = len(text)
    if m > n:
        raise PatternLength("Pattern P is longer than text T. Insert shorter pattern")
    t = 0
    h = 1
    i = 0
    j = 0

    for i in range(m-1):
        h = (h*d) % q

    # Calculate hash value for pattern and text
    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q

    # Find the match
    for i in range(n-m+1):
        if p == t:
            for j in range(m):
                if text[i+j] != pattern[j]:
                    break

            j += 1
            if j == m:
                print("Pattern is found at position: " + str(i+1))

        if i < n-m:
            t = (d*(t-ord(text[i])*h) + ord(text[i+m])) % q

            if t < 0:
                t = t+q


text, pattern = sys.argv[1], sys.argv[2]

fasta_text = SeqIO.parse(open(text), 'fasta')
fasta_pattern = SeqIO.parse(open(text), 'fasta')
pattern = str(fasta_pattern.seq)
sequence = str(fasta_text.seq)



pattern = "CDD"
Rabin_Karp(pattern, sequence, q=13, d=10)
