# RaKa_algorithm
## This repository contains code related to Rabin-Karp's algorithm of aligning of sequence P to text T.
### Group members:
1. Igamberdiev Lev
2. Tarasov Vadim
3. Kaplun Antoniy
4. Zuev Andrey
###### All group members are studying at [MIPT](https://mipt.ru/), 3rd year of BcS.
#### Task description
```
Implement Rabin-Karp's algorithm to substring P alignment to string T. Input format - FASTA. Ouput format - locus of P in T.
Example:

Input: P.fa, T.fa
Output: Alignment locus

Sample_input:
>P
ATC
>T
AATCCG
Sample_output:
2
```
### Algorithm description
Rabin and Karp proposed a string-matchind algorithm that performs well in practice. The Rabin-Karp algorithm uses O(m) preprocessing time (m = P.length), and its worst-case running time is O((n - m +1)m) (n = T.length). 

Given a pattern P[1...m], let p denote its corresponding decimal value. In a similar manner, given a text T[1...n], let t_{s} denote the decimal value of the length-m substring T[s+1...s+m], for s=0,1,...,n-m. Certainly, t_{s} = p if and only if T[s+1...s+m]=P[1..m]; thus, s is a valid shift if and only if  t_{s}=p

### Pseudocode
```
n = T.length
m = P.length
h = d^{m-1}mod(q)
p = 0
t_{0} = 0
for i = 1 to m #preprocessing
  p = (dp + P[i])mod(q)
  t_{0} = (dt_{0} + T[i])mod(q)
for s = 0 to n-m
  if p==t_{s}
    if P[1...m]==T[s+1...s+m]
      print "Pattern occurs with shift" s
  if s < n-m
    t_{s+1} = (d(t_{s} - T[s+1]h) + T[s+m+q])mod(q)
```
