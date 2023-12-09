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

Given a pattern P[1...m], let p denote its corresponding decimal value. In a similar manner, given a text T[1...n], let t_{s} denote the decimal value of the length-m substring T[s+1...s+m], for s=0,1,...,n-m. Certainly, t_{s} = p if and only if T[s+1...s+m]=P[1..m]; thus, s is a valid shift if and only if  t_{s}=p. If we could compute p in time O(m) and all the t_{s} values in a total of O(n-m+1) time, then we could determine all valid shift s in time O(m) + O(n-m+1) = O(n) by compairing p with each of the t_{s} values. We can compute p in time O(m) using Horner's rule:
```
p = P[m] + 10(P[m-1] + 10(P[m-2]+...+10(P[2] + 10P[1])...))
```
Similarly, we can compute t_{0} from T[1...m] in time O(m).
To compute the remaining values t1, t2, ..., t(n-m), we observe that we can compute t_{s+1} from t_{s} in constant time, since
t_{s+1} = 10(t_{s} - 10^(m-1)T[s+1]) + T[s+m+1] (recurence equation)

To reduce time we also need to reasonobly asuume that each arithmetic operation on p (which is m digits long) takes "constant time". To solve this: compute p and the t_{s} values modulo a suitable modulus q. We can compute p modulo q in O(m) time and all the t_{s} values modulo q in O(n - m +1) time. If we choose the modulus q as a prime such that 10q just fits within one computer word, then we can perform all the the necessary computations with single-precision arithmetic. In general, with a d-ary alphabet {0,1,...,d-1}, we choose q so that dq fits within a computer word and adjust the recurence equation to work modulo q, so that it becomes
t_{s+1} = (d(t_{s} - T[s+1]h) + T[s+m+1])mod(q)
where h = d^(m-1) (mod(q)) is the value of the digit "1" in the high-order position of an m-digit text window. 
If q is large enough, then we hope that spurious hits occur infrequently enough that the cost of the extra checking is low. 

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
