# input python to terminal and insert the following code

# read the single record
from Bio import SeqIO

seq = SeqIO.read("one.fasta", "fasta")
print(type(seq))
print(seq)

for s in seq:
    print(type(s))
    print(s)

# result
'''
ID: one
Name: one
Description: one
Number of features: 0
Seq('ATCGTACGATCGATCGATCGCTAGACGTATCG', SingleLetterAlphabet())
'''

# read the multiple records
from Bio import SeqIO

seqs = SeqIO.parse("two.fasta", "fasta") # Using .read() in multiple records cause an error
print(type(seqs))
print(seqs)

for s in seqs:
    print(type(s))
    print(s)

# result
''''
<class 'Bio.SeqRecord.SeqRecord'>
ID: one
Name: one
Description: one
Number of features: 0
Seq('AAAACCCCGGGGTTTTACGTACGTACGTACGT')
<class 'Bio.SeqRecord.SeqRecord'>
ID: two
Name: two
Description: two
Number of features: 0
Seq('AAAACCCCGGGGTTTTACGTACGTACGTACGT')'
'''

# 출처: https://korbillgates.tistory.com/202 [생물정보학자의 블로그:티스토리] 
# [바이오파이썬] 5.1.1 SeqIO 모듈로 서열 파일 읽기 - FASTA
# 모두 et cetra 폴더로 옮김.