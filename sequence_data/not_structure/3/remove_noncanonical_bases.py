from Bio import SeqIO, SeqRecord
import sys, os

def remove_sequences_with_nonbases(infile):
    bases = ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g']
    in_handle = open(infile, "r")
    records = list (SeqIO.parse(in_handle, "fasta"))
    new_records = []
    for record in records:
        keep = True
        for base in record.seq:
            if base not in bases:
                keep = False
                break
        if keep:
            new_records.append(record)

    SeqIO.write(new_records, 'cleaned_nucleotide.fasta', 'fasta')

    return 0

def main():
    infile = sys.argv[1]
    remove_sequences_with_nonbases(infile)

main()
