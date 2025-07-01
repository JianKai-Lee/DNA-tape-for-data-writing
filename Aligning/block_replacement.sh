#!/bin/bash
echo "Program starts"
echo "Processing reads with cutadapt in four steps: F, R, RC(F), RC(R)"

FASTQ="merged.fastq"
CUTADAPT_F_FASTQ="merged_cut_F.fastq"  
CUTADAPT_FR_FASTQ="merged_cut_FR.fastq"  
CUTADAPT_RC_F_FASTQ="merged_cut_RC_F.fastq"  
CUTADAPT_RC_FR_FASTQ="merged_cut_RC_FR.fastq"  
CUTADAPT_FINAL_FASTQ="merged_cutprimers.fastq"
ANNOTATED_FASTQ="merged_cut_annotated.fastq"  
BACK_REPLACED_FASTQ="Back_Replaced.fastq"  
PRIMER_F="TTCTCTTATT"
PRIMER_R="ATATCCCATA"
POOL="pool.fasta"
OUTPUT_DIR="barcode_fastq_files"

mkdir -p $OUTPUT_DIR

PRIMER_F_RC=$(echo $PRIMER_F | rev | tr ATGC TACG)
PRIMER_R_RC=$(echo $PRIMER_R | rev | tr ATGC TACG)

echo "Step 1: Trimming 5' primer (F)..."
cutadapt -g "$PRIMER_F" -e 0 -O 10 -m 0 --untrimmed-output untrimmed_reads_F.fastq -o $CUTADAPT_F_FASTQ $FASTQ > cutadapt_F_log.txt

echo "Step 2: Trimming 3' primer (R) from F-trimmed file..."
cutadapt -a "$PRIMER_R" -e 0 -O 10 -m 0 --discard-untrimmed -o $CUTADAPT_FR_FASTQ $CUTADAPT_F_FASTQ > cutadapt_R_log.txt

echo "Step 3: Trimming 3' reverse-complement primer (RC(F))..."
cutadapt -a "$PRIMER_F_RC" -e 0 -O 10 -m 0 --untrimmed-output untrimmed_reads_RC_F_.fastq -o $CUTADAPT_RC_F_FASTQ $FASTQ > cutadapt_RC_F_log.txt

echo "Step 4: Trimming 5' reverse-complement primer (RC(R)) from RC(F)-trimmed file..."
cutadapt -g "$PRIMER_R_RC" -e 0 -O 10 -m 0 --discard-untrimmed -o $CUTADAPT_RC_FR_FASTQ $CUTADAPT_RC_F_FASTQ > cutadapt_RC_R_log.txt

echo "Step 5: Converting reverse-complement sequences to forward sequences..."
python3 - <<EOF
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

with open("$CUTADAPT_RC_FR_FASTQ", "r") as infile, open("merged_cut_RC_FR_forward.fastq", "w") as outfile:
    lines = infile.readlines()
    for i in range(0, len(lines), 4):
        header = lines[i]
        seq = lines[i+1].strip()
        qual = lines[i+3].strip()
        
        forward_seq = reverse_complement(seq)
        outfile.write(header)
        outfile.write(forward_seq + "\n")
        outfile.write("+\n")
        outfile.write(qual + "\n")
EOF

echo "Step 6: Merging forward and reverse-complement processed files..."
cat $CUTADAPT_FR_FASTQ merged_cut_RC_FR_forward.fastq > $CUTADAPT_FINAL_FASTQ

echo "Step 7: Replacing repetitive O and E sequences with strict matching rules and outputting as FASTQ..."

python3 - <<EOF
import Levenshtein

def is_valid_match(seq, pattern):
    distance = Levenshtein.distance(seq, pattern)
    return distance == 0 or distance == 1

def replace_OE(seq, qual):
    O_seq = "ATCATTTTTC"
    E_seq = "ATTACACTAC"
    modified_seq = list(seq) 
    modified_qual = list(qual) 
    seq_len = len(seq)

    i = 0
    has_OE = False
    while i <= seq_len - len(O_seq):
        sub_seq = seq[i:i + len(O_seq)]

        if is_valid_match(sub_seq, O_seq):
            modified_seq[i] = "O"
            modified_qual[i] = "I"
            for j in range(1, len(O_seq)):
                modified_seq[i + j] = ""
                modified_qual[i + j] = ""
            i += len(O_seq)
            has_OE = True
            continue

        if is_valid_match(sub_seq, E_seq):
            modified_seq[i] = "E"
            modified_qual[i] = "I"  
            for j in range(1, len(E_seq)):  
                modified_seq[i + j] = ""
                modified_qual[i + j] = ""
            i += len(E_seq)
            has_OE = True
            continue

        i += 1

    modified_seq = "".join([c for c in modified_seq if c])
    modified_qual = "".join([c for c in modified_qual if c])

    return has_OE, modified_seq, modified_qual

with open("$CUTADAPT_FINAL_FASTQ", "r") as infile, open("$ANNOTATED_FASTQ", "w") as outfile:
    lines = infile.readlines()
    for i in range(0, len(lines), 4):
        header = lines[i]
        seq = lines[i+1].strip()
        qual = lines[i+3].strip()
        has_OE, modified_seq, modified_qual = replace_OE(seq, qual)

        if has_OE:
            outfile.write(header)
            outfile.write(modified_seq + "\n")
            outfile.write("+\n")
            outfile.write(modified_qual + "\n")
EOF

echo "Step 8: Replacing O/E back to original sequences with highest quality scores..."

python3 - <<EOF
def is_valid_OE_pattern(seq):
    if not seq:
        print("Empty sequence")
        return False

    if seq[0] not in {"O", "E"}:
        print(f"Invalid start character: {seq[0]} in sequence: {seq}")
        return False

    if seq[-1] not in {"O", "E"}:
        print(f"Invalid end character: {seq[-1]} in sequence: {seq}")
        return False

    prev_char = seq[0]
    prev_index = 0

    for i, char in enumerate(seq):
        if char in {"O", "E"}:
            if i == 0:
                continue

            if char == prev_char:
                print(f"Non-alternating characters: {char} at position {i} in sequence: {seq}")
                return False

            distance = i - prev_index
            if distance < 5 or distance > 7:
                print(f"Invalid distance: {distance} between {prev_char} and {char} in sequence: {seq}")
                return False

            prev_char = char
            prev_index = i

    print(f"Valid sequence: {seq}")
    return True

def replace_back(seq, qual):
    O_seq = "ATCATTTTTC"
    E_seq = "ATTACACTAC"
    modified_seq = []
    modified_qual = []
    i = 0
    while i < len(seq):
        if seq[i] == "O":
            modified_seq.append(O_seq)
            modified_qual.append("I" * len(O_seq))
            i += 1
        elif seq[i] == "E":
            modified_seq.append(E_seq)
            modified_qual.append("I" * len(E_seq))
            i += 1
        else:
            modified_seq.append(seq[i])
            modified_qual.append(qual[i])
            i += 1
    return "".join(modified_seq), "".join(modified_qual)

with open("$ANNOTATED_FASTQ", "r") as infile, open("$BACK_REPLACED_FASTQ", "w") as outfile:
    lines = infile.readlines()
    for i in range(0, len(lines), 4):
        header = lines[i]
        seq = lines[i+1].strip()
        qual = lines[i+3].strip()

        if is_valid_OE_pattern(seq):
            modified_seq, modified_qual = replace_back(seq, qual)
            outfile.write(header)
            outfile.write(modified_seq + "\n")
            outfile.write("+\n")
            outfile.write(modified_qual + "\n")
EOF

echo "Step 9: Adding flanking sequences to the beginning and end of each read..."

python3 - <<EOF
def add_flanking(seq, qual):

    flanking_F = "GGACGGATTGACAGTCGGATTTCTCTTATT"
    flanking_R = "ATATCCCATAACACACTGCGTCAGACTTCG"

    new_seq = flanking_F + seq + flanking_R

    new_qual = "I" * len(flanking_F) + qual + "I" * len(flanking_R)
    return new_seq, new_qual

with open("$BACK_REPLACED_FASTQ", "r") as infile, open("temp.fastq", "w") as outfile:
    lines = infile.readlines()
    for i in range(0, len(lines), 4):
        header = lines[i]
        seq = lines[i+1].strip()
        qual = lines[i+3].strip()
        new_seq, new_qual = add_flanking(seq, qual)
        outfile.write(header)
        outfile.write(new_seq + "\n")
        outfile.write("+\n")
        outfile.write(new_qual + "\n")


import os
os.replace("temp.fastq", "$BACK_REPLACED_FASTQ")
EOF

echo "Results:"
echo "  - Intermediate FASTQ file after F cut: $CUTADAPT_F_FASTQ"
echo "  - Intermediate FASTQ file after R cut: $CUTADAPT_FR_FASTQ"
echo "  - Intermediate FASTQ file after RC(F) cut: $CUTADAPT_RC_F_FASTQ"
echo "  - Intermediate FASTQ file after RC(R) cut: $CUTADAPT_RC_FR_FASTQ"
echo "  - Final FASTQ file after merging: $CUTADAPT_FINAL_FASTQ"
echo "  - Annotated FASTQ file: $ANNOTATED_FASTQ (only sequences containing O/E kept)"
echo "  - Back-replaced FASTQ file with flanking sequences: $BACK_REPLACED_FASTQ"

echo "Program ends!"
