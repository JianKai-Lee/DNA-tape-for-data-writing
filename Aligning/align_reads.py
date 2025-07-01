# -*- coding: utf-8 -*-

from Bio import SeqIO
import os

def read_fasta(fasta_file):

    ref_sequences = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_sequences[record.id] = str(record.seq).replace("\n", "")
    return ref_sequences

def process_fastq(fastq_file, ref_sequences, output_dir):

    os.makedirs(output_dir, exist_ok=True)

    failed_file = os.path.join(output_dir, "failed.fastq")
    failed_handle = open(failed_file, "w")

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            matched_flag = False

            for ref_id, ref_seq in ref_sequences.items():
                if seq == ref_seq:

                    output_file = os.path.join(output_dir, f"{ref_id}.fastq")
                    with open(output_file, "a") as output_handle: 
                        SeqIO.write(record, output_handle, "fastq")
                    matched_flag = True
                    break 

            if not matched_flag:
                SeqIO.write(record, failed_handle, "fastq")

    failed_handle.close()

def count_sequences_in_fastq(file_path):

    count = 0
    with open(file_path, "r") as handle:
        for _ in SeqIO.parse(handle, "fastq"):
            count += 1
    return count

def write_statistics(output_dir, input_count):

    stats_file = os.path.join(output_dir, "statistics.txt")
    matched_counts = {}
    failed_count = 0

    for file_name in os.listdir(output_dir):
        if file_name.endswith(".fastq") and file_name != "failed.fastq":
            ref_id = file_name.replace(".fastq", "")
            file_path = os.path.join(output_dir, file_name)
            matched_counts[ref_id] = count_sequences_in_fastq(file_path)

    failed_file = os.path.join(output_dir, "failed.fastq")
    if os.path.exists(failed_file):
        failed_count = count_sequences_in_fastq(failed_file)

    with open(stats_file, "w") as stats_handle:
        stats_handle.write(f"Input reads: {input_count}\n")
        stats_handle.write("\nMatched reads:\n")
        for ref_id, count in matched_counts.items():
            stats_handle.write(f"{ref_id}: {count}\n")
        stats_handle.write(f"\nUnmatched reads: {failed_count}\n")

        total_matched = sum(matched_counts.values())
        total_output = total_matched + failed_count
        stats_handle.write(f"\nMatched reads + Unmatched reads: {total_output}\n")
        if total_output == input_count:
            stats_handle.write("Check: No repetitions or omissions.\n")
        else:
            stats_handle.write("Check: Warning! Repetitions or omissions found.\n")

def main():
    fasta_file = "ref.fasta"  # Replace with your input reference file
    fastq_file = "B4.fastq"  # Replace with your input Bk.fastq file
    output_dir = "output"  # Output folder name

    ref_sequences = read_fasta(fasta_file)
    
    process_fastq(fastq_file, ref_sequences, output_dir)

    input_count = count_sequences_in_fastq(fastq_file)

    write_statistics(output_dir, input_count)

if __name__ == "__main__":
    main()
