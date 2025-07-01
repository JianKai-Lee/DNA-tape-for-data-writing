import os
from math import floor

def classify_reads(input_fastq):

    stats = {}


    with open(input_fastq, "r") as infile:
        lines = infile.readlines()
        for i in range(0, len(lines), 4):
            header = lines[i].strip()
            seq = lines[i+1].strip()
            qual = lines[i+3].strip()


            seq_len = len(seq)


            K = round((seq_len - 70) / 15)
            output_file = f"B{K}.fastq"


            if output_file not in stats:
                stats[output_file] = 0


            with open(output_file, "a") as outfile:
                outfile.write(header + "\n")
                outfile.write(seq + "\n")
                outfile.write("+\n")
                outfile.write(qual + "\n")


            stats[output_file] += 1


    with open("read_stats.txt", "w") as statfile:
        for file, count in stats.items():
            statfile.write(f"{file}: {count} reads\n")


if __name__ == "__main__":
    input_fastq = "Back_Replaced.fastq"
    classify_reads(input_fastq)
    print("Sorted! Stores in read_stats.txt。")
