# DNA-tape-for-data-writing

## `1. block_replacement`

`block_replacement.sh` filters and cleans raw nanopore FASTQ reads, **retaining only sequences that satisfy the project’s coding rules**.  
It automatically

* detects **odd/even structural blocks** and the **payload** in every read;
* corrects any errors located in those structural blocks,

so that **all remaining sequencing errors reside solely in the payload region**.  
Because of this correction step, **this sctipt is not used for whole‐read error-rate statistics**.

---

###  1.1 Preparing the input

1. **Rename** the sequencing FASTQ file to `merged.fastq`

2. If a single run produced *multiple* FASTQ files, concatenate them first, then save the merged file as `merged.fastq`.

---

###  1.2 Running the script
Run the following commands in your terminal:

```bash
cd /path/to/working_directory     # 1 Navigate to the folder
chmod +x block_replacement.sh     # 2 Make the script executable
./block_replacement.sh            # 3 Start the analysis
```
---

## `2. length_filter`

This script further **classifies the cleaned reads by payload length** and writes each group to its own FASTQ file: 

`Bk.fastq   # where k = payload length in blocks`

---

### 2.1 Usage

```bash
python3 length_filter.py
```
Running the command inside the working directory will generate one file per unique payload length (e.g., B21.fastq, B25.fastq, …).


