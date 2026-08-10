import argparse
import itertools
import pandas as pd
import matplotlib.pyplot as plt

##methonds: count_records -> int, validate_fastq_batched -> list, validate_record -> bool

def count_records(path) -> int:
    count = 0
    with open(path) as file:
        for line in file:
            count += 1
    if count % 4 == 0:
        return count // 4
    
    raise ValueError(f"File has {count} lines, not divisible by 4")

def validate_fastq_batched(path):
    count_records(path)
    with open(path, 'r') as file:
        broken_records = []
        for i, record in enumerate(itertools.batched(file, 4)):
            record_list = []
            for line in record:
                record_list.append(line)
            if not validate_record(record_list):
                broken_records.append(i)
    return broken_records
            

def validate_record(lines: list[str]) -> bool:
    if (lines[0][0] != "@" 
        or len(lines[1]) != len(lines[3]) 
        or lines[2][0] != "+"):
        return False
    return True


def qc_summary(path) -> dict:
    """"
    a function that gives "mini fastQC stats to determine quality of read.
    It will compute num_reads, min_len, max_len, mean_len, gc_content

    Why use it: does not keep reads in memmory unlike read_stats(path)
    """
    qc_summary_stats = {}

    qc_summary_stats["num_reads"] = count_records(path)
    
    with open(path, 'r') as file:
        min_len = float('inf')
        max_len = -float('inf')
        total_len = 0
        G_content, C_content = 0, 0
        for record in itertools.batched(file, 4):
            current_seq = record[1].strip()
            total_len += len(current_seq.strip())
            min_len = min(len(current_seq), min_len)
            max_len = max(len(current_seq,), max_len)

            #now lets count GC content:
            for nuc in current_seq.lower():
                if nuc == "g":
                    G_content += 1
                if nuc == "c":
                    C_content += 1

    qc_summary_stats["min_len"] = min_len
    qc_summary_stats["max_len"] = max_len

    mean_len = total_len / qc_summary_stats["num_reads"]
    qc_summary_stats["mean_len"] = mean_len

    if total_len == 0:
        qc_summary_stats["gc_content"] = None
    else:
        qc_summary_stats["gc_content"] = round((G_content + C_content) * 100 / total_len, 3)

    return qc_summary_stats

def read_stats(path):
    #returns a dataframe of: lengths, gc content
    with open(path, 'r') as file:
        rows = []
        for record in itertools.batched(file, 4):
            length_record = len(record[1].strip())

            #count gc
            g_content, c_content = 0, 0
            for bp in record[1].strip().lower():
                if bp == "g":
                    g_content += 1
                elif bp == "c":
                    c_content += 1
            gc_content = round(((g_content + c_content) *100 / length_record), 3)

            #add the stats as doctionary to rows
            rows.append({"length": length_record, "gc_content": gc_content})

        df = pd.DataFrame(rows)
        return df

def plot_gc_content(df):
    ax = df["gc_content"].hist(bins=10)
    ax.set_xlabel("GC content %")
    ax.set_ylabel("Number of reads")
    ax.set_title("GC Content")
    plt.savefig("plots/gc_content.png")
    plt.close()

def plot_length_distribution(df):
    ax = df["length"].hist(bins=30)
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.set_title("Reads Length Distribution")
    plt.savefig("plots/length_distribution.png")
    plt.close()



def main():
    parser = argparse.ArgumentParser(description="Count FASTQ records")
    parser.add_argument("path", help="Path to fastq file")
    #when running the python file, calling the function, the first argument
    #will be the path (from the CL)
    args = parser.parse_args()
    path = args.path
    
    #now I can call a function with the argument from argparse:

    #validate_fastq_assume_starts_with_at(path)
    #validate_fastq_assume_4_lines_record(path)
    #print(f"number of records:{count_records(path)}")
    #print(f"broken records:{validate_fastq_batched(path)}")
    print(qc_summary(path))
    df = read_stats(path)
    #print(df.head())
    #print(df.describe())
    plot_gc_content(df)
    plot_length_distribution(df)


if __name__ == "__main__": #dunder milfin
    main()
    