import argparse

def count_records(path) -> int:
    count = 0
    with open(path) as file:
        for line in file:
            count += 1
    return count // 4

def main():
    parser = argparse.ArgumentParser(description="Count FASTQ records")
    parser.add_argument("path", help="Path to fastq file")
    #when running the python file, calling the function, the first argument
    #will be the path (from the CL)
    args = parser.parse_args()
    print(count_records(args.path))
    #now I call the function with the argument from argparse

if __name__ == "__main__":
    main()