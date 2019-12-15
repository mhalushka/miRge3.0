#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import gzip
from cutadapt.modifiers import AdapterCutter, QualityTrimmer, UnconditionalCutter, QualityTrimmer
import cutadapt

#Custom miRge libraries 
from libs.parse import parseArg
from libs.miRgeEssential import check_dependencies

#with gzip.open('/home/joe/file.txt.gz', 'rb') as f:

def main():
    #WORKING 
    args = parseArg()
    check_dependencies(args)
    
    #TESTING 
    samples = args.samples
    file_exts = ['.txt', '.csv']
    file_list = samples[0].split(',')
    fastq_fullPath = [] #Initiating an empty list to collect all the fastq or fastq.gz files based on user arguments!
    if Path(file_list[0]).exists() and Path(file_list[0]).suffix in file_exts: # READ TXT OR CSV FILE HERE
        with open(file_list[0]) as file:
            for lines in file: 
                line = lines.strip()
                filetype = ''.join(Path(line).suffixes) if Path(line).suffix == ".gz" else Path(line).suffix
                if Path(line).exists() and (filetype == ".fastq" or filetype == ".fastq.gz"):
                    #print(f"{line} and {filetype}")
                    fastq_fullPath.append(str(Path(line).resolve()))
                else:
                    print(f"\nWARNING: File {line} does not exists!") if not Path(line).exists() else print(f"\nWARNING: File {line} is neither fastq or fastq.gz format!")
                    print(f"Omitting file {line}")
    else:
        for files in file_list: # READ FASTQ OR FASTQ.gz FILES HERE
            filetype = ''.join(Path(files).suffixes) if Path(files).suffix == ".gz" else Path(files).suffix
            if Path(files).exists() and (filetype == ".fastq" or filetype == ".fastq.gz"):
                fastq_fullPath.append(str(Path(files).resolve()))
            else:
                print(f"\nWARNING: File {files} does not exists!") if not Path(files).exists() else print(f"\nWARNING: File {files} is neither fastq or fastq.gz format!")
                print(f"Omitting file {files}")
    #Validating files in the list where the list should not be empty: 
    if not fastq_fullPath:
        print("\nERROR!: No valid input files were available!\nPlease verify miRge -s arguments")
        exit()


if __name__ == '__main__':
    main()
    
