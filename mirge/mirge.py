#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import gzip

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
    if Path(file_list[0]).suffix in file_exts: # READ TXT OR CSV FILE HERE
        val = Path(file_list[0]).exists()
        #val = Path(file_list[0]).is_file()
        print(val)
        pass
    else:
        for files in file_list: # READ FASTQ OR FASTQ.gz FILES HERE
            #print(Path(samples[0]).suffix)
            print(Path(files).suffix)



if __name__ == '__main__':
    main()
    
