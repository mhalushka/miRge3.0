#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import gzip

#Custom miRge libraries 
from libs.parse import parseArg
from libs.miRgeEssential import check_dependencies

#with gzip.open('/home/joe/file.txt.gz', 'rb') as f:

args = parseArg()
check_dependencies(args)

samples = args.samples

print(Path(samples).suffix)
