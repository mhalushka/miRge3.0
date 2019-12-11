#!/usr/bin/env python


# Import miRge libraries 
from libs.parse import parseArg
from libs.miRgeEssential import check_dependencies

args = parseArg()
check_dependencies(args)
