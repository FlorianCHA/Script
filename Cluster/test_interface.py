#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package script_pangenome.py
# @author Florian Charriat

import argparse
from gooey import Gooey
 
@Gooey
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input", help="Input File",type=argparse.FileType('r'))
	parser.add_argument("output", help="Output File",type=argparse.FileType('w'))
	parser.add_argument("-t", help="Print time",action='store_true')
	parser.add_argument("--path", help="Path to software",default=softPath)
	parser.parse_args()
