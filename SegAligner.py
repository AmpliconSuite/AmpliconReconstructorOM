import os
import sys
import time
import bisect
import argparse
import subprocess

def run_SegAligner(contig_file,ref_file,arg_list):
    argstring = " ".join(arg_list)
    subprocess.call(" ".join([os.environ['SA_SRC'] + "/SegAligner",ref_file,contig_file,argstring]), shell=True)