#!/usr/bin/env python3


###################################
# script:       SubsetSeq.py
# date:         Mon May 30 13:42:56 CEST 2022
# author:       Abhijeet Singh

###################################
#
version = ': version (1.0)'
program = 'SubsetSeq'
citation = '''\n\nCitation: Singh, Abhijeet. SubsetSeq: a python utility to multisequence file based on 
identifiers from external text file. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.34956.59528, 
Available at GitHub: https://github.com/abhijeetsingh1704/SubsetSeq'''
angle = u'\u2514'
pipe = u'\u2506'
dash = u'\u2500'
L1 = pipe + "\n" + angle + dash
L2 = "\n   " + pipe + "\n   " + angle + dash
program_description = """
Subset multisequence based on identifiers from text file
"""
###################################

##############################################################
import os
import sys
import datetime
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
#
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except  ImportError:
    print("Biopython missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.78"])
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

##############################################################
# Help
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description = program_description + citation)
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--input", dest='input', required=True, type=str, help = L1 + ' input file\n\n')
required.add_argument("-s", "--search", dest='search', required=True, type=str, help = L1 + ' text file with identifiers to query input file\n\n')
required.add_argument("-t", "--type", dest='type', required=True, type=str, help= '\n' + L1 + ' identifier type to match in input file' + L2 + ' 1 = id (sequence identifier before first empty character)' + L2 + ' 2 = description (complete description of sequence according to input file)\n\n')
optional.add_argument("-m", "--match_out", dest='match', required=False, type=str, help = L1 + ' output file with matched sequences (default: <time>-MATCH_<input_file>)\n\n')
optional.add_argument("-u", "--unmatch_out", dest='unmatch', required=False, type=str, help = L1 + ' output file with unmatched sequences (default: <time>-UNMATCH_<input_file>)\n\n')
optional.add_argument("-f", "--input_format", dest='input_format', required=False, default='1', type=int, help = L1 + ' format of input file (default: 1 (fasta))' + L2 + ' 1 = fasta' + L2 + ' 2 = genbank\n\n')
optional.add_argument("-F", "--output_format", dest='output_format', required=False, default='1', type=int, help = L1 + ' format of output file (default: 1 (interleaved fasta))' + L2 + ' 1 = interleaved fasta' + L2 + ' 2 = fasta-2line' + L2 + ' 3 = tab-delimited\n\n')
optional.add_argument("-v", "--verbose", dest='verbose', metavar='Y/y or N/n', default='Y', type=str, help= L1 + ' print progress to the terminal (default: verbose)\n\n')
optional.add_argument('-V', '--version', action='version', version='%(prog)s'+ str(version), help = ' show program\'s version number and exit\n\n')
parser._action_groups.append(optional)
args = parser.parse_args()
##############################################################
#
class colours:
    MATCHheader ='\033[1;37;41m'
    MATCHID = '\033[1;32;40m'
    UNMATCHheader ='\033[1;37;46m'
    UNMATCHID = '\033[1;33;40m'
    RESTORE = '\033[0m'

##############################################################
#
def filenames(filename):
    absname = os.path.abspath(filename)
    basename = os.path.basename(filename)
    return absname, basename

##############################################################
# time
date_var_1 = datetime.datetime.now()
date_var_2 = str(date_var_1.strftime("%Y-%m-%d %H:%M:%S"))
date_var_3 = str(date_var_1.strftime("%H%M%S"))

##############################################################
# mandatory options
input_file = args.input
identifier_fileabs = os.path.abspath(args.search)
search_var = args.type

##############################################################
# optional options match
if args.match is not None:
    match_out_1 = args.match
    matchout = os.path.abspath(match_out_1)
else:
    matchout = os.path.join(os.getcwd(), date_var_3 + "-MATCH_" + filenames(input_file)[1])
    
# optional options unmatch
if args.unmatch is not None:
    unmatch_out_1 = args.unmatch
    unmatchout = os.path.abspath(unmatch_out_1)
else:
    unmatchout = os.path.join(os.getcwd(), date_var_3 + "-UNMATCH_" + filenames(input_file)[1])

# optional options input format
if args.input_format is not None:
    INformat_var = args.input_format
    if type(INformat_var) == int:
        if INformat_var == 1:
            INformat = "fasta"
        if INformat_var == 2:
            INformat = "genbank"
        if INformat_var >= 3:
            sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
    else:
        sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
else:
    INformat = "fasta"

# optional options output format
if args.output_format is not None:
    OUTformat_var = args.output_format
    if type(OUTformat_var) == int:
        if OUTformat_var == 1:
            OUTformat = "fasta"
        if OUTformat_var == 2:
            OUTformat = "fasta-2line"
        if OUTformat_var == 3:
            OUTformat = "tab"
            matchout = os.path.splitext(matchout)[0]+".tab"
            unmatchout = os.path.splitext(unmatchout)[0]+".tab"
        if OUTformat_var > 3:
            sys.exit("[ERROR]: output format flag requires interger value\n\t use -F or --output_format with value 1, 2 or 3")
    else:
        sys.exit("[ERROR]: output format flag requires interger value\n\t use -f or --output_format with value 1, 2 or 3")
else:
    OUTformat = "fasta"

# optional options verbosity
if args.verbose is not None:
    verbosity_var = (args.verbose).lower()
    if verbosity_var == 'n':
        verbosity = "NO"
    elif verbosity_var == 'y':
        verbosity = "YES"
    
##############################################################
#
fasta_list = []
match_list = []
match_out_list = []
unmatch_out_list = []

##############################################################
#
print("-"*80)
print("# [Program]\t: SubsetSeq")
print("# [Date]\t: " + str(date_var_2))
print("-"*80)

##############################################################
#
print("# [Indexing input file]")
with ProcessPoolExecutor() as executor:
    record = SeqIO.index(filenames(input_file)[0], INformat)
    for key, value in record.items():
        if search_var == "1":
            task = fasta_list.append(value.id)
        elif search_var == "2":
            task = fasta_list.append(value.description)
        executor.submit(task)

# ##############################################################
#
print("# [Preparing matching index]")
with ProcessPoolExecutor() as executor:
    with open(identifier_fileabs,"r") as header_file:
        for header in header_file.readlines():
            header = header.replace(">","").strip()
            executor.submit(header)
            task = match_list.append(header)
            executor.submit(task)

# ##############################################################
#
print("# [Preparing unmatching index]")
with ProcessPoolExecutor() as executor:
    unmatch_list = [unmatch for unmatch in fasta_list if unmatch not in match_list]
    executor.submit(unmatch_list)
#
print("-"*80)

##############################################################
#
def def_matching():
    if verbosity == "YES":
        print(f"{colours.MATCHheader}# MATCHING SEQUENCES{colours.RESTORE}")
    with ProcessPoolExecutor() as executor:
        for key, value in record.items():
            if search_var == "1":
                DESCRIPTION = value.id
            elif search_var == "2":
                DESCRIPTION = value.description
            executor.submit(DESCRIPTION)
            for header in match_list:
                if DESCRIPTION == header:
                    Id = value.id
                    Description = value.description
                    Name = value.name
                    Sequence = value.seq
                    if verbosity == "YES":
                        print("# [MATCH]",f"{colours.MATCHID}",Id,f"{colours.RESTORE}")
                    if OUTformat == 'tab':
                        Id = Id + ' ' + Description
                    matching_sequences = (SeqRecord(Seq(str(Sequence)), id=Id, name=Name, description=Description))
                    executor.submit(matching_sequences)
                    task = match_out_list.append(matching_sequences)
                    executor.submit(task)
        if verbosity == "YES":
            print("-"*80)
        #
        task = SeqIO.write(match_out_list, matchout, OUTformat)
        executor.submit(task)

##############################################################
#
def def_unmatching():
    if verbosity == "YES":
        print(f"{colours.UNMATCHheader}# UNMATCHING SEQUENCES{colours.RESTORE}")
    with ProcessPoolExecutor() as executor:
        for key, value in record.items():
            if search_var == "1":
                DESCRIPTION = value.id
            elif search_var == "2":
                DESCRIPTION = value.description
            executor.submit(DESCRIPTION)
            for header in unmatch_list:
                if DESCRIPTION == header:
                    Id = value.id
                    Description = value.description
                    Name = value.name
                    Sequence = value.seq
                    if verbosity == "YES":
                        print("# [UNMATCH]",f"{colours.UNMATCHID}",Id,f"{colours.RESTORE}")
                    if OUTformat == 'tab':
                        Id = Id + ' ' + Description
                    unmatching_sequences = (SeqRecord(Seq(str(Sequence)), id=Id, name=Name, description=Description))
                    unmatch_out_list.append(unmatching_sequences)
        if verbosity == "YES":    
            print("-"*80)
        #
        task = SeqIO.write(unmatch_out_list, unmatchout, OUTformat)
        executor.submit(task)

##############################################################
#
def_matching()
def_unmatching()

##############################################################
#
print(f"{colours.MATCHheader}# Output file =",matchout,f"{colours.RESTORE}")
print(f"{colours.UNMATCHheader}# Output file =",unmatchout,f"{colours.RESTORE}")
print("# Sequences Count: Input file\t= " , (len(fasta_list)))
print("# Sequences Count: match file\t= " , len(match_out_list))
print("# Sequences Count: unmatch file\t= " , len(unmatch_out_list))
print("-"*80+"\n##")

##############################################################
# end of script
