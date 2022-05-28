# SubsetSeq
a python utility to subset multisequence file based on identifiers from external text file




```
usage: SubsetSeq [-h] -i INPUT -s SEARCH -t TYPE [-m MATCH] [-u UNMATCH]
                 [-f INPUT_FORMAT] [-F OUTPUT_FORMAT] [-v Y/y or N/n] [-V]

Subset multifasta based on identifiers from text file

Citation: Singh, Abhijeet. SubsetSeq: a python utility to multisequence file based on 
identifiers from external text file. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.34956.59528, 
Available at GitHub: https://github.com/abhijeetsingh1704/SubsetSeq

required arguments:
  -i INPUT, --input INPUT
                        ┆
                        └─ input file
                        
  -s SEARCH, --search SEARCH
                        ┆
                        └─ text file with identifiers to query input file
                        
  -t TYPE, --type TYPE  
                        ┆
                        └─ identifier type to match in input file
                           ┆
                           └─ 1 = id (sequence identifier before first empty character)
                           ┆
                           └─ 2 = description (complete description of sequence according to input file)
                        

options:
  -h, --help            show this help message and exit
  -m MATCH, --match_out MATCH
                        ┆
                        └─ output file with matched sequences (default: <time>-MATCH_<input_file>)
                        
  -u UNMATCH, --unmatch_out UNMATCH
                        ┆
                        └─ output file with unmatched sequences (default: <time>-UNMATCH_<input_file>)
                        
  -f INPUT_FORMAT, --input_format INPUT_FORMAT
                        ┆
                        └─ format of input file (default: 1 (fasta))
                           ┆
                           └─ 1 = fasta
                           ┆
                           └─ 2 = genbank
                        
  -F OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        ┆
                        └─ format of output file (default: 1 (interleaved fasta))
                           ┆
                           └─ 1 = interleaved fasta
                           ┆
                           └─ 2 = fasta-2line
                           ┆
                           └─ 3 = tab-delimited
                        
  -v Y/y or N/n, --verbose Y/y or N/n
                        ┆
                        └─ print progress to the terminal (default: verbose)
                        
  -V, --version          show program's version number and exit
                        
```

# Citation

Singh, Abhijeet. SubsetSeq: a python utility to multisequence file based on 
identifiers from external text file. ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.34956.59528, 
Available at GitHub: https://github.com/abhijeetsingh1704/SubsetSeq
