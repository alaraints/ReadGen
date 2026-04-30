# ReadGen
fastq read simulator
## License

This software is dual-licensed:

- Open-source under the GNU GPL v3 for academic and non-commercial use
- Commercial license required for proprietary, industrial, or SaaS usage

For commercial licensing, contact: Alar.Aints@gmail.com


ReadGen © Alar Aints 2026
a program for generating simulated fastq sequence reads

Options: 

	-l Read Length  [100]
	-n number of reads to generate [1000000]
	-m mutation probability % per base [0]
	-o outfiles (optional)
	-p paired output, no argument. Single by default
	-O where the outfiles go (Optional. By default the outfiles go to the same folder as input files.)
	-t title tag for title lines (optional)
	-L readlimit. 0 reads the whole files [0]
	-i input file(s). Many files can be entered, gap separated. Fasta format only
	-F only forward reads output
	-R only reverse reads output
	-B output reads both ways mixed  
 

System requirements

Disk space: 
using default settings, one outfile with 1M reads of 100 bp takes about 250 MB.

Memory:
Largest infile size + 1 MB

Dependencies: 
ReadGen has no dependencies. ReadGen is written in C.

File formats:
Input: fasta format only. Output: fastq format. Single or paired reads can be generated, -F forward, -R reverse, -B mixed orientation. 

MacOS notice:
ReadGen runs on MacOS and Linux. 

Network connectivity:
ReadGen has no network connectivity functions. All files must be available locally. 

Multihreading: 
ReadGen runs on a single thread

Behaviour
Input files: -i FOLDER/*.* can be used. Only fasta files will be processed. 
Output files: -o none, one or two files can be specified. Outfile names will be completed automatically. Outfiles are written in fastq format. If -p is set, two output files with paired reads will be generated. If -t is specified, read names will contain this title tag. If not specifed, title tag is derived from filename. Output path -O can be specified, outfiles will be written in this folder. If the folder does not exist, it will be created. Number of reads -n applies per infile. Readlimit -L applies per infile. Output from successive infiles will be appended to the outfiles. 
