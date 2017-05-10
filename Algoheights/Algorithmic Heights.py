# Define console output Colors
class Color:
    Purple = '\033[95m'
    Cyan = '\033[96m'
    DarkCyan = '\033[36m'
    Blue = '\033[94m'
    Green = '\033[92m'
    Yellow = '\033[93m'
    Red = '\033[91m'
    Bold = '\033[1m'
    Underline = '\033[4m'
    End = '\033[0m'

def readable_print(list, number):
    count = 0
    readout = []
    for i in list:
        readout.append(i)
        count += 1
        if count == number:
            print(readout)
            readout = []
            count = 0
    print('')

def file_len_v2(file):  # more memory efficient?
    for i, l in enumerate(file):
            pass
    return i + 1

def fasta_dict(file):
    sum_lines = file_len_v2(file)

    array = {}

    for i in range(sum_lines):
        y = file.readline()
        if '>' in y:
            k = y.strip('>')
            n = k.strip('\n')
            array[n] = ''
        else:
            b = y.strip('\n')
            array[n] += b
    return array


print(Color.Green + '''Rosalind Armory: prefix arm_
Anaconda:   Installed via windows installer.
            Access Anaconda terminal and type> conda list to see available packages
            *How to identify the import name and convention?
BioPython:  This package can be fetched and installed from the native Anaconda command prompt.
            *How to list the modules in Biopython?
            ''' + Color.End)

print(Color.Blue + '''arm_ini
Given:      A DNA string s of length at most 1000 bp.
Return:     Four integers (separated by spaces) representing the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
Package:    from Bio.Seq import Seq
Conv:       Biopython = Bio
            Seq = Seq, a class that includes many sequence manipulation methods, see help() for algorithm limitations''' + Color.End)
print(Color.Red + '''Questions:  What is class really? How do they work?
            What do double underscore lines mean? i.e. __init__?
            How does 'from Bio.Seq import Seq' differ compared to 'from Bio import Seq'?''' + Color.End)

from Bio.Seq import Seq

with open('/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/rosalind_ini.txt') as y:
    z = y.read()
    f = z.strip('\n')

my_seq=Seq(f)
print(my_seq.count('A'), my_seq.count('C'), my_seq.count('G'), my_seq.count('T'))
print('')

print(Color.Blue + '''arm_dbpr
Given:      The UniProt ID of a protein.
Return:     A list of biological processes in which the protein is involved (biological processes are found in a subsection of the protein's "Gene Ontology" (GO) section).
Package:    from Bio import ExPASy
            from Bio import SwissProt
Conv:       ExPASy = ExPASy
            SwissProt = SwissProt
            ExPASy.get_sprot_raw= returns http used by SwissProt, only good for 1 arg?
            SwissProt.read(x)
            SwissProt.parse(x, y, z)
            dir = use to query SwissProt record''' + Color.End)
print(Color.Red + '''Questions:  What types are the items listed by dir? Are they classes?
            How is the webpage information organised into the format retrieved?''' + Color.End)

from Bio import ExPASy
from Bio import SwissProt

with open('/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/rosalind_dbpr.txt') as y:
    z = y.read()
    f = z.strip('\n')

handle = ExPASy.get_sprot_raw(f)
record = SwissProt.read(handle)

print ('Protein ID: ' +f)
print(handle)
print(record)
print('')
readable_print(dir(record),6)
# now calling a method? class? listed in dir
readable_print(record.cross_references[0:], 3)
# find GO terms?
for i in record.cross_references[0:]:
    if 'GO' in i:
        print(i)
print('')

print(Color.Blue + '''arm_frmt
Given: A collection of n (n≤10  ) GenBank entry IDs.
Return: The shortest of the strings associated with the IDs in FASTA format.
Package:    from Bio import Entrez
Conv:       Entrez = Entrez
            Entrez.efetch''' + Color.End)
print(Color.Red + '''Questions:  I don't really understand the format of the SeqIO parse record''' + Color.End)

with open('/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/rosalind_frmt.txt') as y:
    z = y.read()
    f = z.strip('\n')

from Bio import Entrez
Entrez.email = "jeffsmith8@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=["BT149870, JX469983, NM_001197168, HM595636, JX317624, FJ817486, JX069768, JQ342169, NM_002133"], rettype="fasta")
records = handle.read()
print(records)
frmt_str=''
for i in records.split('\n'):
    if '>' in i:
        frmt_str+='\n'+i+'\n'
    else:
        frmt_str+=i
print (frmt_str)
print('')

# # from Bio import Entrez
from Bio import SeqIO
Entrez.email = "jeffsmith8@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=["BT149870, JX469983, NM_001197168, HM595636, JX317624, FJ817486, JX069768, JQ342169, NM_002133"], rettype="fasta")
records = list (SeqIO.parse(handle, "fasta")) #we get the list of SeqIO objects in FASTA format
print(records)
# print(records.format("fasta"))

for i in records:
    print (i + '\n'*2)
print (readable_print(dir(records[0]),6))
print (records[0].id)  #first record id
frmt_len={}
# for i in

print (len(records[-1].seq))  #length of the last record
print('')
print(Color.Blue + '''arm_meme
Given: A set of protein strings in FASTA format that share some motif with minimum length 20.
Return: Regular expression for the best-scoring motif.
Package:    from Bio import Entrez
Conv:       Entrez = Entrez
            Entrez.efetch''' + Color.End)
print(Color.Red + '''Questions:  The queue time for my job exceeds the Rosalind time''' + Color.End)

print('')
print(Color.Blue + '''arm_need
An online interface to EMBOSS's Needle tool for aligning DNA and RNA strings can be found here.
Use:
            The DNAfull scoring matrix; note that DNAfull uses IUPAC notation for ambiguous nucleotides.
            Gap opening penalty of 10.
            Gap extension penalty of 1.
For our purposes, the "pair" output format will work fine; this format shows the two strings aligned at the bottom of the output file beneath some statistics about the alignment.
Given:      Two GenBank IDs.
Return:     The maximum global alignment score between the DNA strings associated with these IDs. ''' + Color.End)
print(Color.Red + '''Questions:  Can the score be at all meaningful empirically?
            Or is it only meaningful relatively- i.e. for choosing the highest scoring matches?''' + Color.End)

Entrez.email = "jeffsmith8@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=["NM_001251956.1 JQ712981.1"], rettype="fasta")
records = handle.read()
print(records)

print(Color.Blue + '''arm_tfsq
Problem:    Sometimes it's necessary to convert data from FASTQ format to FASTA format. For example, you may want to perform a BLAST search using reads in FASTQ format obtained from your brand new Illumina Genome Analyzer.
Links:      A FASTQ to FASTA converter can be accessed from the Sequence conversion website
            A free GUI converter developed by BlastStation is available here for download or as an add-on to Google Chrome.
            There is a FASTQ to FASTA converter in the Galaxy web platform. Note that you should register in the Galaxy and upload your file prior to using this tool.
            Biopython SeqIO will also do it.
Given:      FASTQ file
Return:     Corresponding FASTA records
Package:    SeqIO.convert(input file path, input format, output file path, output format)
            Notes: SeqIO can take either fastq or txt as an input file but the input format still needs to reflect the convention i.e. fastq
            Output file path and format convention, however, must match ''' + Color.End)
print(Color.Red + '''Questions:  ''' + Color.End)

SeqIO.convert("/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/rosalind_tfsq.fastq", "fastq-illumina", "/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/output-rosalind_tfsq.fasta", "fasta")
print('')

# print(fasta_dict('/media/jeff/Data/MEGA/x- Jeff- Programs & Modules/Rosalind/Algoheights/Algoheights test files/rosalind_ptra.txt'))