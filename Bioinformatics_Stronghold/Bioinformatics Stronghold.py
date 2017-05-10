# Import packages
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import ExPASy
from Bio import SwissProt
from Bio import Entrez
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import itertools
from math import factorial


# Define OS environment
windows='D:/'
linux='/media/jeff/Data/'
project_directory='MEGA/x- Jeff- Programs & Modules/Rosalind/Bioinformatics Stronghold/str test files/'
os=windows

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

# Define special variables
dna_aa_table = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'],
                'M': ['ATG'],
                'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'],
                '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
                'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
                'C': ['TGT', 'TGC'], 'W': ['TGG'],
                'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']}

rna_aa_table = {'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'I': ['AUU', 'AUC', 'AUA'],
                'M': ['AUG'],
                'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                'U': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'],
                '*': ['UAA', 'UAG', 'UGA'], 'H': ['CAU', 'CAC'],
                'Q': ['CAA', 'CAG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'],
                'C': ['UGU', 'UGC'], 'W': ['UGG'],
                'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']}

# Define special functions
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

def txt_file_len(file):  # more memory efficient?
    with open (file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def multiline_txt_to_fasta_dict(file):
    input_file = os + project_directory + file

    sum_lines = txt_file_len(input_file)

    array = {}
    with open(input_file) as f:
        for i in range(sum_lines):
            y = f.readline()
            if '>' in y:
                k = y.strip('>')
                n = k.strip('\n')
                array[n] = ''
            else:
                b = y.strip('\n')
                array[n] += b
        return array

# str_prtn
# Problem
# To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows:
# [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.
#
# You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into
# http://www.uniprot.org/uniprot/uniprot_id
# Alternatively, you can obtain a protein sequence in FASTA format by following
# http://www.uniprot.org/uniprot/uniprot_id.fasta
# For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.
#
# Given: At most 15 UniProt Protein Database access IDs.
# Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.
#
# Sample Dataset
# A2Z669
# B5ZC00
# P07204_TRBM_HUMAN
# P20840_SAG1_YEAST
#
# Sample Output
# B5ZC00
# 85 118 142 306 395
# P07204_TRBM_HUMAN
# 47 115 116 382 409
# P20840_SAG1_YEAST
# 79 109 135 248 306 348 364 402 485 501 614

def n_glycosylation_motifs(string):
    positions=[]
    for i in range(len(string)-3):
        if string[i] =='N' and string[i+1]!='P' and (string[i+2]=='S' or string[i+2]=='T') and string[i+3]!='P':
            positions.append(i+1)
    return positions

mprt_list = []
sum_lines = txt_file_len(os + project_directory + 'rosalind_mprt.txt')

with open(os + project_directory + 'rosalind_mprt.txt') as f:
    for i in range(sum_lines):
        y = f.readline()
        n = y.strip('\n')
        if n!='':
            mprt_list.append(n)

# for i in mprt_list:
#     handle = ExPASy.get_sprot_raw(i)
#     record = SwissProt.read(handle)
#     i_seq = record.sequence
#     glyc_list = n_glycosylation_motifs(i_seq)
#     print (i)
#     print (*glyc_list, sep=' ')

# str_lia
# Problem
# Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
# Return: The adjacency list corresponding to O3. You may return edges in any order.

lia_dict=multiline_txt_to_fasta_dict('rosalind_grph.txt')

def adjacency_list(dict, overlap):
    adjacency_list=[]
    for i in dict:
        prefix=dict[i][:overlap]
        suffix=dict[i][-overlap:]
        for k in dict:
            if prefix==dict[k][-overlap:] and i!=k:
                x=(k + ' '+ i)
                if x not in adjacency_list:
                    adjacency_list.append(x)
            if suffix==dict[k][:overlap] and i!=k:
                y=(i + ' ' +k)
                if y not in adjacency_list:
                    adjacency_list.append(y)
    return adjacency_list

str_lia=adjacency_list(lia_dict,3)
print(str_lia)
print(*str_lia, sep='\n')

# str_splc
# Problem
# Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
# Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)

splc_dict=multiline_txt_to_fasta_dict('rosalind_splc.txt')
# print(splc_dict)

def splice_premrna(dict):
    mrna_list = []
    for i in splc_dict:
        mrna_list.append(splc_dict[i])

    premrna = max(mrna_list, key=len)
    for i in mrna_list:
        if i in premrna and i != premrna:
            premrna = premrna.replace(i, '')
    return(premrna)

mrna=splice_premrna('rosalind_splc.txt')
# print(mrna)

def ros_prot_all_dnaframes(Text):
    protein = []
    for i in range(len(Text) - 2):
        for key, val in dna_aa_table.items():
            if Text[i:i + 3] in val:
                protein.append(key)
    return protein

def ros_prot_dnaframe(Text, k):
    prot = ros_prot_all_dnaframes(Text)
    prot1 = prot[k::3]
    protein = ''.join(prot1)
    return protein

# print(ros_prot_dnaframe(mrna,0))

# str_perm
# Problem
# Given: A positive integer n≤7.
# Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).

def str_perm(number, perm_len):
    integer_range=[]
    for i in range(number):
        integer_range.append(i+1)

    perm_list=list(itertools.permutations(integer_range, perm_len))
    return perm_list

# perm=str_perm(6,6)
# perm_string=''
# for i in range(len(perm)):
#     perm_string+='\n'
#     for j in range(len(perm[i])):
#         perm_string+=str(perm[i][j])+' '
# print(len(perm))
# print(perm_string)

# str_revp
# Problem:
# Given: A DNA string of length at most 1 kbp in FASTA format.
# Return: The position and length of every reverse palindrome in the string having length between 4 and 12.
# You may return these pairs in any order.

def complement(Nucleotide):
    comp=''
    for i in Nucleotide:
        if i=='A':
            comp+='T'
        elif i=='T':
            comp += 'A'
        elif i=='G':
            comp += 'C'
        elif i=='C':
            comp += 'G'
        else:
            comp += 'n'
    return comp

def reverse(Pattern):
    result=''
    index=len(Pattern)-1
    while index>=0:
        result+=Pattern[index]
        index-=1
    return result

def reverse_complement(string):
    revcomp=complement(reverse(string))
    return revcomp

def find_palindrome(string, k):
    pal_find={}
    for i in range (len(string)-k+1):
        pattern=string[i:i+k]
        if pattern==reverse_complement(pattern):
            pal_find[i+1]=k
    return pal_find

def str_revc(string, min, max):
    pattern_list=[]
    for i in range (max-min+1):
        pattern_list.append(min+i)

    for i in pattern_list:
        x=find_palindrome(string, i)
        for keys, values in x.items():
            print (str(keys) + ' ' + str(values))

# str_revp_dict=multiline_txt_to_fasta_dict('rosalind_revp.txt')
# print(str_revp_dict)
#
# str_input='''TATATA '''
# str_revc(str_input, 4, 12)

# str_pper
# Problem
# Given: Positive integers n and k such that 100≥n>0 and 10≥k>0.
# Return: The total number of partial permutations P(n,k), modulo 1,000,000.

# Note that ORDER COUNTS when calculating the number of permutaations of a set vs the number of combinations of a set.
# See http://stattrek.com/online-calculator/combinations-permutations.aspx for a good explanation

def str_pper(n, r):
    return (factorial(n)/factorial(n-r))%1000000

# str_pmch
# Problem: http://rosalind.info/problems/pmch/
# Given: An RNA string s of length at most 80 bp having the same number of occurrences of 'A' as 'U' and the
# same number of occurrences of 'C' as 'G'.
# Return: The total possible number of perfect matchings of basepair edges in the bonding graph of s.

pmch_test=multiline_txt_to_fasta_dict('rosalind_pmch.txt')
for key, values in pmch_test.items():
    pmch_count=values
pmch_a=pmch_count.count('A')
pmch_c=pmch_count.count('C')

str_pmch = factorial(pmch_a)*factorial(pmch_c)
# print(str_pmch)
# print(factorial(pmch_a))
# print(factorial(pmch_c))

with open(os + project_directory + 'output-rosalind_pmch.txt', 'w') as output_data:
    output_data.write(str(str_pmch))

# str_tree
# Given: A positive integer n (n≤1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.
# Return: The minimum number of edges that can be added to the graph to produce a tree.

def multiline_txt_to_list(file):
    input_file=os + project_directory + file

    sum_lines = txt_file_len(input_file)

    list=[]
    with open (input_file) as f:
        for i in range(sum_lines):
            y = f.readline()
            k=y.strip('\n')
            list.append(k)
        return list

def multiline_txt_to_listof_intlists(file):
    input_file=os + project_directory + file

    sum_lines = txt_file_len(input_file)

    list=[]
    with open (input_file) as f:
        for i in range(sum_lines):
            m=[]
            y = f.readline()
            k=y.strip('\n')
            j=k.split(' ')
            for n in j:
                a=int(n)
                m.append(a)
            list.append(m)
        return list

# x=multiline_txt_to_listof_intlists('rosalind_tree.txt')
x=multiline_txt_to_list('rosalind_tree.txt')
print(x, 'hi')

def tree_missing_nodes(list):
    merged_list = []
    for i in range(1,len(list)):
        merged_list+=(list[i])

    missing_nodes=[]
    for i in range(1,1+list[0][0]):
        if i not in merged_list:
            missing_nodes.append(i)
    return missing_nodes

with open (os+project_directory+'rosalind_tree.txt') as f:
    edges=f.read().strip().split('\n')
    n=int(edges.pop(0))
    edges=[map(int,edge.split())for edge in edges]

connected_nodes=[{i} for i in range(1,n+1)]
print('edges: ', list(edges))
print('connected nodes: ',connected_nodes)

# for i in edges:
#     print (list(map(edges(i))))

# for edge in edges:
#     temp_nodes=set()
#     del_nodes=[]
#     for nodes in connected_nodes:


#     adj_list=list[1:]
#     node_groups=[adj_list[0]]
#     first_nodes=[]
#     second_nodes=[]
#     for i in range(len(adj_list)):
#         first_nodes.append(adj_list[i][0])
#         second_nodes.append(adj_list[i][1])
#     print(adj_list)
#     print(first_nodes)
#     print(second_nodes)
#
#     for i in range(len(first_nodes)):
#         x=first_nodes[i]
#         matched_pairs = []
#         for k in range(len(adj_list)):
#             if x in adj_list[k]:
#                 matched_pairs.append(adj_list[k])
#
#             for k in range(0,len(node_groups)):
#                 if x in node_groups[k]:
#                     node_groups[k].append(x)
#                 else:
#                     node_groups.append(x)
#     return node_groups

# def node_groups(list):
#     adj_list=list[1:]
#     node_groups=[[]]
#     first_nodes=[]
#     second_nodes=[]
#     for i in range(len(adj_list)):
#         first_nodes.append(adj_list[i][0])
#         second_nodes.append(adj_list[i][1])
#
#     for i in range(len(adj_list)):
#         x=adj_list[i][0]
#         y=adj_list[i][1]
#         int_list=[]
#         for m in range(len(first_nodes)):
#             if x==first_nodes[m] or y==first_nodes[m]:
#                 int_list.append(adj_list[m])
#         for n in range(len(second_nodes)):
#             if x==second_nodes[n] or y==second_nodes[n]:
#                 int_list.append(adj_list[n])
#         return int_list


# def node_group(list):
#     adj_list = list[1:]
#     node_index=[]
#     node_group=[]
#
    # for i in range(len(adj_list)):
    #     x=adj_list[i][0]
    #     y=adj_list[i][1]
#         new_list=adj_list[:i]+adj_list[i+1:]
#         if adj_list[i] not in node_index:
#             node_index.append(adj_list[i])
#             node_group.append([adj_list[i]])
#         for k in new_list:
#             if x or y in k:


# def node_group(list):
#     adj_list = list[1:]
#     node_index=[]
#     node_group={}
#
#     for i in range(len(adj_list)):
#         x=adj_list[i][0]
#         y=adj_list[i][1]
#         new_list = adj_list[:i] + adj_list[i + 1:]
#         if adj_list[i] not in node_index:
#             node_index.append(adj_list[i])
#             node_group[i]=[adj_list[i]]
#         for k in new_list:
#             if x or y in k:
#                 node_group()
#

# tree_missing_nodes(x)

# print(node_groups(x))