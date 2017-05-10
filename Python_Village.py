import numpy as np


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


print(Color.Blue + '''ini2
Given: Two positive integers a and b, each less than 1000.
Return: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths a and b.''' + Color.End)


def sqr_hypotenuse(a, b):
    hyp = (a ** 2) + (b ** 2)
    return hyp


print(sqr_hypotenuse(805, 933))
print('')

print(Color.Blue + '''#ini3
#Given: A string of length at most 200 letters and four integers a, b, c and d
#Return: The slice of this string from indices through a to b and c through d (with space in between), inclusively.''' + Color.End)


def slice_string(string, a, b, c, d):
    ab = string[a:b + 1]
    cd = string[c:d + 1]
    return '%s %s' % (ab, cd)


print(slice_string('ihavegoatlolly', 0, 0, 5, 8))
print('')

print(Color.Blue + '''ini4
Given: Two positive integers a and b (a<b<10000)
Return: The sum of all odd integers from a through b, inclusively.''' + Color.End)


def sum_odd(a, b):
    sum_a = 0
    for i in range(a, b + 1):
        if i % 2 != 0:
            sum_a += i
    return sum_a


print(sum_odd(4877, 9292))
print('')

print(Color.Cyan + '''mycode
Given: A file to read
Return: Read contents from default path D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/''' + Color.End)


def read(file):
    path = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    f = open(path, 'r')
    print(f.read())
    f.close()


# read('rosalind_ini5.txt')
print('')

print(Color.Cyan + '''mycode
Given: A file to open or close
Return: Automatically open a file from the default path D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/''' + Color.End)


def access(file):
    path = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    f = open(path, 'r')
    return f


print('')

print(Color.Cyan + '''mycode
Given: A file
Return: Count lines from file in default path D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/''' + Color.End)


def file_len_v1(file):  # memory inefficient?
    f = access(file)
    x = len(f.readlines())
    f.close()
    return x


print(file_len_v1('rosalind_ini5.txt'))
print('')

print(Color.Cyan + '''mycode- alternative
Given: A file
Return: Count lines in file from default path D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/''' + Color.End)


def file_len_v2(file):  # more memory efficient?
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


print(file_len_v2('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/rosalind_ini5.txt'))
print('')

print(Color.Blue + '''ini5
Given : A file containing at most 1000 lines.
Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.''' + Color.End)


def ini5(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file
    x = file_len_v2(input_file)
    y = open(output_file, 'w')
    for i in range(x):
        if i % 2 != 0:
            with open(input_file, 'r') as f:
                y.write(f.readlines()[i])
    y.close()


ini5('rosalind_ini5.txt')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_ini5.txt') as g:
    print(g.read())
print('')

print(Color.Blue + '''ini6
Given: A string of length at most 10000 letters.
Return: How many times any word occurred in string. Each letter case (upper or lower) in word matters.
Lines in output can be in any order.''' + Color.End)


def dict_count(string_input):
    ls = list(string_input.split(' '))
    x = {}
    for i in ls:
        x[i] = ls.count(i)
    return x


def dict_col(string_input):
    x = dict_count(string_input)
    l = ''
    for key, value in x.items():
        l += key + ' ' + str(value) + '\n'
    return l


def ini6(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    l = dict_col(t)

    f = open(output_file, 'w')
    f.write(l)
    f.close()


ini6('rosalind_ini6.txt')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_ini6.txt') as g:
    print(g.read())

print(Color.Blue + '''ros_dna
Given: A string of length at most 10000 letters.
Return: How many times any word occurred in string. Each letter case (upper or lower) in word matters.
Lines in output can be in any order.''' + Color.End)


def ros_dna(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    n = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for i in t:
        n[i] += 1

    return n
    # tally=N_count(t)          #replace last line of ros_DNA with these two for printed string output ordered ACGT.
    # print(tally['A'], tally['C'], tally['G'], tally['T'])


print(ros_dna('rosalind_dna.txt'))
print('')

print(Color.Blue + '''ros_rna
Given: A DNA string t having length at most 1000 nt.
Return: The transcribed RNA string of t.''' + Color.End)


def ros_rna(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    rna = t.replace('T', 'U')

    f = open(output_file, 'w')
    f.write(rna)
    f.close()


ros_rna('rosalind_rna.txt')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_rna.txt') as g:
    print(g.read())
print('')

print(Color.Blue + '''ros_complement
Given: A DNA string s of length at most 1000 bp.
Return: The reverse complement sc of s.''' + Color.End)


def ros_revc(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    rev_t = t[::-1]

    comp = ''
    for i in rev_t:
        if i == 'A':
            comp += 'T'
        elif i == 'T':
            comp += 'A'
        elif i == 'G':
            comp += 'C'
        elif i == 'C':
            comp += 'G'
        else:
            comp += 'n'

    f = open(output_file, 'w')
    f.write(comp)
    f.close()


ros_revc('rosalind_revc.txt')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_revc.txt') as g:
    print(g.read())
print('')

print(Color.Blue + '''ros_iprb
Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms:
k individuals are homozygous dominant for a factor, m are heterozygous, and n  are homozygous recessive.
Return: The probability that two randomly selected mating organisms will produce an individual possessing
a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.''' + Color.End)


def ros_iprb(k, m, n):
    pop = k + m + n
    pop1 = pop - 1
    pk = k / pop
    pm = m / pop
    pn = n / pop
    p1k = k / pop1
    p1m = m / pop1
    p1n = n / pop1
    p1k1 = (k - 1) / pop1
    p1m1 = (m - 1) / pop1
    p1n1 = (n - 1) / pop1
    p = np.array([[1, 0, 0], [.5, .5, 0], [0, 1, 0], [.5, .5, 0], [.25, .5, .25], [0, .5, .5], [0, 1, 0], [0, .5, .5],
                  [0, 0, 1]])
    h = np.array(
        [[pk * p1k1], [pk * p1m], [pk * p1n], [pm * p1k], [pm * p1m1], [pm * p1n], [pn * p1k], [pn * p1m], [pn * p1n1]])
    pr_homo_dom = h[0][0] * p[0][0] + h[1][0] * p[1][0] + h[3][0] * p[3][0] + h[4][0] * p[4][0]
    pr_hetero = h[1][0] * p[1][1] + h[2][0] * p[2][1] + h[3][0] * p[3][1] + h[4][0] * p[4][1] + h[5][0] * p[5][1] + \
                h[6][0] * p[6][1] + h[7][0] * p[7][1]
    pr_homo_rec = h[4][0] * p[4][2] + h[5][0] * p[5][2] + h[7][0] * p[7][2] + h[8][0] * p[8][2]
    phd = 'Pr[homozygous dominant]= %s' % pr_homo_dom
    ph = 'Pr[heterozygous]= %s' % pr_hetero
    phr = 'Pr[homozygous recessive]= %s' % pr_homo_rec
    print(
        '''Array "h" (below) represents the Pr(Parent 1 x Parent 2) of each possible pairing in order: \nGG.GG, GG.Gg, GG.gg, Gg.GG, Gg.Gg, Gg.gg, gg.GG, gg.Gg, gg.gg \nwith respect to k, m, n values''')
    print(h)
    print(
        '''Array "p" (below) represents the absolute Pr(genotype) where columns match for homo dom, hetero, and homo rec. \nRows represent the possible pairings of Parent 1 x Parent 2 \nValues k, m, n are not considered''')
    print(p)
    print(phd + '\n' + ph + '\n' + phr)
    print(pr_homo_dom + pr_hetero)


ros_iprb(19, 29, 25)
print('')

print(Color.Blue + '''ros_fib
Problem: Write an equation for fibonacci sequence in which 1 or more breeding pairs are produced.
Given: Positive integers n ≤ 40 and k ≤ 5.
Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair
and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs
(instead of only 1 pair).''' + Color.End)


def ros_fib(n, k):
    a, b = 0, 1
    for i in range(n):
        c = k * a
        a, b = b, c + b
    return a


print(ros_fib(23, 1))
print('')

print(Color.Blue + '''ros_gc
Given: Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.
Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated''' + Color.End)


def fasta_dict(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file

    sum_lines = file_len_v2(input_file)

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


def ros_unsortedgc(file):
    array = fasta_dict(file)
    gc_count = {}
    for i in array:
        y = ((array[i].count('G') + array[i].count('C')) / len(array[i])) * 100
        gc_count[i] = y
    return gc_count


def ros_gc(file):
    gc_count = ros_unsortedgc(file)
    y = max(gc_count, key=lambda i: gc_count[i])
    r = y.strip('>')
    print(r)
    print(gc_count[y])


print(fasta_dict('rosalind_gc.txt'))
print(ros_unsortedgc('rosalind_gc.txt'))
ros_gc('rosalind_gc.txt')
print('')

print(Color.Blue + '''ros_prot
Given: An RNA string s  corresponding to a strand of mRNA (of length at most 10 kbp).
Return: The protein string encoded by s .''' + Color.End)

# Where 'X' is a stop codon
rna_aa_table = {'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'I': ['AUU', 'AUC', 'AUA'],
                'M': ['AUG'],
                'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'],
                '*': ['UAA', 'UAG', 'UGA'], 'H': ['CAU', 'CAC'],
                'Q': ['CAA', 'CAG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'],
                'C': ['UGU', 'UGC'], 'W': ['UGG'],
                'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']}

print(rna_aa_table)


def ros_prot_all_rnaframes(Text):
    protein = []
    for i in range(len(Text) - 2):
        for key, val in rna_aa_table.items():
            if Text[i:i + 3] in val:
                protein.append(key)
    return protein


def ros_prot_rnaframe(Text, k):  # WHATS WRONG HERE????
    prot = ros_prot_all_rnaframes(Text)
    prot1 = prot[k::3]
    protein = ''.join(prot1)
    return protein


def ros_prot(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    zero = ros_prot_rnaframe(t, 0) + '\n'
    plus1 = ros_prot_rnaframe(t, 1) + '\n'
    plus2 = ros_prot_rnaframe(t, 2) + '\n'

    f = open(output_file, 'w')
    f.write('Frameshift: 0')
    f.write('\n')
    f.write(zero)
    f.write('\n')
    f.write('Frameshift: +1')
    f.write('\n')
    f.write(plus1)
    f.write('\n')
    f.write('Frameshift: +2')
    f.write('\n')
    f.write(plus2)
    f.write('\n')
    f.close()


ros_prot('rosalind_prot.txt')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_prot.txt') as g:
    print(g.read())
print('')

print(Color.Blue + '''ros_subs
Given: Two DNA strings s and t(each of length at most 1 kbp).
Return: All locations of t as a substring of s. ''' + Color.End)
print(
    Color.Red + '''Notes: This code returns positions AS FOR A 1-BASED NUMBERING SYSTEM. Patterncount uses a 0-based system.''' + Color.End)


def ros_subs(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            positions.append(i + 1)

    return positions


# Below modifies ros_subs to reflect the correct input for the Rosalind exercise:

def ros_subs_format(Pattern, Genome):
    x = ros_subs(Pattern, Genome)
    new_string = ''
    for i in x:
        new_string += str(i) + ' '
    print(new_string)


ros_subs_format('AAGGATCAA',
                'TGCAAGGATCAAAAGGATCTAAGGATCGCTAAGGATCAGTAAAGGATCAAGGATCGAAGGATCGAAGGATCCAAAGGATCTTGACAAGGATCCTAAGGATCAGAAGGATCCAAGGATCGTGCCTCAAAGGATCAAGGATCAAGGATCAAGAAGGATCAAAGGATCCTCAAGGATCAAAAGGATCAAGGATCAAGGATCTATTTGTAAGGATCGCTTTAAGGATCAAGGATCGTAAGGATCGAAATAGAAGGATCGCAAAGGATCTGTTTTTTGCAAGGATCTGTTAAGGATCGGATTGCACTGAAGGATCAAGGATCCTTAAAGGATCTGATGTTTGAAGGATCAAGGATCAAGGATCAAGGATCGAAAGGATCTGCAAAGGATCCAAGGATCAAAGGATCAAGGATCATCAAGGATCATAAGGATCAAAGGATCAAGGATCAAGGATCAAGGATCGAAAGGATCCATCAAGGATCATAAGGATCACACGATAACCAAGGATCAAAGGATCAGAAGGATCAAGGATCGTTTCTAAGGATCCAAAAGGATCATGCACTAAGGATCAAAGGATCAAAGGATCTCAGAAGGATCGTCTACGAAGGATCGGTCCACAAGGATCAAGGATCGAAGGATCGAAGGATCAAATCAAGGATCTCCAGATAAGGATCAAGGATCTCAGTAAACCAAGGATCAGTCCTTCGAAGGATCCAAAAGGATCAAAGGATCAAGGATCAAGGATCAACACAAGGATCAAGCTACGAGAATATCTTTGTGCAAGGATCACCTTAAGGATCGTAAAGGATCGGAAGGATCAAGGATCAAGGATCGAAGGATCTCGGCTAAGGATCAAGGATCGACAAGGATCATTACAAAGGATCCCGTACGGAAAAGGATCTGCAAGGATCGTGCCAAGGATCAAGGATCAAGGATCTCCCAAGGATCGGAAGGATC')
print('')

print(Color.Blue + '''ros_hamm
Given: Two DNA strings s and t(each of length at most 1 kbp).
Return: All locations of t as a substring of s. ''' + Color.End)


def ros_hamm(p, q):
    hamdist = 0
    for i in range(0, len(p)):
        if p[i] is not q[i]:
            hamdist += 1
    return hamdist


print('')

print(Color.Blue + '''ros_iev
Given: Given: Six positive integers, each of which does not exceed 20,000. The integers correspond to the number of
couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent
the number of couples having the following genotypes:

1.    AA-AA
2.    AA-Aa
3.    AA-aa
4.    Aa-Aa
5.    Aa-aa
6.    aa-aa

Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the
assumption that every couple has exactly two offspring. ''' + Color.End)


# input k,m,n,x,y,z,a,b,c refers to no. of matings for the pairings given by array h respectively
# input o refers to no. offspring born to each parent

def iev(k, m, n, x, y, z, a, b, c, o):
    pop = k + m + n + x + y + z + a + b + c
    pop_offspring = o * pop
    pk = k / pop
    pm = m / pop
    pn = n / pop
    px = x / pop
    py = y / pop
    pz = z / pop
    pa = a / pop
    pb = b / pop
    pc = c / pop

    p = np.array([[1, 0, 0], [.5, .5, 0], [0, 1, 0], [.5, .5, 0], [.25, .5, .25], [0, .5, .5], [0, 1, 0], [0, .5, .5],
                  [0, 0, 1]])
    h = np.array(
        [[pk], [pm], [pn], [px], [py], [pz], [pa], [pb], [pc]])

    pr_homo_dom = h[0][0] * p[0][0] + h[1][0] * p[1][0] + h[3][0] * p[3][0] + h[4][0] * p[4][0]
    pr_hetero = h[1][0] * p[1][1] + h[2][0] * p[2][1] + h[3][0] * p[3][1] + h[4][0] * p[4][1] + h[5][0] * p[5][1] + \
                h[6][0] * p[6][1] + h[7][0] * p[7][1]
    pr_homo_rec = h[4][0] * p[4][2] + h[5][0] * p[5][2] + h[7][0] * p[7][2] + h[8][0] * p[8][2]

    no_offspring_homdom = pr_homo_dom * pop_offspring
    no_offspring_hetero = pr_hetero * pop_offspring
    no_offspring_homrec = pr_homo_rec * pop_offspring

    phd = 'no. offspring(homozygous dominant)= %s' % no_offspring_homdom
    ph = 'no. offspring(heterozygous)= %s' % no_offspring_hetero
    phr = 'no. offspring(homozygous recessive)= %s' % no_offspring_homrec
    print(
        '''Array "h" (below) represents the ratio/fraction (Parent 1 x Parent 2) of each pairing in order: \nGG.GG, GG.Gg, GG.gg, Gg.GG, Gg.Gg, Gg.gg, gg.GG, gg.Gg, gg.gg \nwith respect to k, m, n values''')
    print(h)
    print(
        '''Array "p" (below) represents the absolute Pr(genotype) where columns match for homo dom, hetero, and homo rec. \nRows represent the possible pairings of Parent 1 x Parent 2 \nValues k, m, n are not considered''')
    print(p)
    print(phd + '\n' + ph + '\n' + phr)
    # print(no_offspring_homdom + no_offspring_hetero)


print(iev(18383, 16480, 18112, 0, 16410, 17402, 0, 0, 19456, 2))

print(Color.Blue + '''ros_fibd
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2 and
assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.
Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed
number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).
Given: Positive integers n≤100 and m≤20 .
Return: The total number of pairs of rabbits that will remain after the n -th month if all rabbits live for m months.''' + Color.End)

# I DON'T UNDERSTAND CODE BELOW!!! See http://stackoverflow.com/questions/17310051/fibonacci-rabbits-dying-after-arbitrary-of-months

# run for n months, rabbits die after m months.
# CODE START
# n, m = input("Enter months to run, and how many months rabbits live, separated by a space ").split()
# n, m = int(n), int(m)
# generations = [1, 1] #Seed the sequence with the 1 pair, then in their reproductive month.
# def fib(i, j):
#     count = 2
#     while (count < i):
#         if (count < j):
#             generations.append(generations[-2] + generations[-1]) #recurrence relation before rabbits start dying (simply fib seq Fn = Fn-2 + Fn-1)
#         elif (count == j or count == j+1):
#             print ("in base cases for newborns (1st+2nd gen. deaths)") #Base cases for subtracting rabbit deaths (1 death in first 2 death gens)
#             generations.append((generations[-2] + generations[-1]) - 1)#Fn = Fn-2 + Fn-1 - 1
#         else:
#             generations.append((generations[-2] + generations[-1]) - (generations[-(j+1)])) #Our recurrence relation here is Fn-2 + Fn-1 - Fn-(j+1)
#         count += 1
#     return (generations[-1])
#
# print (fib(n, m))
# print ("Here's how the total population looks by generation: \n" + str(generations))
# CODE END

print(Color.Blue + '''ros_mrna
Given: A protein string of length at most 1000 aa.
Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000.
(Don't neglect the importance of the stop codon in protein translation.) ''' + Color.End)


def ros_mrna_combinations(string):
    f = 1
    string_with_stop = string + 'X'
    for i in string_with_stop:
        if i in rna_aa_table:
            f *= len(rna_aa_table[i])
    return f


def ros_mrna(file, modulo):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    x = ros_mrna_combinations(t)
    rna_modulo = x % modulo
    return rna_modulo


print(ros_mrna('rosalind_mrna.txt', 1000000))
print('')

print(Color.Blue + '''ros_mrna
Given: A protein string P  of length at most 1000aa.
Return: The total weight of P . Consult the monoisotopic mass table.''' + Color.End)

monoisotypic_mass_table = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333,
}


def ros_prtm(file):
    input_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/' + file

    with open(input_file) as y:
        z = y.read()
        t = z.strip('\n')

    monoisotypic_mass = 0
    for i in t:
        if i in monoisotypic_mass_table:
            monoisotypic_mass += monoisotypic_mass_table[i]
    return round(monoisotypic_mass, 2)


print(ros_prtm('rosalind_prtm.txt'))
print('')

print(Color.Blue + '''ros_cons
Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.
Return: A consensus string and profile matrix for the collection. (If several possible consensus
strings exist, then you may return any one of them.)
NOTE: These functions have been adapted from Bioinformatics for Beginners to accommodate Rosalind
input/output formats. The code itself still functions in the same way''' + Color.End)

def len_list(file):
    x = fasta_dict(file)
    key_len = []
    for i in x:
        key_len.append(len(x[i]))
    return key_len


def initialise_dna_matrix(file):
    x = min(len_list(file))
    count = {}
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(x):
            count[symbol].append(0)
    return count


def count(file):
    x = fasta_dict(file)
    y = initialise_dna_matrix(file)
    k = len(y['G'])
    for i in x:
        for j in range(k):
            symbol = x[i][j]
            y[symbol][j] += 1  # for key[base] look at list index [j] in each row and add +1
    return k, y

def ros_cons(file):
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file
    x = count(file)
    min_len = x[0]
    count_di = x[1]
    consensus = ''
    for j in range(min_len):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count_di[symbol][j] > m:
                m = count_di[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    f = open(output_file, 'w')
    f.write(consensus)
    f.write('\n')
    for k in count_di:
        f.write('%s: ' % k)
        for v in count_di[k]:
            f.write('%s ' % v)
        f.write('\n')
    f.close()

print(count('rosalind_cons.txt'))
print('')
with open('D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-rosalind_cons.txt') as g:
    print(g.read())
print('')

print(Color.Blue + '''ros_orf
Given: A DNA string s  of length at most 1 kbp in FASTA format.
Return: Every distinct candidate protein string that can be translated from ORFs of s . Strings can be returned in any order''' + Color.End)

dna_aa_table = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'],
                'M': ['ATG'],
                'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'],
                'X': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
                'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
                'C': ['TGT', 'TGC'], 'W': ['TGG'],
                'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']}

def ros_prot_all_dnaframes(Text):
    protein = []
    for i in range(len(Text) - 2):
        for key, val in dna_aa_table.items():
            if Text[i:i + 3] in val:
                protein.append(key)
    return protein


def ros_prot_dnaframe(Text, k):  # WHATS WRONG HERE????
    prot = ros_prot_all_dnaframes(Text)
    prot1 = prot[k::3]
    protein = ''.join(prot1)
    return protein

def orf_string(string):
    M_index=[]
    X_index=[]
    for i in range(len(string)):
        if string[i]=='M':
            M_index.append(i)
        if string[i]=='X':
            X_index.append(i)

    orfs=[]
    for i in M_index:
        for j in X_index:
            if i<j:
                orfs.append(string[i:j])
            # else:
            #     orfs.append(string[i:])

    orfs_clean=[]
    for i in orfs:
        if ('X' not in i):
            orfs_clean.append(i)

    return orfs_clean

def ros_revc_string(string):
    rev_string = string[::-1]

    comp = ''
    for i in rev_string:
        if i == 'A':
            comp += 'T'
        elif i == 'T':
            comp += 'A'
        elif i == 'G':
            comp += 'C'
        elif i == 'C':
            comp += 'G'
        # else:
        #     comp += 'n'

    return comp

def protein_list(protein_array):
    protein_list=[]
    for x in protein_array:
        for y in protein_array[x]:
            klist=protein_array[x]
            for z in range(len(klist)):
                for h in range(len(klist[z])):
                    if klist[z][h] not in protein_list:
                        protein_list.append(klist[z][h])
    return protein_list

def ros_orf(file):
    fasta_map=fasta_dict(file)
    output_file = 'D:/MEGA/x- Jeff- Programs & Modules/Rosalind/Test files/output-' + file

    frame_array={}
    for i in fasta_map:
        frames=[]
        frames.append(ros_prot_dnaframe(fasta_map[i],0))
        frames.append(ros_prot_dnaframe(fasta_map[i], 1))
        frames.append(ros_prot_dnaframe(fasta_map[i], 2))
        frames.append(ros_prot_dnaframe(ros_revc_string(fasta_map[i]),0))
        frames.append(ros_prot_dnaframe(ros_revc_string(fasta_map[i]), 1))
        frames.append(ros_prot_dnaframe(ros_revc_string(fasta_map[i]), 2))
        frame_array[i]=frames

    protein_array={}
    for j in frame_array:
        protein_array[j]=[]
        for k in frame_array[j]:
            protein_array[j].append(orf_string(k))

    # print(protein_array)

    protein_frames=''
    for m in protein_array:
        protein_frames+=m + '\n'
        protein_frames+='pos strand, frame +0: ' + ', '.join(protein_array[m][0]) +'\n'
        protein_frames+='pos strand, frame +1: ' + ', '.join(protein_array[m][1]) +'\n'
        protein_frames+='pos strand, frame +2: ' + ', '.join(protein_array[m][2]) +'\n'
        protein_frames+='neg strand, frame +0: ' + ', '.join(protein_array[m][3]) +'\n'
        protein_frames+='neg strand, frame +1: ' + ', '.join(protein_array[m][4]) +'\n'
        protein_frames+='neg strand, frame +2: ' + ', '.join(protein_array[m][5]) +'\n'
        protein_frames+='\n'

    print(protein_frames)

    new_protein_array={}
    for a in protein_array:
        protein_list = []
        for b in protein_array[a]:
            klist=protein_array[a]
            for z in range(len(klist)):
                for c in range(len(klist[z])):
                    if klist[z][c] not in protein_list:
                        protein_list.append(klist[z][c])
                        new_protein_array[a]=protein_list

    # print(new_protein_array)

    all_protein_string=''
    for y in new_protein_array:
        all_protein_string+=y + ':\n'
        zlist=new_protein_array[y]
        all_protein_string+='\n'.join(zlist)
        all_protein_string += '\n'*2

    print(all_protein_string)

    f = open(output_file, 'w')
    f.write('All candidate proteins, defined as Met-STOP, for input sequences' +'\n'*2)
    f.write(protein_frames)
    f.write('Total candidate proteins without duplicates' +'\n'*2)
    f.write(all_protein_string)
    f.close()

ros_orf('rosalind_orf1.txt')
print('')


#cannot generate motif list of all motifs- takes FOREVER!
#this coding approach has been abandoned

# def all_motifs(string):
#     all_motifs=[]
#     for a in range(len(string)):
#         x=string[a:]
#         b=len(string)-a
#         for c in range(b):
#             y=string[a:-c]
#             if x !='' and x not in all_motifs:
#                 all_motifs.append(x)
#             if y !='' and y not in all_motifs:
#                 all_motifs.append(y)
#     return all_motifs
#
# def ros_lcsm(file):
#     string_map= fasta_dict(file)
#
#     keys=[]
#     for i in string_map:
#         keys.append(i)
#
#     motifs = all_motifs(string_map[keys[0]])
#     return motifs
#
#
# string_map= fasta_dict('rosalind_lcsm.txt')
# for i in string_map:
#     print(len(string_map[i]))

#New coding approach is to produce the check the longest possible strings against other sequences in every loop.

def long_str(string, string2):
    for a in range (len(string)):
        m=len(string)
        c = blist(a)
        print('\n main loop: ' + str(a))
        for b in c:
            k=b+1
            x = string[b:m - c[-k]]
            if x == string2:
                print(x)
                print(len(x))
                break

def blist(string):
    blist = []
    for b in range(string):
        blist.append((b + 1) - 1)
    return blist

def ros_lcsm(file):
    s_map = fasta_dict(file)

    s_map_values={}
    for i in s_map:
        s_map_values[i]=len(s_map[i])

    minkey= min(s_map_values, key=s_map_values.get)
    motif= s_map[minkey]
    entries=len(s_map)

    for a in range (len(motif)):
        m=len(motif)
        c = blist(a)
        for b in c:
            k=b+1
            x = motif[b:m - c[-k]]
            count=0
            for y in s_map:
                if x in s_map[y]:
                    submotif=x
                    count+=1
                    if count==entries:
                        print(submotif)
                        print(len(submotif))
                        break

#PROBLEM- its slow and doesnt break.
# ros_lcsm('rosalind_lcsm.txt')

    # for a in range (len(string)):
    #     m=len(string)
    #     c = blist(a)
    #     for b in c:
    #         k=b+1
    #         x = string[b:m - c[-k]]
    #         count=0
    #         for y in string_map:
    #             if x in string_map[y]:
    #                 substring=x
    #                 count+=1
    #                 if count==entries:
    #                     print(substring)
    #                     print(len(substring))
    #                     break

