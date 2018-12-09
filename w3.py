'''
penicillium fungus kills staphylococcus

antibiotic agentcs are compounds that kill bacteria

russian biologists gerogy gause and maria brazhnikova discovered the bacteria bacillus brevis,
which killed staphylococcus aureus and penicillin could be isolated from it

americans found a moldy cantaloupe that could grow penicillin colonies


peptides - mini proteins, short amino acid sequences
antibiotics are peptides (Richard Synge)

pharmaceutical companies have continued to develop antibiotics since WWII,
but pathogenic bacteria continue to evolve resistance to them

MRSA is methicillin-resistant staph aureus, and the leading cause of infections death in hospitals
MRSA death rate is higher than AIDS in the USA

***
    developing new antibiotics is a central challenge of modern medicine
    a difficult problem in this research is sequencing newly discovered antibioitics
        (determining the order of amino acids making up the antibiotic peptide)
***

Bacillus brevis bacteria creates many antibiotics,
including Tyrocidine B1 (10 amino acid long sequence below)
Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr
V   K   L   F   P   W   F   N   Q   Y
^ this is actually a cyclic peptide tho, so DNA searches should look for 10 diff strings starting at each position
^^ ALSO, this peptide isnt sythesized in ribosomes of Bacillus brevis! its non-ribosomal peptide (NRP)!
    tyrocidines and gramicidins are non-ribosomal peptides (NRPs) created by a giant protein called NRP Synthetase
    the giant enzyme pieces together antibiotic peptides without any RNA or genetic code
    NRP sythetase assembles peptides by growing them one amino acid at a time - the final product can form a cyclic final pattern after release from the NRP enzyme/factory
    NRPs often hae pharmaceutical applications due to molecular evolution of bacteria and fungi using thse NRPs
    NRPs can have antibacterial properties, or even anti-tumor or immunosuppressor properties! some are used for bacteria-bacteria communication ("quorum sensing")

    since NRPs are assembled without DNA, have to use mass spec to reconstruct and sequence them (not DNA sequencing to search for the parts of the DNA genome that sequence a peptide)
    also, NRPs are often circular / cyclic so its harder to do linear-based searches


"Central Dogma of Molecular Biology" - DNA makes RNA makes proteins
Gene is transcribed into RNA, which is translated into amino acid sequence of a protein
    transcription: replace T with U...
    translation: partition RNA into 3-mers called "codons"
    convert codons into one of 20 amino acids via the "genetic code"
    now you get an amino acid string
    64 RNA codons -> amino acids [3 codons are "stop"]

rna is ACGU, dna is ACGT


***
thousands of different DNA 3*x-mers can code an x-amino acid long peptide
in DNA, the sequence can start in 3 different places, and can translate either foward or backwards strand...
leading to 6 different reading frames per sing-stranded DNA string!

dont forget about the DNA string's reverse complement!!!

'''
def line_concat(result_list):
    #convert a list of integers into multiple lines
    for x in range(len(result_list)):
        print (result_list[x])

def string_into_codons(rna):
    if (len(rna)%3)!=0: return []
    return [rna[i:i+3] for i in range(0,len(rna),3)]

def create_codon_dict(filename="RNA_codon_table_1.txt", lines_to_skip = 0):
    results,inverse = {},{}
    line_num = 0
    with open(filename,'r') as f:
        for line in f:
            line_num += 1
            read_line = line.rstrip('\n')
            if line_num <= lines_to_skip:
                continue
            elif line_num >= lines_to_skip+1:
                codon = read_line.split(" ")
                if len(codon[1])!=1: continue
                results[codon[0]] = codon[1]
                if codon[1] in inverse:
                    inverse[codon[1]].append(codon[0])
                else:
                    inverse[codon[1]] = [codon[0]]
            else:
                break
        f.close()
    return results,inverse

def protein_translation(rna):
    print (rna)
    codon_dict,inverse_codon_dict = create_codon_dict()
    rna_codons = string_into_codons(rna)
    results = [codon_dict[i] for i in rna_codons]
    return ''.join(results)

'''a,b = create_codon_dict()
print (a)
print (b)

print ("")
for i in 'LEADER':
    print (len(b[i]))
'''
'''
print (protein_translation('CCCAGUACCGAGAUGAAU'))
print (protein_translation('CCUCGUACAGAAAUCAAC'))
print (protein_translation('CCUCGUACUGAUAUUAAU'))
print (protein_translation('CCCAGGACUGAGAUCAAU'))'''
'''rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
#print (protein_translation(rna))

input = []
line_num = 0
lines_to_skip = 0
with open('dataset_96_4.txt','r') as f:
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif line_num >= lines_to_skip+1:
            input = read_line
        else:
            break
    f.close()
print (protein_translation(input))'''


def peptide_to_rna(peptide):
    codon_dict,inverse_codon_dict = create_codon_dict()
    rna_enumerated = [inverse_codon_dict[i] for i in peptide]
    results = []
    for x in range(len(rna_enumerated)):
        if x==0:
            #first amino acid, just add directly to results, no concatenation
            for y in rna_enumerated[0]:
                results.append(y)
        else:
            subresults = []
            for y in results:
                for z in rna_enumerated[x]:
                    subresults.append(y+z)
            results = subresults
    return (results)

def rna_to_dna(rna_string):
    return rna_string.replace('U','T')

def dna_to_reverse_complement(dna_string):
    reverse = dna_string[::-1] #reverse the string
    complement = {'A':'T', 'C':'G', 'G': 'C', 'T':'A'}
    return ''.join([complement[x] for x in reverse])

def peptide_encoding_search(rna, peptide, debug = False):
    '''We say that a DNA string Pattern encodes an amino acid
    string Peptide if the RNA string transcribed from either
    Pattern or its reverse complement Pattern translates into
    Peptide.
    '''
    rna_list = peptide_to_rna(peptide)
    if debug: print (rna_list)
    target_dna = {rna_to_dna(x):'' for x in rna_list}
    target_complements = {dna_to_reverse_complement(x):'' for x in target_dna.keys()}
    for x in target_complements.keys():
        if x in target_dna:
            continue
        else:
            target_dna[x] = ''
    results = []
    peptide_length = len(peptide)*3
    search_length = len(rna) - peptide_length + 1
    for i in range(search_length):
        if rna[i:i+peptide_length] in target_dna:
            results.append(rna[i:i+peptide_length])
    return (results)
'''
input_rna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
input_peptide = 'MA'
#print (peptide_encoding_search(input_rna, input_peptide))

input = ''
peptide = ''
line_num = 0
lines_to_skip = 0
with open('dataset_96_7.txt','r') as f:
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif line_num == lines_to_skip+1:
            input = read_line
        elif line_num == lines_to_skip+2:
            peptide = read_line
        else:
            break
    f.close()
line_concat(peptide_encoding_search(input,peptide))
'''

def number_of_subpeptides_from_cyclic_length_n(n):
    return n * (n-1)

def create_integer_mass_dict(filename="integer_mass_table.txt", lines_to_skip = 0):
    results,inverse = {},{}
    line_num = 0
    with open(filename,'r') as f:
        for line in f:
            line_num += 1
            read_line = line.rstrip('\n')
            if line_num <= lines_to_skip:
                continue
            elif line_num >= lines_to_skip+1:
                codon = read_line.split(" ")
                if len(codon[0])!=1: continue
                results[codon[0]] = int(codon[1])
                if int(codon[1]) in inverse:
                    inverse[int(codon[1])].append(codon[0])
                else:
                    inverse[int(codon[1])] = [codon[0]]
            else:
                break
        f.close()
    return results,inverse

def prefix_masses(peptide):
    mass_dict,inverse_dict = create_integer_mass_dict()
    running_total = 0
    results = [0]
    for i in range(len(peptide)):
        running_total += mass_dict[peptide[i]]
        results.append(running_total)
    return results

import numpy as np

def linear_spectrum(peptide):
    #print (peptide)
    #return all possible sub masses in peptide
    prefix_mass_list = prefix_masses(peptide)
    results = [0]
    for i in range(len(peptide)+1):
        for j in range(i+1,len(peptide)+1,1):
            #print ("i,j",i," ",j, " -- ",prefix_mass_list[j], " minus ",prefix_mass_list[i]," ",peptide[i:j])
            results.append(prefix_mass_list[j]-prefix_mass_list[i])
    return np.sort(np.asarray(results, dtype=int)).tolist()

def cyclic_spectrum(peptide):
    #print (peptide)
    prefix_mass_list = prefix_masses(peptide)
    peptide_total_mass = prefix_mass_list[-1]
    results = [0]
    #blah = []
    for i in range(len(peptide)+1):
        for j in range(i+1,len(peptide)+1,1):
            #print ("i,j",i," ",j, " -- ",prefix_mass_list[j], " minus ",prefix_mass_list[i]," ",peptide[i:j])
            #blah.append(peptide[i:j])
            results.append(prefix_mass_list[j]-prefix_mass_list[i])
            if i > 0 and j < len(peptide):
                #add wrap-around sub peptides to account for the peptide being cyclic
                results.append(peptide_total_mass - (prefix_mass_list[j]-prefix_mass_list[i]))
                #blah.append(peptide[j:]+peptide[:i])
    return np.sort((np.asarray(results, dtype=int))).tolist()
    #print (blah)
    #return results
'''
a,b = create_integer_mass_dict()
print (a)
print (b)
'''
#print (prefix_masses('NQEL'))
#print (linear_spectrum('NQEL'))
#print (linear_spectrum('KWWHVWLIMDNANFFACTMRPSGLQFPKVWDTREAEYE'))

#print (linear_spectrum('VAQ'))
#print (cyclic_spectrum('DPDFARFVNIPKGG')) #dataset_98_4

#How many subpeptides does a linear peptide of given length n have? (Include the empty peptide and the entire peptide.)
'''
x = 12648 #answer is 79992277!
results = 1
while x > 0:
    results += x
    x -= 1
print (x)
print (results)
'''
def consistent_check(linear_spectrum, full_spectrum):
    intermediate = [x for x in full_spectrum]
    for i in linear_spectrum:
        if i not in intermediate:
            return False
        else:
            intermediate.remove(i)
    return True

def generate_pseudo_peptide(integer_spectrum):
    #turn a list of peptide weights into a string that can be input to the spectrum generators
    mass_dict,inverse_mass_dict = create_integer_mass_dict()
    dummy_peptide = ''
    for i in integer_spectrum:
        dummy_peptide += inverse_mass_dict[i][0]
    return dummy_peptide

def generate_pseudo_cyclo_spectrum(integer_spectrum):
    return cyclic_spectrum(generate_pseudo_peptide(integer_spectrum))

def generate_pseudo_linear_spectrum(integer_spectrum):
    #print ("   ",integer_spectrum,'->',dummy_peptide,'->>',x)
    return linear_spectrum(generate_pseudo_peptide(integer_spectrum))

def CyclopeptideSequencing(spectrum):
    peptide_total_mass = int(np.max(np.asarray(spectrum)))
    print ("spectrum, " , spectrum)
    print ("full weight, ",peptide_total_mass)
    mass_dict,inverse_mass_dict = create_integer_mass_dict()
    possible_masses = inverse_mass_dict.keys()
    candidate_peptides = []
    results = [np.asarray([x], dtype=int) for x in possible_masses] #hard code the process for the 1st expansion
    len_results = 0
    counter = 0
    #while counter < 4 and not (counter>0 and len_results==0): #len_results>0 or counter==0: #
    while not (counter>0 and len_results==0):
        print ("counter #",counter," - ",len_results," starting candidates")#," -- ",results)
        counter += 1
        if counter > 1: #only expand the candidate list AFTER the first iteration
            new_results = []
            for x in results:
                for y in possible_masses:
                    new_results.append(np.append(x,y))
            results = [x for x in new_results]
        remove_candidates = []
        for i in range(len(results)):
            candidate = results[i]
            if np.sum(candidate) == peptide_total_mass:
                #print ("    candidate w same mass! ", candidate)
                #print ("    *",generate_pseudo_cyclo_spectrum(candidate))
                if generate_pseudo_cyclo_spectrum(candidate) == spectrum:
                    candidate_peptides.append(candidate)
                remove_candidates.append(i)
            else:
                dummy_linear_spectrum = generate_pseudo_linear_spectrum(candidate)
                if not consistent_check(dummy_linear_spectrum, spectrum):
                    remove_candidates.append(i)
        #print ("remove candidates ",remove_candidates)
        results = [results[x] for x in range(len(results)) if x not in remove_candidates]
        len_results = len(results)
        #print (" ending results... ",len_results)
    print ("counter #",counter," - ",len_results," candidates left. end search")
    return ['-'.join([str(y) for y in x.tolist()]) for x in candidate_peptides]
    '''
        Peptides ← a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ← Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    **if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                **else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides
    '''
#print (CyclopeptideSequencing([0, 113, 128, 427]))
#print (CyclopeptideSequencing([0, 113, 128, 186, 241, 299, 314, 427]))
#print ('186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186') #correct answer

#input = [0,71,87,99,113,114,114,115,131,137,147,158,208,212,213,218,227,245
    #,252,261,262,289,295,323,326,326,332,344,374,376,399,403,410,426,431,440
    #,457,470,473,489,502,513,540,541,544,557,571,584,587,588,615,626,639,655
    #,658,671,688,697,702,718,725,729,752,754,784,796,802,802,805,833,839,866
    #,867,876,883,901,910,915,916,920,970,981,991,997,1013,1014,1014,1015,1029
    #,1041,1057,1128]
input = [0,87,99,103,113,128,128,129,129,186,212,215,216,216,232,257,257,285
    ,314,315,344,344,344,345,360,398,413,443,444,447,472,473,473,501,526,530
    ,542,560,572,576,601,629,629,630,655,658,659,689,704,742,757,758,758,758
    ,787,788,817,845,845,870,886,886,887,890,916,973,973,974,974,989,999,1003
    ,1015,1102]
line_concat(CyclopeptideSequencing(input)) #dataset_100_6.txt
