'''
Lecture 4 notes... [code below]
    in reality, mass spec does not give perfect results
        mass spectrometers generate "noisy" spectra that are far from ideal — they are characterized by having both false masses and missing masses
        A false mass is present in the experimental spectrum but absent from the theoretical spectrum;
        a missing mass is present in the theoretical spectrum but absent from the experimental spectrum
    suprising new masses will appear or not appear unexpectedly
    so we have to SCORE spectra to grade how similar they are to one another

    score(a,b) = number of integers in common in list a and b

    "cut" - minimum score needed to not eliminate a candidate Spectrum

    ***another branch and bound algorithm....

    leaderboard cyclo peptide sequencing:
        at each iteration, still expand every candidate 18 ways (same)
        cut low-scoring peptides (keep top-N, with ties)
        update current_leader if there is a higher scoring candidate (with mass=parent mass)
        eliminate all candidates with mass > parent mass
        repeat til leaderboard is empty

    warning:
    this is a heuristic method, sacrificing precision - may miss the ultimate solution by
    eliminating it early on (if the optimal solution starts out badly)

    ^works tho, for 90% match... but not for 75% experimental match

    another wrench:
    NRPs contain more non standard amino acids, because they are free from the Central Dogma
        example: Ornithine is a non-standard amino acid (not used by ribosomal proteins)
    so have to assume ANY integer between 57-200 can act as the mass of an amino acid!!!
        instead of an alphabet of 18 weights, we run algo with 144 potential masses.
    ^this approach may be WORSE off for 90% experimental matches...
        highest scoring candidate may incorrectly use non standard amino acids :(


    spectral convolution!
        goal: reduce number of potential amino acids to consider.
        better approach: subtracting the different integers in the spectrum,
            will give potential single-amino-acid weights
            example: spectrum has 484 and 355, but not 129. but 484-355=129!
                so the 129 IS in the mass spec results, just not explicitly...
                so 129 should be in the algo for this spec as a potential candidate
            the M-most frequent differences here are likely amino acids' weights that should be considered
            ^by narrowing the alphabet of potential amino acids to only M possibilities,
            this approach can sequence peptides even at just 75% experimental accuracy!
            ^tho real world mass spec is way less accurate than 75% [33% possible... need better algo adjustments then!]
'''
import numpy as np

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

def expand_spectrum(peptide, cyclic = True):
    #print (peptide)
    #return all possible sub masses in peptide
    prefix_mass_list = prefix_masses(peptide)
    peptide_total_mass = prefix_mass_list[-1]
    results = [0]
    for i in range(len(peptide)+1):
        for j in range(i+1,len(peptide)+1,1):
            #print ("i,j",i," ",j, " -- ",prefix_mass_list[j], " minus ",prefix_mass_list[i]," ",peptide[i:j])
            results.append(prefix_mass_list[j]-prefix_mass_list[i])
            if cyclic and i > 0 and j < len(peptide):
                #add wrap-around sub peptides to account for the peptide being cyclic
                results.append(peptide_total_mass - (prefix_mass_list[j]-prefix_mass_list[i]))
    return np.sort((np.asarray(results, dtype=int))).tolist()

def spectrum_list_to_dict(spec_list):
    results = {}
    for i in spec_list:
        if i in results:
            results[i] += 1
        else:
            results[i] = 1
    return results

def PeptideScoring(peptide, experimental_spectrum, cyclic = True):
    theoretical_spectrum = expand_spectrum(peptide, cyclic)
    theoretical_spectrum_dict = spectrum_list_to_dict(theoretical_spectrum)
    experimental_spectrum_dict = spectrum_list_to_dict(experimental_spectrum)
    result = 0
    for i in theoretical_spectrum_dict.keys():
        if i in experimental_spectrum_dict:
            #i is in both spectra, add the minimum
            theo_score = theoretical_spectrum_dict[i]
            expr_score = experimental_spectrum_dict[i]
            if theo_score > expr_score:
                result += expr_score
            else:
                result += theo_score
    return result
print (PeptideScoring('PEEP',[0,97,97,97,100,129,194,226,226,226,258,323,323,355,393,452],False)) #toy example
'''
input = "NQEL"
#print (input)
#print (cyclic_spectrum(input))

spectra = [0,99,113,114,128,227,257,299,355,356,370,371,484]
print (PeptideScoring(input,spectra, cyclic = True))
print (PeptideScoring(input,spectra, cyclic = False))

input = ''
spectra = ''
line_num = 0
lines_to_skip = 0
with open('dataset_102_3.txt','r') as f:
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif line_num == lines_to_skip+1:
            input = read_line
        elif line_num == lines_to_skip+2:
            spectra = [int(x) for x in read_line.split(' ')]
        else:
            break
    f.close()
print (PeptideScoring(input,spectra, cyclic = True))
'''

def generate_pseudo_peptide(integer_spectrum):
    #turn a list of peptide weights into a string that can be input to the spectrum generators
    mass_dict,inverse_mass_dict = create_integer_mass_dict()
    dummy_peptide = ''
    for i in integer_spectrum:
        dummy_peptide += inverse_mass_dict[i][0]
    return dummy_peptide

def trim_leaderboard(leaderboard, spectrum, N):
    if len(leaderboard)==0: return leaderboard
    print ("  ","start trim with..." , leaderboard.shape[0])
    scores = np.asarray([],dtype=int)
    for i in leaderboard:
        s = PeptideScoring(generate_pseudo_peptide(i.tolist()), spectrum, False)
        #if s>0:
        scores = np.append(scores,s)
    if len(scores)<N: return leaderboard
    cutoff = scores[np.argsort(scores)[-N]] #nth highest score #np.sort(scores*-1)[N]*-1
    results = leaderboard[scores>=cutoff]
    print ("  ","cutoff: ",cutoff,", min: ",np.min(scores),", max: ",np.max(scores)," ... after trim: ",results.shape[0])
    '''
    if leaderboard.shape[0]==324:
        print (leaderboard.tolist())
        print (scores)
        print (np.sort(scores*-1))
        print ("better cutoff???",np.sort(scores*-1)[:N]*-1)
    '''
    return results

def LeaderboardCyclopeptideSequencing(spectrum, N):
    peptide_total_mass = int(np.max(np.asarray(spectrum)))
    print ("~~~~~~~~~~spectrum, " , spectrum)
    print ("~~~~~~~~~~full weight, ",peptide_total_mass)
    mass_dict,inverse_mass_dict = create_integer_mass_dict()
    possible_masses = inverse_mass_dict.keys()
    leaderboard = np.asarray([np.asarray([x], dtype=int) for x in possible_masses])
    leader_peptides = [[]]
    leader_score = 0
    counter = 0
    print_more = ""
    while not (counter>0 and leaderboard.size==0):
        expected_size = leaderboard.shape[0]*len(possible_masses) + 0
        if (counter>0): print_more = " -> " + str(expected_size) + " expected new candidates"
        print ("counter #",counter," - ",leaderboard.shape[0]," starting candidates",print_more)#," -- ",results)
        counter += 1
        if counter > 1: #only expand the candidate list AFTER the first iteration
            expanded_board = []
            for x in leaderboard:
                for y in possible_masses:
                    new = np.append(x,y)
                    new_mass = np.sum(new)
                    #remove any candidates from the leaderboard that already have grown to weigh more than the entire peptide
                    if new_mass < peptide_total_mass or new_mass == peptide_total_mass:
                        expanded_board.append(new)
            leaderboard = np.asarray([x for x in expanded_board])
            if expected_size != leaderboard.shape[0] and leaderboard.shape[0] > 0:
                print ("   removed ",expected_size-leaderboard.shape[0]," too-large candidates. new len = ",leaderboard.shape[0])
        #leaderboard = np.asarray([x for x in leaderboard if (np.sum(x) < (peptide_total_mass-50) or np.sum(x) == peptide_total_mass)])
        for candidate in leaderboard:
            if np.sum(candidate) == peptide_total_mass: #if candidate mass is the same as the target peptide...
                candidate_score = PeptideScoring(generate_pseudo_peptide(candidate.tolist()), spectrum, True) #check with cyclo score, tho we will trim with linear score
                if candidate_score > leader_score:
                    print ("     new leader! (score ",candidate_score,")")
                    leader_peptides = [candidate]
                    leader_score = candidate_score
                elif candidate_score == leader_score and candidate_score > 0:
                    print ("     candidate ties current leader (",len(leader_peptides), " co-leaders)")
                    leader_peptides.append(candidate)
        leaderboard = trim_leaderboard(leaderboard, spectrum, N)
        if leader_score>0: print ("  *current leader score: ",leader_score,", ",leaderboard.shape[0]," candidates left")
    #print ("counter #",counter," - ",leaderboard.size," candidates left. end search. best score: ",leader_score)
    #return ['-'.join([str(y) for y in x.tolist()]) for x in leaderboard]
    #print (leader_peptides)
    return [x.tolist() for x in leader_peptides]

'''
    LeaderboardCyclopeptideSequencing(Spectrum, N)
        Leaderboard ← set containing only the empty peptide
        LeaderPeptide ← empty peptide
        while Leaderboard is non-empty
            Leaderboard ← Expand(Leaderboard)
            for each Peptide in Leaderboard
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                        LeaderPeptide ← Peptide
                else if Mass(Peptide) > ParentMass(Spectrum)
                    remove Peptide from Leaderboard
            Leaderboard ← Trim(Leaderboard, Spectrum, N)
        output LeaderPeptide
'''
#spectrum = [0,71,113,129,147,200,218,260,313,331,347,389,460]
#print (LeaderboardCyclopeptideSequencing(spectrum, 10))
#print ("real answer: 113-147-71-129")


spectrum = []
N = 0
line_num = 0
lines_to_skip = 0
with open('dataset_102_8.txt','r') as f:
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif line_num == lines_to_skip+1:
            N = int(read_line)
        elif line_num == lines_to_skip+2:
            spectrum = [int(x) for x in read_line.split(' ')]
        else:
            break
    f.close()

#print ('-'.join([str(x) for x in LeaderboardCyclopeptideSequencing(spectrum, N)]))

spectrum = []
N = 0
line_num = 0
lines_to_skip = 0
with open('dataset_102_10.txt','r') as f:
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif line_num == lines_to_skip+1:
            N = int(read_line)
        elif line_num == lines_to_skip+2:
            spectrum = [int(x) for x in read_line.split(' ')]
        else:
            break
    f.close()


#for x in LeaderboardCyclopeptideSequencing(spectrum, N):
#    print ('-'.join([str(y) for y in x]))


#sepctrum = [0,97,99,114,128,147,147,163,186,227,241,242,244,260,261,262,283,291,333,340,357,385,389,390,390,405,430,430,447,485,487,503,504,518,543,544,552,575,577,584,632,650,651,671,672,690,691,738,745,747,770,778,779,804,818,819,820,835,837,875,892,917,932,932,933,934,965,982,989,1030,1039,1060,1061,1062,1078,1080,1081,1095,1136,1159,1175,1175,1194,1194,1208,1209,1223,1225,1322]
#N = 1000
#for x in LeaderboardCyclopeptideSequencing(spectrum, N):
#    print ('-'.join([str(y) for y in x]))



def convolution(spectrum):
    spectrum = np.sort(np.asarray(spectrum)).tolist()
    results = []
    for i in range(len(spectrum)):
        for k in spectrum[i+1:]:
            diff = spectrum[i]-k
            if diff<0: diff *= -1
            if diff>0: results.append(diff)
    return results


#print (' '.join([str(x) for x in convolution([0, 137, 186, 323])]))

spectrum = [841,97,0,630,114,212,744,743,128,536,211,468,487,792,581,97,342,827,419,613,128,615,113,353,602,291,699,256,439,516,842,227,728,228,840,727,858,374,664,325,163,827,567,225,260,631,340,730,516,324,453,439,955,695,502,388,858,115]
print (spectrum)
print (' '.join([str(x) for x in convolution(spectrum)]))

spectrum = [0,57,118,179,236,240,301]
print (spectrum)
print (np.sort(np.asarray(convolution(spectrum)).tolist()))
