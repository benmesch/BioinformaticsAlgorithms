#rewrite W2.py to work with strings, instead of integers

import numpy as np

def parse_graph_dict(graph_list):
    #turn an input adjacancy list into a dictionary that can be handled easier in python
    results = {}
    for i in graph_list:
        if " -> " in i:
            _pos = i.find(" -> ")
            _start = str(i[:_pos])
            _finish = i[_pos+4:]
            _finish_list = [str(x) for x in _finish.split(",")]
            results[_start] = _finish_list
    return results

def EulerianCycle(graph_dict, debug=False):
    '''form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
    while there are unexplored edges in Graph
        select a node newStart in Cycle with still unexplored edges
        form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking
        Cycle ← Cycle’
    return Cycle
    '''
    results = []
    if debug: print ("starting...",graph_dict)
    nodes = [x for x in graph_dict.keys()]
    current_node = str(np.random.choice(nodes,1)[0])
    #current_node = '00'
    #current_node = 6
    #current_node = 1
    #current_node = 1954
    if debug: print ("initial: ",current_node)
    current_cycle = []
    while len(graph_dict)>0:
        if debug: print ("current_node",current_node," , ",len(graph_dict)," , current cycle: ",current_cycle)
        if current_node in graph_dict:
            current_cycle.append(current_node)
            edges = graph_dict[current_node]
            if len(edges)==1:
                next_node = edges[0]
                graph_dict.pop(current_node)
                current_node = next_node
                if len(graph_dict)==0:
                    if debug: print ("   subresult: ",current_cycle)
                    results += current_cycle
            else:
                next_node = str(np.random.choice(edges,1)[0])
                graph_dict[current_node].remove(next_node)
                current_node = next_node
        else:
            results += current_cycle
            if debug: print ("   subresult: ",current_cycle)
            if debug: print ("   remaining: ",graph_dict)
            current_node = None
            for candidate in results:
                if candidate in graph_dict:
                    current_node = candidate
                    break
            current_cycle = []
            if current_node is not None:
                current_cycle = []
                #print (results.index(current_node)) #reset results to START with new starter...
                if debug: print ("original ",results)
                results = results[results.index(current_node):] + results[:results.index(current_node)]
                if debug: print ("new      ",results)
            else:
                break
    results.append(results[0])
    return [str(v) for v in results]


def EuclerianCycle_from_list(cycle_list, debug=False):
    #wrapper function for using a list input instead of a dictionary
    return EulerianCycle(parse_graph_dict(cycle_list), debug)

def FormatEuclerianCycle_from_list(cycle_list, debug=False):
    #wrapper function for using a list input (not a dict) and formatting for printing
    return "->".join(EulerianCycle(parse_graph_dict(cycle_list), debug))

'''
input = [
'0 -> 3'
,'1 -> 0'
,'2 -> 1,6'
,'3 -> 2'
,'4 -> 2'
,'5 -> 4'
,'6 -> 5,8'
,'7 -> 9'
,'8 -> 7'
,'9 -> 6']
print (EuclerianCycle_from_list(input, True))

input = []
with open('dataset_203_2.txt','r') as f: #dataset_203_2
    for line in f:
        read_line = line.rstrip('\n')
        input.append(read_line)
    f.close()
print (FormatEuclerianCycle_from_list(input))
'''




def node_degree(graph_dict):
    results = {}
    for node,edges in graph_dict.items():
        if node in results:
            results[node] -= len(edges) #out degree (negative)
        else:
            results[node] = -1 * len(edges)
        for out_connection in edges:
            if out_connection in results:
                results[out_connection] += 1
            else:
                results[out_connection] = 1
    return results

def find_unbalanced_nodes(node_degrees):
    results = {}
    for node,degree in node_degrees.items():
        if degree != 0:
            results[node] = degree
    return results

def EuclerianPath(graph_dict, debug=False):
    results = []
    node_degrees = node_degree(graph_dict)
    unbalanced_nodes = find_unbalanced_nodes(node_degrees)
    if len(unbalanced_nodes) != 2:
        print ("Graph must be nearly balanced to find a Euclerian Path")
        return None
    #if it makes it this far, the input graph is nearly balanced
    missing_edge = [None, None]
    for k,v in unbalanced_nodes.items():
        if v==-1:
            missing_edge[1]=k
        elif v==1:
            missing_edge[0]=k
    if missing_edge[0] is None or missing_edge[1] is None:
        print ("Graph is not nearly balanced")
        return None
    #add the missing edge to the graph...
    if missing_edge[0] in graph_dict:
        graph_dict[missing_edge[0]].append(missing_edge[1])
    else:
        graph_dict[missing_edge[0]] = [missing_edge[1]]
    if debug: print ("missing edge...",missing_edge)
    #graph_dict is now a fully balanced graph!
    results = EulerianCycle(graph_dict, debug)
    #re-organize cycle to end at the missing node...
    if debug: print ("original results...", results)
    for i in range(len(results)-1):
        if (str(results[i]) == str(missing_edge[0])) and (str(results[i+1]) == str(missing_edge[1])):
            results = results[i+1:] + results[1:i+1]
            #^doesnt apply to edge case where the missing edge is actually the final one in original results
            break
    return results

def EuclerianPath_from_list(input_list, debug=False):
    #wrapper function for using a list input instead of a dictionary
    return EuclerianPath(parse_graph_dict(input_list), debug)

def FormatEuclerianPath_from_list(input_list, debug=False):
    #wrapper function for using a list input (not a dict) and formatting for printing
    return "->".join(EuclerianPath(parse_graph_dict(input_list), debug))

'''
input = ['0 -> 2'
     ,'1 -> 3'
     ,'2 -> 1'
     ,'3 -> 0,4'
     ,'6 -> 3,7'
     ,'7 -> 8'
     ,'8 -> 9'
     ,'9 -> 6']
print (EuclerianPath_from_list(input))

input = []
with open('dataset_203_6.txt','r') as f: #dataset_203_2
    for line in f:
        read_line = line.rstrip('\n')
        input.append(read_line)
    f.close()
print (FormatEuclerianPath_from_list(input))
'''

#from bioinformatics II, w1...
def StringSpelledbyGenomePathProblem(kmers):
    k = len(kmers[0])
    result = kmers[0]
    for i in range(1,len(kmers)):
        if kmers[i][:-1] == kmers[i-1][1:]:
            result += (kmers[i][-1])
    return result

def DeBruijnGraphfromKmers(patterns):
    results = {}
    for kmer in patterns:
        prefix = kmer[:-1]
        if prefix in results:
            results[prefix].append(kmer[1:])
        else:
            results[prefix] = [kmer[1:]]
    return results

def EuclerianPath_from_patterns(k, pattern_list, debug=False):
    if debug: print (pattern_list)
    db_graph = DeBruijnGraphfromKmers(pattern_list)
    if debug: print (db_graph)
    #wrapper function for using a list input (not a dict) and formatting for printing
    if debug: print ("FINAL CYCLE...",results)
    return EuclerianPath(db_graph,debug)

def Reconstruct_String_from_patterns(k, pattern_list, debug=False):
    x = EuclerianPath_from_patterns(k, pattern_list, debug)
    return StringSpelledbyGenomePathProblem(x)

def FormatEuclerianPath_from_patterns(k, pattern_list, debug=False):
    x = EuclerianPath_from_patterns(k, pattern_list, debug)
    return "->".join(x)

'''
patterns = ['CTTA',
     'ACCA',
     'TACC',
     'GGCT',
     'GCTT',
     'TTAC']
#print (FormatEuclerianPath_from_patterns(4, patterns, True))

k=0
input = []
line_num = 0
lines_to_skip = 0
with open('dataset_203_7.txt','r') as f: #dataset_203_2
    for line in f:
        line_num += 1
        read_line = line.rstrip('\n')
        if line_num <= lines_to_skip:
            continue
        elif read_line=='Output:':
            break
        elif line_num == 1:
            #print (read_line)
            k = int(read_line)
        elif line_num >= lines_to_skip+1:
            input.append(read_line)
        else:
            break
    f.close()
print (Reconstruct_String_from_patterns(k, input))
'''
import itertools
def binary_kmer_space(k):
    #enumerate all possible binary strings of length k
    return ["".join(seq) for seq in itertools.product("01", repeat=k)]

def find_universal_circular_string(k):
    kmers = binary_kmer_space(k)
    db_graph = DeBruijnGraphfromKmers(kmers)
    cyc = EulerianCycle(db_graph)
    print (cyc)
    chop = -1 * (k-1)
    return StringSpelledbyGenomePathProblem(cyc)[:chop]

'''
k = 9
print (binary_kmer_space(k))
print (find_universal_circular_string(k))
'''

#kmers = ['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT','CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA']
#print (Reconstruct_String_from_patterns(4,kmers))
a=['ACC','ACT','ATA','ATT','CAC','CCG','CGA','CTG','CTG','GAA','GAT','GAT','TAC','TCT','TGA','TGA','TTC']
b=['ATA','ATT','TGA','TGA','GAT','TAC','ACT','AGC','TTC','CTT','CTG','CTG','GAT','AAG','GCT','TCT','GAA']

print (Reconstruct_String_from_patterns(3,a))
print (Reconstruct_String_from_patterns(3,b))
