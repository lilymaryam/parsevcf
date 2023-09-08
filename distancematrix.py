import argparse
import gzip
import bte 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree', required=True, type=str,help='path to MAT')
parser.add_argument('-s', '--samples', required=True, type=str,help='comma separated list of samples')

args = parser.parse_args()
tree = args.tree
samps = args.samples.split(',')


def dist_matrix(tree, samples):
    samp_ancs = {}
    samp_dist = {}
    t = bte.MATree(tree)
    #for each input sample, find path to root and branch lengths
    for s in samples:
        samp_ancs[s] = []
        samp_dist[s] = []
        for n in t.rsearch(s):
            samp_ancs[s].append(n.id)
            samp_dist[s].append(n.branch_length)

    #print(samp_ancs)
    #print(samp_dist)
    #matrix = []
    matrix = np.empty([len(samples),len(samples)])
    for i in range(len(samples)):
        s = samples[i]
        #matrix.append([])
    
        for j in range(len(samples)):
            os = samples[j]
            #add catch to prevent reiteration of already checked pairs 
            #don't check own array 
            if os == s:
                matrix[i][j] = '0'
            if os != s:
                #print('S',s)
                #print('os', os)

                #ancestors_s2 = [n.id for n in t.rsearch(os)]
                #print(ancestors, other_ancestors)

                #find lca, add up branch lengths? 
                s_path = 0 
                os_path = 0 
                for a in samp_ancs[s]:
                    #print('a', a)
                    #print("NODE",t.get_node(a))
                    #print(len(t.get_node(a).mutations))
                    #s_path += t.get_node(a).branch_length
                    s_path += len(t.get_node(a).mutations)
                    if a in samp_ancs[os]:
                        #print('s_path', s_path)
                        #print('lca?', a)
                        lca = a
                        #print(t.get_node(a))
                        break
                for a in samp_ancs[os]:
                    #os_path += t.get_node(a).branch_length
                    os_path += len(t.get_node(a).mutations)
                    if a == lca:
                        #print('os_path', os_path)
                        break
                #print(s,os)
                #print(s_path, os_path)
                matrix[i][j] = str(s_path + os_path)
    
    
    print(matrix)
    #print(f'sample\t'+'\t'.join(samples))
    #for s in range(len(samples)):
    #    string = '\t'.join(matrix[s])
    #    print(f'{samples[s]}\t{string}')
    return samples, matrix                    
    

samps, mat = dist_matrix(tree, samps)

#for m in range(len(mat)):
#    print(mat[m])

print(f'sample\t'+'\t'.join(samps))
for s in range(len(samps)):
    padding = ['-']*(s+1)
    print(len(samps))
    print(len(padding))
    print(len(padding)+len(mat[s]))
    string = '\t'.join(mat[s])
    print(f'{samps[s]}\t{string}')