import argparse
import gzip
import bte 
import numpy as np
import ete3

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree', required=True, type=str,help='path to MAT')
parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
#parser.add_argument('-s', '--samples', required=True, type=str,help='comma separated list of samples')


args = parser.parse_args()
tree = args.tree
t = ete3.Tree(tree, format=1)
if args.samples:
    samps = args.samples.split(',')
else:
    samps = sorted([leaf.name for leaf in t])

    


def path_to_root(ete_tree, node):
    # Browse the tree from a specific leaf to the root
    #print(node)
    node = ete_tree.search_nodes(name=node)[0]
    path = [node]
    #print(node)
    while node:
        node = node.up
        path.append(node)
    #print(path)
    return path


def dist_matrix(tree, samples):
    samp_ancs = {}
    samp_dist = {}
    #t = bte.MATree(tree)
    #t = ete3.Tree(tree, format=1)
    #for each input sample, find path to root and branch lengths
    for s in samples:
        s_ancs = path_to_root(tree, s)
        samp_ancs[s] = s_ancs
    #create matrix for samples
    matrix = np.empty([len(samples),len(samples)])
    for i in range(len(samples)):
        s = samples[i]

        for j in range(len(samples)):
            os = samples[j]
            #Future goal: add catch to prevent reiteration of already checked pairs 
            if os == s:
                matrix[i][j] = '0'
            if os != s:
                
                #find lca, add up branch lengths
                s_path = 0 
                os_path = 0 
                for a in samp_ancs[s]:
                    
                    s_path += a.dist
                    if a in samp_ancs[os]:
                        
                        lca = a
                        s_path -= a.dist
                        #print(t.get_node(a))
                        break
                for a in samp_ancs[os]:
                    
                    os_path += a.dist
                    if a == lca:
                        #print('os_path', os_path)
                        os_path -= a.dist
                        break
                #print(s,os)
                #print(s_path, os_path)
                matrix[i][j] = int(s_path + os_path)
    #print(matrix)
    return samples, matrix                    
    

samps, mat = dist_matrix(t, samps)

#for m in range(len(mat)):
#    print(mat[m])


print(f'sample\t'+'\t'.join(samps))
for i in range(len(samps)):
    #strng = np.array2string(mat[i], separator='\t')[1:-1]
    line = [ str(int(count)) for count in mat[i]]
    print(f'{samps[i]}\t' + '\t'.join(line))
    