import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--VCF', required=True, type=str,help='path to VCF to be processed')
parser.add_argument('-d', '--working_directory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-tbmf', '--tb_maskfile', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-cf', '--bed_coverage_file', required=False, type=str, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")
parser.add_argument('-cd', '--coverage_depth', required=False, default=10, type=int, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")

#NOTICE: Do not use this script if you have space limitations
args = parser.parse_args()
vcf = args.VCF
wd = args.working_directory
tbmf = args.tb_maskfile
cf = args.bed_coverage_file
cd = args.coverage_depth

print('vcf', vcf)
print('wd', wd)
print('tbmf', tbmf)
print('cf', cf)
print('cd', cd)

#makes sure input path wont cause error
if wd[-1] != '/':
    wd = wd+'/'

#Functions

#for lines where len(ref)==len(alt), look for snps instead of processing as one large chunk                        
def find_snps(line):
    ref = line[3]
    alt = line[4]
    end = len(ref)
    snps = []
    for i in range(end):
        if ref[i] != alt[i]:
            snps.append(i)
    #generate new lines for diff file
    lines = []
    for s in range(len(snps)):
        lines.append([alt[snps[s]], str(int(line[1])+snps[s]), '1'])

    return lines

#for lines where len(ref) > 1 AND len(alt) == 1, these deletions are easily processed 
#currently deletions and missing data are all converted to '-'
def process_dels(line):
    ref = line[3]
    alt = line[4]
    end = len(ref)
    assert len(alt) == 1

    #make sure remaining alt nucleotide is the same as the corresponding ref nuc
    if alt[0] == ref[0]:
        l = ['-', str(int(line[1])+1), str(end-1)]

    elif alt[0]  == '-':
        #this is for missing data
        l = ['-', line[1], str(end)]

    else:
        #this is a scenario that could be represented by a snp at the first ref position
        l = ['-', line[1], str(len(line[3]))]
        
    return l


#for sitations where reference and alt do not align
#will mask entire reference 
def process_others(line):
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    return l

#read TB mask file (different than coverage file) and generate sites to be masked
def mask_TB(tbmf):
    tb_sites = {}
    #count = 0 
    with open(tbmf) as file:
        for line in file:
            #count += 1
            line=line.strip().split()
            tb_sites[int(line[1])] = int(line[2])
    #print('regions', tb_sites)
    return tb_sites

#read coverage file and generate sites to be masked 
#if coverage does not have HR37c reference it will throw an error (this can be changed)
#old version

def mask_low_depth(cf, cd):
    ld_sites = {}
    #count = 0 
    prev = None
    with open(cf) as cf:
        for line in cf:
            #count += 1
            #for currect bed coverage file 
            if line.startswith('NC_000962.3'):
                line = line.strip().split()
                if int(line[3]) < cd:
                    #print('prev', prev)
                    #print(line)
                    if prev == None:
                        ld_sites[int(line[1])] = int(line[2])
                        prev = [int(line[1]), int(line[2])]
                    else:
                        if int(line[1]) == prev[1]:
                            ld_sites[prev[0]] = int(line[2])
                            print('squish', 'prev', prev, 'line', line)
                            prev[1] = int(line[2])

                        else:
                            print('no squish')
                            ld_sites[int(line[1])] = int(line[2])
                            prev = [int(line[1]), int(line[2])]
                    

    if ld_sites == {}:
        raise Exception('coverage file has incorrect reference')
    #print('regions', count)
    #print('ld_sites', len(ld_sites))
    #for l in ld_sites:
    #    print(l, ld_sites[l])
    return ld_sites

def condense_mask_regions(cf,cd,tbmf):
    ld_sites = mask_low_depth(cf, cd)
    tb_sites = mask_TB(tbmf)
    tb_keys = sorted(tb_sites.keys())
    ld_keys = sorted(ld_sites.keys())
    all_sites = {}
    #print('ld',ld_keys)
    #print('tb',tb_keys)
    tb_keys_ind = 0
    ld_keys_ind = 0
    cont = 0
    while tb_keys_ind < len(tb_keys) or ld_keys_ind < len(ld_keys):
        if len(all_sites) > 1:
            all_sites_keys = sorted(all_sites.keys())
            print(all_sites_keys[-1])
        if tb_keys_ind < len(tb_keys) and ld_keys_ind < len(ld_keys): 
            tb_start = tb_keys[tb_keys_ind]
            tb_end =  tb_sites[tb_keys[tb_keys_ind]]
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]

            if ld_start < tb_start:
                if ld_end < tb_start:
                    print(f'ld{ld_keys_ind} is below tb{tb_keys_ind}')
                    print('ld', ld_start, ld_end)
                    print('tb', tb_start, tb_end)
                    all_sites[ld_start] = ld_end
                    ld_keys_ind += 1
                elif ld_end >= tb_start:
                    print('left overlap')
                    print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                    all_sites[ld_start] = tb_end
                    tb_keys_ind += 1
                    ld_keys_ind += 1
            
            elif ld_start <= tb_end and ld_end > tb_end:
                print('right over lap')
                #print('left overlap')
                print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                all_sites[tb_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            elif ld_start > tb_end:
                print('no overlap')
                print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                print('ld', ld_start, ld_end)
                print('tb', tb_start, tb_end)
                all_sites[tb_start] = tb_end
                tb_keys_ind += 1

            elif ld_start >= tb_start and ld_end <= tb_end:
                print('full overlap: ld inside')
                print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                print('ld', ld_start, ld_end)
                print('tb', tb_start, tb_end)
                all_sites[tb_start] = tb_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            
            elif ld_start <= tb_start and ld_end >= tb_end:
                print('full overlap: tb inside')
                print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                all_sites[ld_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1
                print('ld', ld_start, ld_end)
                print('tb', tb_start, tb_end)

            else:
                print('other', 'what else could happen?')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                #all_sites[tb_start]


                #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)
            #print('tb',tb_keys[tb_keys_ind], tb_sites[tb_keys[tb_keys_ind]])
            #print('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
        elif tb_keys_ind >= len(tb_keys) and ld_keys_ind < len(ld_keys):
            print('no more tb masks, ld only')
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]
            all_sites[ld_start] = ld_end
            #tb_keys_ind += 1
            ld_keys_ind += 1


        elif ld_keys_ind >= len(ld_keys) and tb_keys_ind < len(tb_keys):
            print('no more ld masks, tb only')
            tb_start = tb_keys[tb_keys_ind]
            tb_end = tb_sites[tb_keys[tb_keys_ind]]
            all_sites[tb_start] = tb_end
            #tb_keys_ind += 1
            tb_keys_ind += 1 

        cont += 1
        #print('tb ind', tb_keys_ind)
        #print('ld ind', ld_keys_ind)
        #print('len tb', len(tb_keys))
        #print('len ld', len(ld_keys))
        #if cont == 10000:
        #    print(all_sites)
        #    break
    return all_sites

    
'''
def mask_low_depth(cf, cd):
    ld_sites = {}
    #count = 0 
    #prev = None
    start = None
    end = None 
    with open(cf) as cf:
        for line in cf:
            #count += 1
            #for currect bed coverage file 
            if line.startswith('NC_000962.3'):
                line = line.strip().split()
                if start == None:
                    if int(line[3]) < cd:
                        print('FIRST')
                        start = int(line[1])
                        end = int(line[2])
                        print(start,end)
                        ld_sites[start] = end
                    #prev = line
                elif int(line[1]) == end and int(line[3]) <cd:
                    print('COMPRESS?')
                    print('start', start)
                    print('end', end)
                    print('line', line)
                    #ld_sites[start] = end
                    #print('prev', prev)
                    #print('line',line)
                    #print('low coverage')
                    #print('compres?','prev', prev,'line',line)
                if int(line[1]) == end and int(line[3]) >= cd:
                    start = int(line[1])
                    end = int(line[2])
                    print('dont compress')
                    print('start', line[1])
                    print('end', line[2])
                    print('line', line)
                    #ld_sites[start] = end
                else:
                    pass
                    #print('cant compress')
                    #print('prev', prev)
                    #print('line',line)
                    #print('correct covereage?','prev', prev,'line',line)
                #prev = line
                
    print(ld_sites)
    if ld_sites == {}:
        raise Exception('coverage file has incorrect reference')
    #for ld in ld_sites:

    return ld_sites'''



#run this function on all generated diff functions to mask full regions
#note that these regions are not considered missing data as they are universally masked
'''
def mask_TBsites(diff_file, tbmf, sample):
    #print('TB mask regions')
    tb_sites =  mask_TB(tbmf)
    #keys = sorted(tb_sites.keys())
    print('tb sites',len(tb_sites))
    #for i in tb_sites:
    #    print(i)
    #print('MASK TB SITES', diff_file)
    #print(sample)
    with open(diff_file) as df:
        with open(f'{wd}{sample}.masked.diff', 'w') as out:
            for line in df:
                keys = sorted(tb_sites.keys())
                if not line.startswith('>'):
                    line = line.strip().split()
                    print('position', line[1])
                    prevkey = None
                    for k in keys:
                        
                        
                        #print('k', k, 'tb_sites', tb_sites[k])
                        #check if greater than first value in coverage range
                        if int(line[2]) == 1:
                            print('SNP')
                            if int(line[1]) >= k:
                                #print('k', line[1], k)
                                #check if less than final value in coverage range
                                #indicates that the position should be masked
                                #break after to save time
                                #print()
                                
                                if int(line[1]) < tb_sites[k]:
                                    #print(f'in a range: MASK {line[1]}')
                                    #line[4] = '-'
                                    #line[-1] = '1'
                                    #should i do this?
                                    #missing += 1
                                    #print('start', k)
                                    #print('pos', line[1])
                                    #print('end', ld_sites[k])
                                    #print
                                    break

                                
                                #if greater than first and last value, remove range from dict
                                #time saver
                                else:
                                    #print('delete?')
                                    #print('start', k)
                                    #print('pos', line[1])
                                    #print('end', tb_sites[k])
                                    distance = int(tb_sites[k])-k
                                    #print('dists', distance)
                                    out.write('\t'.join(['-',str(k), str(distance)])+'\n')
                                    
                                    dist = tb_sites[k]-k
                                    #prevkey = [k, dist]
                                    #print('update','k', k, 'prevkey', prevkey)
                                    tb_sites.pop(k)
                                    #prevkey = (k, )
                                    #print('len dict', len(tb_sites))
                                    
                                    #need to set up value for masks
                            
                            elif int(line[1]) < k:
                                #print('range', k, tb_sites[k])
                                #print('no mask?', line[1])
                                #print('\t'.join(line))
                                
                                out.write('\t'.join(line)+'\n')
                            
                                #alts = line[4].split(',')
                                #alt = alts[int(var)-1]
                                #line[4] = alt
                                #line[-1] = '1'
                                
                                break
                        
                        else:
                            print('init','k', k, 'prevkey', prevkey)

                            print('not snp!!!!!', line)
                            #print('k', line[1], k)


                            #the indel positions are all less than the mask region
                            if int(line[1])+int(line[2]) < k and int(line[1]) < k:
                                print('no overlap, too soon')
                                #print('range', k, tb_sites[k])
                                #print('no mask?', line[1])
                                #print('\t'.join(line))
                                
                                out.write('\t'.join(line)+'\n')
                                break

                            #if greater than first and last value, remove range from dict
                            elif int(line[1]) > tb_sites[k]:
                                print('no overlap, too late')
                                print('non-SNP remove from lib')
                                print('delete?')
                                print('start', k)
                                print('pos', line[1])
                                print('end', tb_sites[k])
                                distance = int(tb_sites[k])-k
                                #print('dists', distance)
                                out.write('\t'.join(['-',str(k), str(distance)])+'\n')
                                #dist = tb_sites[k]-k
                                #prevkey = [k,dist]
                                
                                #print('update', prevkey )
                                tb_sites.pop(k)
                                #print('len dict', len(tb_sites))
                            

                            #the end of the indel overlaps with the beginning of the mask region
                            elif int(line[1])+int(line[2]) >= k and int(line[1]) < k:
                                print('OVERLAP!!!!! end of snp and beginnin of mask', k, line, tb_sites[k])
                                print('end of snp', int(line[1])+int(line[2]))
                                print('start of mask', k)
                                print( 'start of snp', int(line[1]))
                                print('end of mask', tb_sites[k])
                                #if errors emerge add another conditional
                                if line[0] == '-':
                                    print('combine', '-',line[1], k)
                                    print('\t'.join(['-',line[1], str(int(k)-int(line[1]))]))
                                    out.write('\t'.join(['-',line[1], str(int(k)-int(line[1]))])+'\n')

                                #dist = tb_sites[k]-k
                                #prevkey = [k,dist]
                                
                                #print('update', prevkey )
                                #print(f'in a range: MASK {line[1]}')
                                #line[4] = '-'
                                #line[-1] = '1'
                                #should i do this?
                                #missing += 1
                                #print('start', k)
                                #print('pos', line[1])
                                #print('end', ld_sites[k])
                                #print
                                break
                                #check if less than final value in coverage range
                                #indicates that the position should be masked
                                #break after to save time
                                #print()
                                
                            elif int(line[1]) > k and int(line[1])+int(line[2])<tb_sites[k]:
                                print(f'fully overlap: ignore?')
                                #line[4] = '-'
                                #line[-1] = '1'
                                #should i do this?
                                #missing += 1
                                print('start', k)
                                print('pos', line[1])
                                print('end', tb_sites[k])
                                #print
                                break

                                
                            #if greater than first and last value, remove range from dict
                            #time saver
                            else:
                                print('right side overlap?', 'prevkey', prevkey)
                                print('start', k)
                                print('pos', line[1])
                                print('end', tb_sites[k])
                                print('OVERLAP!!!!! end of snp and beginnin of mask', k, line, tb_sites[k])
                                #print('end of snp', int(line[1])+int(line[2]))
                                print('start of mask', k)
                                print( 'start of snp', int(line[1]))
                                print('end of mask', tb_sites[k])
                                #if errors emerge add another conditional
                                #there should be no other possibilities...
                                if line[0] == '-':
                                #    print('combine', '-',line[1], k)
                                #    print('\t'.join(['-',line[1], str(int(k)-int(line[1]))]))
                                    
                                    print('mask up to position')
                                    print('\t'.join(['-',str(k), str(int(line[1])-k)]))
                                    out.write('\t'.join(['-',str(k), str(int(line[1])-k)])+'\n')
                                    print('position after mask')
                                    print('\t'.join(['-',line[1], line[2]])+'\n')
                                    out.write('\t'.join(['-',line[1], line[2]])+'\n')
                                    dist = tb_sites[k]-k
                                    prevkey = [k,dist]
                                
                                    print('update', prevkey )
                                    tb_sites.pop(k)
                                else:
                                    raise ValueError(f"Indel needs to be masked at {line[1]}")
                                #distance = int(tb_sites[k])-k
                                #print('dists', distance)
                                #out.write('\t'.join(['-',str(k), str(distance)])+'\n')
                                #tb_sites.pop(k)
                                #print('len dict', len(tb_sites))
                                
                                #need to set up value for masks
                            
                        
                            elif int(line[1]) < k:
                                #print('range', k, tb_sites[k])
                                #print('no mask?', line[1])
                                #print('\t'.join(line))
                                
                                out.write('\t'.join(line)+'\n')
                                
                                #alts = line[4].split(',')
                                #alt = alts[int(var)-1]
                                #line[4] = alt
                                #line[-1] = '1'
                                
                                break
                            
                        #dist = tb_sites[k]-k
                        #prevkey = [k, dist]
                else:
                    out.write(line)
                    '''
                    

'''
def mask_LDsites(cf, cd, diff,sample):
    #print('TB mask regions')
    ld_sites = mask_low_depth(cf, cd)
    #print(ld_sites)
    #keys = sorted(tb_sites.keys())
    #print('tb sites',len(tb_sites))
    #for i in tb_sites:
    #    print(i)
    #print('MASK TB SITES', diff_file)
    #print(sample)
    
    with open(diff) as df:
        with open(f'{wd}{sample}.masked.diff', 'w') as out:
            for line in df:
                #print(line)
                keys = sorted(ld_sites.keys())
                if not line.startswith('>'):
                    line = line.strip().split()
                    #print('position', line[1])

                    for k in keys:
                        #print('k', k)
                        #print('k', k, 'tb_sites', tb_sites[k])
                        #check if greater than first value in coverage range
                        if int(line[1]) >= k:
                            #print('k', line[1], k)
                            #check if less than final value in coverage range
                            #indicates that the position should be masked
                            #break after to save time
                            #print()
                            
                            if int(line[1]) < ld_sites[k]:
                                #print(f'in a range: MASK {line[1]}')
                                #line[4] = '-'
                                #line[-1] = '1'
                                #should i do this?
                                #missing += 1
                                #print('start', k)
                                #print('pos', line[1])
                                #print('end', ld_sites[k])
                                #print
                                break
                            
                            #if greater than first and last value, remove range from dict
                            #time saver
                            else:
                                #print('delete?')
                                #print('start', k)
                                #print('pos', line[1])
                                #print('end', tb_sites[k])
                                distance = int(ld_sites[k])-k
                                #print('dists', distance)
                                out.write('\t'.join(['-',str(k), str(distance)])+'\n')
                                ld_sites.pop(k)
                                #print('len dict', len(tb_sites))
                                
                                #need to set up value for masks
                        
                        elif int(line[1]) < k:
                            #print('range', k, tb_sites[k])
                            #print('no mask?', line[1])
                            out.write('\t'.join(line)+'\n')
                            #alts = line[4].split(',')
                            #alt = alts[int(var)-1]
                            #line[4] = alt
                            #line[-1] = '1'
                            break
                        
        
                else:
                    out.write(line)
                    '''
                    
                    
                
#option to write output to file
def vcf_to_diff(vcf_file, output):
    #takes a single sample vcf and converts to diff format 
    lines = []
    with open(vcf_file, 'rt') as v:
        #with open(output, 'w') as o:
        missing = 0
        total = 0
        for line in v:
            if not line.startswith('##'):
                if line.startswith('#'):
                    line = line.strip().split()
                    sample = line[-1]
                    #print('sample', sample)
                    #o.write(f'>{sample}\n')
                    lines.append([f'>{sample}'])
                else:
                    total += int(len(line[3]))
                    line = line.strip().split()
                    var = line[-1]
                    
                    if var != '0/0':
                        
                        if var == './.':
                            #print('missing', line)
                            #missing += int(len(line[3]))
                            line[4] = '-'
                            line [-1] = '1'

                        else: 
                            
                            var = var.split('/')
                            var = var[0]
                            alts = line[4].split(',')
                            alt = alts[int(var)-1]
                            line[4] = alt
                            line[-1] = '1'

                        #print('line', line)
                        assert type(line[4]) == str
                        if len(line[3]) == 1:

                            if len(line[4]) == 1:
                                #print(line)
                                #o.write('\t'.join([line[4],line[1], '1'])+'\n')
                                lines.append([line[4],line[1], '1'])
                            #elif len(line[4]) > 1:
                                #print('insertion')


                        elif len(line[3]) > 1:
                            if len(line[4]) == len(line[3]):
                                #print(line)
                                newlines = find_snps(line)
                                for n in newlines:
                                    #o.write('\t'.join(n)+'\n')
                                    lines.append(n)

                            elif len(line[4]) == 1:
                                #print('deletion')
                                newline = process_dels(line)
                                #print(newline)
                                #o.write('\t'.join(newline)+'\n')
                                lines.append(newline)

                            else:
                                newline = process_others(line)
                                #o.write('\t'.join(newline)+'\n')
                                lines.append(newline)
    #figure out how to count missing data and add stats here 
    return lines

                                    #print('indel', line)
                                #print(line)
            #print('missing', missing, 'total', total)

    #return sample, missing, total
    #return sample, missing, total

def squish(file):
    #condenses diff entries that can be 
    with open(file) as f:
        with open(f'{file}squish', 'w') as o:

            prev = None
            for line in f:
                if not line.startswith('>'):

                    #process line
                    line = line.strip().split()
                    #print(line, prev)
                    
                    if prev != None:
                        #print('line', line)
                        #print('prev', prev)
                        if prev[0] == line[0]:
                            #print('prev', prev, 'line', line)
                            if int(prev[1]) == int(line[1])-int(prev[2]):
                                #print('prev', prev, 'now', line)
                                prev[2] = str(int(prev[2])+int(line[2]))
                                #print('new prev', prev)
                            else:
                                o.write('\t'.join(prev)+'\n') 
                                prev = line
                        else:
                            o.write('\t'.join(prev)+'\n') 
                            prev = line
                        
                    
                    #the first line of file becomes prev variable 
                    else:
                        prev = line
                        #print('first prev!', prev)
                        
                        
            

                #write header to new file
                else:
                    o.write(line)

def make_files(lenRow, samps,wd):
    files = {}
    for i in range(len(samps)):
        #print(i)
        s = samps[i]
        #print(s)
        if '/' not in s:
            file = open(f'{wd}{s}.vcf','w')
        else:
            newname = s.replace('/', '-')
            #print(newname)
            file = open(f'{wd}{newname}.vcf','w')
        #files.append(file)
        files[i+9] = file
        #print(files)
    return files

                            
def count_samples(vcf):
    #open file to indentify number of samples 
    with open(vcf, 'r') as v:
        for line in v:
            #if the line is not part of the heading
            if not line.startswith('##'):

                #if the line contains the column names 
                #note that this assumes 9 meta data columns, if you need to account for more/less change this 
                if line.startswith('#'):
                    line = line.strip().split('\t')
                    samps = line[9:]
                    lenRow = len(line)

                    #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                    # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                    #print('lenRow',lenRow)
                    break
    
    return lenRow, samps                           


def count_samples_bin(vcf):
    #open file to indentify number of samples 
    with gzip.open(vcf, 'rt') as v:
        for line in v:
            #if the line is not part of the heading
            if not line.startswith('##'):

                #if the line contains the column names 
                #note that this assumes 9 meta data columns, if you need to account for more/less change this 
                if line.startswith('#'):
                    line = line.strip().split('\t')
                    samps = line[9:]
                    lenRow = len(line)

                    #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                    # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                    #print('lenRow',lenRow)
                    break
    
    return lenRow, samps

def read_VCF(vcf, files):
    #count = 0 
    with open(vcf) as v:
        for line in v:
            if not line.startswith('##'):
                line = line.strip().split()
                position = line[0:9] 
                #print('position', position)
                #print(len(line))
                #count = 1
                for f in files:
                    #print('f',f)
                    #print(files[f])
                    #count += 1
                    parcel = [line[f]]
                    #print('parcel',parcel)
                    #print(parcel)
                    newline = position+parcel
                    #print(newline)
                    #print(len(newline))
                    files[f].write('\t'.join(newline)+'\n')
                    #print(position[1],f,files[f][0])
                    #print(f)
                    
                #count += 1

            
                
            
            else:
                #print(line)
                for f in files:
                    files[f].write(line)

def read_VCF_bin(vcf, files):
    #count = 0 
    with gzip.open(vcf, 'rt') as v:
        for line in v:
            if not line.startswith('##'):
                line = line.strip().split()
                position = line[0:9] 
                #print('position', position)
                #print(len(line))
                #count = 1
                for f in files:
                    #print('f',f)
                    #print(files[f])
                    #count += 1
                    parcel = [line[f]]
                    #print('parcel',parcel)
                    #print(parcel)
                    newline = position+parcel
                    #print(newline)
                    #print(len(newline))
                    files[f].write('\t'.join(newline)+'\n')
                    #print(position[1],f,files[f][0])
                    #print(f)
                    
                #count += 1

            
                
            
            else:
                #print(line)
                for f in files:
                    files[f].write(line)




#SCRIPT STARTS HERE

binary = True
with gzip.open(vcf, 'r') as test:
    try:
        test.read(1)
    except OSError:
        #print('not binary')
        binary = False

if binary == True:
    lenRow, samps = count_samples_bin(vcf)
else:
    lenRow, samps = count_samples(vcf)

#be careful w dictionaries!!!
files = make_files(lenRow, samps, wd)

if binary == True:
    read_VCF_bin(vcf, files)
    
else:
    read_VCF(vcf, files)
    
masks = mask_TB(tbmf)
#this is not parallelized, the more samples in the vcf the longer this will take
#print(files)
#print(masks)

for f in files:
    files[f].close()
    #print(f, files[f])
    
    sample = os.path.basename(files[f].name)[:-4]
    filepath = files[f].name
    #print('sample', sample)
    #print('filepath', filepath)

    
    #note: theres no point to do this since the files are already completely filling the space
    #    os.system(f'bgzip -f {filepath}')
    #    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}.gz")
    #    os.system(f"rm {filepath}.gz")
    
    
    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}")
    os.system(f"rm {filepath}")
    
    #note: theres no point to do this since the files are already completely filling the space
    #os.system(f'bgzip -f {filepath}.filt')
    #vcf_to_diff(f'{filepath}.filt.gz', f'{wd}{sample}.diff')

    #if there is a provided coverage file it will be used to mask low coverage (less than cd) regions 
    #note that only one coverage file can be provided and it will result in an error if the vcf has more samples than coverage files 

    #more work on this later 
    if cf == None:
        lines = vcf_to_diff(f'{filepath}.filt', f'{wd}{sample}.diff')
        os.system(f'rm {filepath}.filt')
        #squish(f'{wd}{sample}.diff')
        #os.system(f'mv {wd}{sample}.diffsquish {wd}{sample}.diff')
    
    #figure out how to count missing samples, mask full low coverage regions

    #take out low depth masks see if it works 
    else:
        print('FILTERING FOR LOW COVERAGE')
        assert len(samps) == 1, f'must only have one sample for each coverage file'
        lines = vcf_to_diff(f'{filepath}.filt', f'{wd}{sample}.diff')
        os.system(f'rm {filepath}.filt')
        #squish(f'{wd}{sample}.diff')
        #os.system(f'mv {wd}{sample}.diffsquish {wd}{sample}.diff')
        #all_masks = condense_mask_regions(cf,cd,tbmf)
        #for m in all_masks:
        #    print(m, all_masks[m])
        #print(len(all_masks))
        #mask_LDsites(cf,cd,f'{wd}{sample}.diff',sample)
        #os.system(f'mv {wd}{sample}.masked.diffsquish {wd}{sample}.diff')
    

    #os.system(f"rm {filepath}.filt")
    #mask_TBsites(f'{wd}{sample}.diff', tbmf, sample)
    #squish(f'{wd}{sample}.masked.diff')
    #os.system(f'mv {wd}{sample}.masked.diffsquish {wd}{sample}.masked.diff')
    #os.system(f'rm {wd}{sample}.diff')
    

#SCRIPT ENDS HERE



#ld = mask_low_depth(cf,cd)

with open('testlines','w') as o:
    for line in lines:
        o.write('\t'.join(line)+'\n')
        #print(line)

    

    








