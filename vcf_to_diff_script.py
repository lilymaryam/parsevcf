import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--VCF', required=True, type=str,help='path to VCF to be processed')
parser.add_argument('-d', '--working_directory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-tbmf', '--tb_maskfile', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
#parser.add_argument('-z', '--zip_output', required=False,type=bool, default = False, help="zip output files (default false")
parser.add_argument('-cf', '--bed_coverage_file', required=False, type=str, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")
parser.add_argument('-cd', '--coverage_depth', required=False, default=10, type=int, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")

#parser.add_argument('-sd', '--scriptsDirectory', type=str, required=True, help='path to directory where scripts are')
#parser.add_argument('-od', '--outDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')

#NOTICE: Do not use this script if you have space limitations


args = parser.parse_args()
vcf = args.VCF
wd = args.working_directory
tbmf = args.tb_maskfile
cf = args.bed_coverage_file
cd = args.coverage_depth
#z = args.zip_output

#print('z', z)
#od = args.outDirectory

#makes sure input path wont cause error
if wd[-1] != '/':
    #print('add a final slash')
    wd = wd+'/'
#if od[-1] != '/':
    #print('add a final slash')
#    od = od+'/'
#sd = args.scriptsDirectory
#if sd[-1] != '/':
    #print('add a final slash')
#    sd = sd+'/'

#t = args.threads

                            
def find_snps(line):
    #for lines where len(ref)==len(alt), look for snps instead of processing as one large chunk

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

def process_dels(line):
    #for lines where len(ref) > 1 AND len(alt) == 1, these deletions are easily processed 
    #currently deletions and missing data are all converted to '-'
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
        #this can be changed 
        l = ['-', line[1], str(len(line[3]))]
        
    return l



def process_others(line):
    #for sitations where reference and alt do not align 
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    return l

#read mask file, send mask regions to vcf-to-diff
def mask_TB(tbmf):
    tb_sites = {}
    with open(tbmf) as file:
        for line in file:
            line=line.strip().split()
            tb_sites[int(line[1])] = int(line[2])
    #print(sites)
    return tb_sites

def mask_low_depth(cf, cd):
    ld_sites = {}
    with open(cf) as cf:
        for line in cf:
            #for currect bed coverage file 
            if line.startswith('NC_000962.3'):
                line = line.strip().split()
                if int(line[3]) < cd:
                    ld_sites[int(line[1])] = int(line[2])

    return ld_sites

def mask_TBsites(diff_file, tbmf, sample):
    print('TB mask regions')
    tb_sites =  mask_TB(tbmf)
    #keys = sorted(tb_sites.keys())
    #print(tb_sites)
    #print('MASK TB SITES', diff_file)
    #print(sample)
    with open(diff_file) as df:
        with open(f'{sample}masked.diff', 'w') as out:
            for line in df:
                keys = sorted(tb_sites.keys())
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
                            
                            if int(line[1]) < tb_sites[k]:
                                print(f'in a range: MASK {line[1]}')
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
                                tb_sites.pop(k)
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
    tb_sites = mask_TB(tbmf)
    keys = sorted(tb_sites.keys())
    print('position', line[1])
    for k in keys:
        #print('line', line)
        print('k', k)

        #check if greater than first value in coverage range
        if int(line[1]) >= k:
            #check if less than final value in coverage range
            #indicates that the position should be masked
            #break after to save time
            if int(line[1]) < ld_sites[k]:
                print(f'in a range: MASK {line[1]}')
                line[4] = '-'
                line[-1] = '1'
                missing += 1
                #print('start', k)
                #print('pos', line[1])
                #print('end', ld_sites[k])
                break
            '''
            



def vcf_to_diff(vcf_file, output):
    #takes a single sample vcf and converts to diff format 

    with open(vcf_file, 'rt') as v:
        with open(output, 'w') as o:
            missing = 0
            total = 0
            for line in v:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.strip().split()
                        #sample = line[-1]
                        #print('sample', sample)
                        o.write(f'>{sample}\n')
                    else:
                        total += 1
                        line = line.strip().split()
                        var = line[-1]
                        
                        if var != '0/0':
                            
                            if var == './.':
                                #print('missing', line)
                                missing += 1
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
                                    o.write('\t'.join([line[4],line[1], '1'])+'\n')
                                #elif len(line[4]) > 1:
                                    #print('insertion')


                            elif len(line[3]) > 1:
                                if len(line[4]) == len(line[3]):
                                    #print(line)
                                    newlines = find_snps(line)
                                    for n in newlines:
                                        o.write('\t'.join(n)+'\n')

                                elif len(line[4]) == 1:
                                    #print('deletion')
                                    newline = process_dels(line)
                                    #print(newline)
                                    o.write('\t'.join(newline)+'\n')

                                else:
                                    newline = process_others(line)
                                    o.write('\t'.join(newline)+'\n')



                                    #print('indel', line)
                                #print(line)
            #print('missing', missing, 'total', total)

    return sample, missing, total

def vcf_to_diff_pipeline(vcf_file, output, cf, cd):
    #takes a single sample vcf and converts to diff format 
    ld_sites = mask_low_depth(cf,cd)
    #tb_sites = mask_TB(tbmf)
    #print('tb sites', tb_sites)
    with open(vcf_file, 'rt') as v:
        with open(output, 'w') as o:
            missing = 0
            total = 0
            for line in v:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.strip().split()
                        #sample = line[-1]
                        #print('sample', sample)
                        o.write(f'>{sample}\n')
                    else:
                        total += 1
                        line = line.strip().split()
                        var = line[-1]
                        
                        if var != '0/0':
                            
                            if var == './.':
                                #print('missing', line)
                                missing += 1
                                line[4] = '-'
                                line [-1] = '1'

                            else: 
                                
                                var = var.split('/')
                                var = var[0]
                                
                                keys = sorted(ld_sites.keys())
                                for k in keys:
                                    #print('line', line)
                                    #print('k', k)

                                    #check if greater than first value in coverage range
                                    if int(line[1]) >= k:
                                        #check if less than final value in coverage range
                                        #indicates that the position should be masked
                                        #break after to save time
                                        if int(line[1]) < ld_sites[k]:
                                            print(f'low coverage: MASK {line[1]}')
                                            line[4] = '-'
                                            line[-1] = '1'
                                            #should i do this?
                                            missing += 1
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
                                            #print('end', ld_sites[k])
                                            ld_sites.pop(k)
                                            
                                            #need to set up value for masks

                                    elif int(line[1]) < k:
                                        #print('range', k, ld_sites[k])
                                        #print('no mask?', line[1])
                                        alts = line[4].split(',')
                                        alt = alts[int(var)-1]
                                        line[4] = alt
                                        line[-1] = '1'
                                        break

                                    #print(line)
                                        

                                        
                                        #break
                                    #    if int(line[1]) <= ld_sites[s]:
                                    #        print('line', line[1], s, ld_sites[s])
                                    #    else:
                                    #        print('out of range', s, ld_sites[s])
                                        #print('alts', alts)
                                        #print('alt', alt)
                                        
                                

                            #print('line', line)
                            assert type(line[4]) == str
                            if len(line[3]) == 1:

                                if len(line[4]) == 1:
                                    #print(line)
                                    o.write('\t'.join([line[4],line[1], '1'])+'\n')
                                #elif len(line[4]) > 1:
                                    #print('insertion')


                            elif len(line[3]) > 1:
                                if len(line[4]) == len(line[3]):
                                    #print(line)
                                    newlines = find_snps(line)
                                    for n in newlines:
                                        o.write('\t'.join(n)+'\n')

                                elif len(line[4]) == 1:
                                    #print('deletion')
                                    newline = process_dels(line)
                                    #print(newline)
                                    o.write('\t'.join(newline)+'\n')

                                else:
                                    newline = process_others(line)
                                    o.write('\t'.join(newline)+'\n')



                                    #print('indel', line)
                                #print(line)
            print('missing', missing, 'total', total)

    return sample, missing, total


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
for f in files:
    files[f].close()
    sample = os.path.basename(files[f].name)[:-4]
    filepath = files[f].name
    
    #note: theres no point to do this since the files are already completely filling the space
    #    os.system(f'bgzip -f {filepath}')
    #    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}.gz")
    #    os.system(f"rm {filepath}.gz")
    
    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}")
    os.system(f"rm {filepath}")
    
    #note: theres no point to do this since the files are already completely filling the space
    #os.system(f'bgzip -f {filepath}.filt')
    #vcf_to_diff(f'{filepath}.filt.gz', f'{wd}{sample}.diff')
    if cf == None:
        vcf_to_diff(f'{filepath}.filt', f'{wd}{sample}.diff')
    else:
        assert len(samps) == 1, f'must only have one sample for each coverage file'
        #mask_low_depth(cf)
        vcf_to_diff_pipeline(f'{filepath}.filt', f'{wd}{sample}.diff',cf,cd)

    squish(f'{wd}{sample}.diff')
    os.system(f'mv {wd}{sample}.diffsquish {wd}{sample}.diff')
    os.system(f"rm {filepath}.filt")
    mask_TBsites(f'{wd}{sample}.diff', tbmf, sample)

#tb = mask_TB(tbmf)
#print(tb)

    

    








