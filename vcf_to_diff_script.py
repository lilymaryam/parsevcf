import os
import sys
import argparse
import gzip
import logging

parser = argparse.ArgumentParser(description="Tool for converting VCF to diff formats.")
parser.add_argument('-v', '--VCF', required=True, type=str,help='path to VCF to be processed')
parser.add_argument('-d', '--working_directory', required=False, default='./', type=str, help='Output directory for storing intermediate and sample-level outputs.')
parser.add_argument('-mf', '--maskfile', required=False, type=str, help='Masking data in BED format.')
parser.add_argument('-cf', '--bed_coverage_file', required=False, type=str, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")
parser.add_argument('-cd', '--coverage_depth', required=False, default=10, type=int, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")
parser.add_argument("-ha", "--haploid", action='store_true', help='Indicate whether the VCF is diploid or haploid. Default is diploid.')
#https://stackoverflow.com/questions/14097061/easier-way-to-enable-verbose-logging
parser.add_argument(
    '-D', '--debug',
    help="Pring debugging statements.",
    action="store_const", dest="loglevel", const=logging.DEBUG,
    default=logging.WARNING,
)
parser.add_argument(
    '-V', '--verbose',
    help="Verbose output.",
    action="store_const", dest="loglevel", const=logging.INFO,
)
parser.add_argument("-r", "--reference_length", default=4411532, type=int, help="Length of the reference genome. Default is 4411532 (length of the TB genome).")
parser.add_argument("-rn", "--reference_name", default='NC_000962.3', help='Name of the reference genome. Default is "NC_000962.3", the TB reference.')
parser.add_argument('-t', "--tag_filter", action='store_true', help='Use to filter GT tags in input vcfs. Requires BCFtools.')
parser.add_argument("-sp", "--split_diff", action='store_true', help='Use to create separate .diff files for each of your input samples.')
parser.add_argument("-mr", "--mask_report", action='store_true', help='Print sample level text reports of masking information.')
parser.add_argument("-od", "--output_diff", help='Name of a compiled output diff.')
args = parser.parse_args()
vcf = args.VCF
wd = args.working_directory
tbmf = args.maskfile
cf = args.bed_coverage_file
cd = args.coverage_depth
refname = args.reference_name
logging.basicConfig(level=args.loglevel)

len_ref = args.reference_length

#makes sure input path wont cause error
if wd[-1] != '/':
    wd = wd+'/'
#create the directory if it does not exist.
if not os.path.exists(wd):
    os.mkdir(wd)

def find_snps(line):
    '''
    for lines where len(ref)==len(alt), look for snps instead of processing as one large chunk 
    Args: 
        line: a list containg the line from the VCF
    Returns:
        lines: a list of lists of lines to be added to the diff file 
    '''
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
    '''
    for lines where len(ref) > 1 AND len(alt) == 1 (currently deletions and missing data are all converted to '-')
    Args: 
        line: a list containing info from a line of the VCF
    Returns:
        l: a list containing the diff-formatted version of the deletion
    '''
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
 
def process_others(line):
    '''
    when reference and alt do not align (are differenct lens), mask entire reference 
    Args:
        line: a list containing info from a line of the VCF
    Input: 
        l: a list containing the diff-formatted version of the line
    '''
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    return l

def mask_TB(tbmf):
    '''
    Opens and reads bed file containing universally ignored TB positions
    note: bed coverage file is 0 indexed, so 1 added to everything 
    Args:
        tbmf: bed file with positions to be ignored (note positions are assumed to be 0 indexed)
    Returns:
        tb_sites: dictionary where key is start of masked region (1 index) and value is end of masked region (not inclusive)
    '''
    tb_sites = {}
    with open(tbmf) as file:
        for line in file:
            line=line.strip().split()
            tb_sites[int(line[1])+1] = int(line[2])+1
    return tb_sites

def mask_low_depth(cf, cd):
    '''
    read bed coverage file and generate sites to be masked 
    note: if coverage does not have HR37c reference it will throw an error (this can be changed)
    note: bed files are 0-indexed in col1 and 1-indexed in col2, i am adding one to both to make them both one indexed
    Args: 
        cf: path to bed coverage file 
        cd: integer indicating coverage depth needed 
    Returns:
        ld_sites: dictionary of low-depth sites needing to be masked 
    '''
    ld_sites = {}
    prev = None
    with open(cf) as cf:
        for line in cf:
            #might need to delete or change this
            #for currect bed coverage file 
            if line.startswith(refname):
                line = line.strip().split()
                #re-index to match VCF
                line[1] = str(int(line[1])+1)
                line[2] = str(int(line[2])+1)
                #if the coverage is below cd
                if int(line[3]) < cd:
                    if prev == None:
                        #for first low-coverage line in file
                        ld_sites[int(line[1])] = int(line[2])
                        prev = [int(line[1]), int(line[2])]
                    else:
                        #for all subsequent low-coverage lines, determine if they can be combined 
                        #if start of line is the same as the end of prev
                        if int(line[1]) == prev[1]:
                            #combine sites and update prev
                            ld_sites[prev[0]] = int(line[2])
                            prev[1] = int(line[2])
                        #create a new site and update prev
                        else:
                            ld_sites[int(line[1])] = int(line[2])
                            prev = [int(line[1]), int(line[2])]
    #this conditional might need to be fixed if there are no low-coverage areas
    if ld_sites == {}:
        raise Exception('coverage file contains no low-coverage areas- check input.')

    return ld_sites

"""
NOT FOR USE WITH UNIVERSAL MASK2REF
def check_prev_mask(prev, line):
    '''
    when merging ld and tb masks, make sure the masks are not overlapping with previously added masks
    Args:
        prev: a list of the start and stop of the previously added mask region 
        line: a list of the start and stop of the to-be-added mask region
    Returns:
        overlap: a boolean meant to indicate if prev and line overlap
        change: a list containing important information for updating the prev value
    '''
    #currently not checking overlap to left of prev bc that indicates a bigger error
    #may need to change?
    overlap = False
    change = None
    
    prev_s = int(prev[0])
    prev_e = int(prev[1])
    line_s = int(line[0])
    line_e = int(line[1])

    '''
    DEBUGGING
    logging.debug('prev', prev_s, prev_e, 'line', line_s, line_e)
    '''

    #NO OTHER CONDITIONALS NEEDED BC LINE CANT BE TO THE LEFT OF PREV AND IF LINE HAS NO OVERLAP W PREV NO TRACKING IS NEEDED 
    #add error checking to make sure line isn't to left of prev?

    #if line is completely contained inside prev, line should not be added to all_sites 
    if line_s >= prev_s and line_e <= prev_e:
        overlap = True
    #if line is overlapping the right end of prev
    elif line_s >= prev_s and line_s <= prev_e and line_e >= prev_e:
        overlap = True 
        prev[1] = line_e
        change = prev
    return overlap, change


def condense_mask_regions(ld_sites,tb_sites):
    '''
    for instances with low-depth masking, combine low-depth masks and universal masks into a single data structure 
    Args: 
        ld_sites: a dictionary containing low-depth sites to be masked
        tb_sites: a dictionary containing universal masking sites 
    Returns:
        all_sites: a dictionary containing all masking sites from both universal and low-depth
    '''
    #editing thought: would likely benefit from being a list rather than a dictionary 
    tb_keys = sorted(tb_sites.keys())
    ld_keys = sorted(ld_sites.keys())
    all_sites = {}
    tb_keys_ind = 0
    ld_keys_ind = 0
    #cont = 0

    #track previous line in all_sites
    prev = None
    #iterate through tb_keys and ld_keys exactly 1 time, track index of keys as you iterate
    while tb_keys_ind < len(tb_keys) or ld_keys_ind < len(ld_keys):
        #if still iterating through both lists
        if tb_keys_ind < len(tb_keys) and ld_keys_ind < len(ld_keys): 
            tb_start = tb_keys[tb_keys_ind]
            tb_end =  tb_sites[tb_keys[tb_keys_ind]]
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]

            if ld_start < tb_start:
                if ld_end < tb_start:
                    #the low-depth site is completely to the left of the universal site
                    if prev != None:
                        #before adding to all_sites, make sure it doesn't overlap prev
                        overlap, change = check_prev_mask(prev, [ld_start, ld_end])
                        #if overlap is detected AND requires prev to be updated
                        if overlap == True and change != None:
                            #update prev in all_sites
                            all_sites[change[0]] = change[1]
                    if prev == None or overlap == False:
                        #add new site to all_sites
                        all_sites[ld_start] = ld_end
                    #update ld index because site was processed
                    ld_keys_ind += 1
                elif ld_end >= tb_start:
                    #this accounts for ld overlapping tb on the left
                    if prev != None:
                        overlap, change = check_prev_mask(prev, [ld_start, tb_end])
                        if overlap == True and change != None:
                            #if new region overlaps, update prev
                            all_sites[change[0]] = change[1]
                    #if there is no overlap, add new region to all_sites
                    if prev == None or overlap == False:
                        all_sites[ld_start] = tb_end
                    #since both regions are added at the same time, update index for both lists
                    tb_keys_ind += 1
                    ld_keys_ind += 1
            
            elif ld_start <= tb_end and ld_end > tb_end:
                #if ld overlaps tb on the right
                if prev != None:
                    #before adding to all_sites, make sure it doesn't overlap prev
                    overlap,change = check_prev_mask(prev, [tb_start, ld_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                #if there is no overlap, add new region to all_sites
                if prev == None or overlap == False:
                    all_sites[tb_start] = ld_end
                #since both regions are added at the same time, update index for both lists
                tb_keys_ind += 1
                ld_keys_ind += 1

            elif ld_start > tb_end:
                #if there is no overlap between ld and tb, and ld is on the right
                if prev != None:
                    #before adding to all_sites, make sure it doesn't overlap prev
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    #if there is no overlap, add new region to all_sites
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1

            elif ld_start >= tb_start and ld_end <= tb_end:
                # if tb and ld fully overlap with ld inside
                if prev != None:
                    #make sure tb doesnt overlap prev 
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != False:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            
            elif ld_start <= tb_start and ld_end >= tb_end:
                #if ld and tb fully overlap with tb inside
                if prev != None:
                    #make sure ld doesnt overlap prev
                    overlap,change = check_prev_mask(prev, [ld_start, ld_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[ld_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1
                
            #covered all 6 possible positions of the two regions
            '''
            DEBUG
            logging.debug('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
            '''

        #possible bug: prev overlaps with one list the first time the other list expires
        elif tb_keys_ind >= len(tb_keys) and ld_keys_ind < len(ld_keys):
            #if reach end of tb_masks process ld only
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]
            all_sites[ld_start] = ld_end
            ld_keys_ind += 1

        #possible bug: prev overlaps with one list the first time the other list expires
        elif ld_keys_ind >= len(ld_keys) and tb_keys_ind < len(tb_keys):
            #if reach end of ld masks process tb only
            tb_start = tb_keys[tb_keys_ind]
            tb_end = tb_sites[tb_keys[tb_keys_ind]]
            all_sites[tb_start] = tb_end
            tb_keys_ind += 1 

        #keep all sites sorted
        if len(all_sites) > 0:
            all_sites_keys = sorted(all_sites.keys())
        #track last site for every iteration
        prev = [all_sites_keys[-1],all_sites[all_sites_keys[-1]]]
    return all_sites
    """
                    
def squish(lines):
    '''
    condenses diff lines that can be compressed into a single line
    Args:
        lines: a list of diff-formatted lines each stored as a list
    Outputs:
        newlines: a list of diff-formatted lines after compression
    '''
    #track previous line in lines
    prev = None
    #create a new list of lines for after compression
    newLines = []
    for line in lines:
        #skip header 
        if not line[0].startswith('>'):
            #if not the first line in the file
            if prev != None:
                #if prev and line have the same nucleotide
                if prev[0] == line[0]:
                    #if end of prev overlaps w beginning of line
                    if int(prev[1]) == int(line[1])-int(prev[2]):
                        #rewrite prev and line into a new prev
                        prev[2] = str(int(prev[2])+int(line[2]))
                    else:
                        #if prev and line don't overlap, add prev to newLines and update prev
                        newLines.append(prev) 
                        prev = line
                else:
                    #if prev and line don't overlap, add prev to newLines and update prev
                    newLines.append(prev) 
                    prev = line
            
            #the first line of file becomes prev variable 
            else:
                prev = line
        #write header to new file
        else:
            newLines.append(line)
    #add last prev to end of newLines
    if prev != None:
        newLines.append(prev)  
    else:
        #prev should probably not be None
        logging.warning('no lines in file?')
    #should i overwrite lines variable for storage consideration?
    return newLines

def vcf_to_diff(vcf_file):
    '''
    takes a single sample vcf and converts to diff format
    NOTE: this function makes the assumption that incoming vcf file is genotyped as diploid
          (it's common practice to call variants on TB as if it were diploid)
    Args: 
        vcf_file: uncompressed single sample vcf 
    Outputs:
        newlines: a list of diff-formatted lines for the file
    ''' 
    lines = []
    with open(vcf_file, 'rt') as v:
        #missing = 0
        total = 0
        for line in v:
            #ignore header lines
            if not line.startswith('##'):
                #find column names
                if line.startswith('#'):
                    line = line.strip().split()
                    #last column name is single-sample vcf will be sample name
                    sample = line[-1]
                    #make diff file header
                    lines.append([f'>{sample}'])
                #all position lines
                else:
                    #total += int(len(line[3]))
                    line = line.strip().split()
                    #genotype
                    var = line[-1]

                    # for all lines with multiple alt alleles
                    #if len(line[4].split(','))>1:
                    #    logging.debug('OPTIONS!!!!', line[4].split(','))

                    #combine ref allele and alt alleles
                    alleles = [line[3]] + line[4].split(',')
                    
                    if args.haploid:
                        #make it a fake homozygous diploid. hack solution.
                        #TODO: something smarter than this.
                        var = var + "/" + var

                    #if genotype is not reference allele
                    if var != '0/0':
                        #logging.debug("Not a reference allele")
                        #logging.debug('line',line)
                        #logging.debug('alleles', alleles)

                        #split genotype to check for heterozygosity
                        genos = var.split('/')
                        
                        if var == './.':
                            #potentially useful to track number of positions with missing info 
                            #missing += int(len(line[3]))
                            line[4] = '-'
                            line[-1] = '1'
                            #logging.debug("Missing info")
                            #logging.debug('missing', line)

                        #assumes diploid genotype
                        #if genotype is heterozygous and reference position is not and indel
                        #NOTE: may need to change this later when indels are not ignored by usher 
                        elif genos[0]!= genos[1] and len(line[3])==1:
                            #logging.debug("Hetero")
                            #logging.debug('HETERO', line)
                            #logging.debug('genos', genos)
                            #logging.debug('alleles', alleles)
                            
                            IUPAC = {
                                'R':['A','G'], 
                                'Y':['C','T'],
                                'S':['C','G'],
                                'W':['A','T'],
                                'K':['G','T'],
                                'M':['A','C']   
                                    }
                            
                            #generated sorted list of both alleles genotyped 
                            vars = sorted([alleles[int(genos[0])], alleles[int(genos[1])]])
                            logging.debug('vars', vars)
                            #if the heterozygous position is a SNP, replace with an IUPAC symbol
                            if len(vars[0])==len(vars[1])==1:
                                logging.debug('SNP')
                                for key in IUPAC:
                                    if IUPAC[key] == vars:
                                        logging.debug('IUPAC key', key)
                                        #alts = line[4].split(',')
                                        #alt = alts[int(var)-1]
                                        line[4] = key
                                        line[-1] = '1'
                                        logging.debug('line after ', line)
                                        break

                            ##if one of the vars is an indel, mask the position       
                            else:
                                logging.info('one of alleles is an indel, mask the ref for clarity')
                                logging.debug(line)
                                line[4] = '-'
                                line [-1] = '1'
                                logging.debug('after', line)

                        #if reference is an indel and/or genotype is homozygous
                        else: 
                            var = var.split('/')
                            var = var[0]
                            alts = line[4].split(',')
                            alt = alts[int(var)-1]
                            line[4] = alt
                            line[-1] = '1'

                        #after above processing there should only be a string with one alt allele
                        assert type(line[4]) == str

                        #if len of ref position and len of alt are both one, process as a SNP
                        if len(line[3]) == 1:
                            if len(line[4]) == 1:
                                lines.append([line[4],line[1], '1'])
                            #if len(line[4]) > 1, the position is an insertion which will not be included in the file

                        elif len(line[3]) > 1:
                            if len(line[4]) == len(line[3]):
                                #if the ref and alt are both longer than 1 but equal to each other,
                                #search through alt for snps
                                newlines = find_snps(line)
                                for n in newlines:
                                    lines.append(n)

                            elif len(line[4]) == 1:
                                #if ref is >1 and alt=1, process line as a simple deletion
                                newline = process_dels(line)
                                lines.append(newline)

                            else:
                                #if len(ref) and len(alt) are both greater than 1 but not the same len as each other
                                newline = process_others(line)
                                lines.append(newline)
    
    #compress adjacent diff lines where possible 
    newLines = squish(lines)
    return newLines

def make_files(samps,wd):
    '''
    Initializes and opens VCF file for each sample in the VCF (will create smaller VCFs for multi-sample VCFs 
    and create a separate editable temp VCF from the original input)
    Args: 
        samps: a list of sample names
        wd: a working directory where all files should be created *make sure directory contains enough space for all uncompressed 
        sampleVCFs 
    outputs:
        files: a dictionary that has the column number as key and the VCF filepath as the value
    '''
    files = {}
    for i in range(len(samps)):
        s = samps[i]
        #replace '/' in sample name with '-'
        if '/' not in s:
            file = open(f'{wd}{s}.vcf','w')
        else:
            newname = s.replace('/', '-')
            file = open(f'{wd}{newname}.vcf','w')
        files[i+9] = file
    return files
                      
def count_samples(vcf):
    '''
    opens VCF and determines how many samples it has (note that this assumes 9cols of metadata)
    *note that VCF is 1-indexed and will remain as such
    Args: 
        vcf: uncompressed VCF containing >=1 samples
    returns:
        lenRow: an int that determines number of columns in VCF (including metadata)
        samps: a list of sample names from the VCF
    '''
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
                    break
    return lenRow, samps                           


def count_samples_bin(vcf):
    '''
    opens VCF and determines how many samples it has (note that this assumes 9cols of metadata)
    *note that VCF is 1-indexed and will remain as such
    Args: 
        vcf: compressed VCF containing >=1 samples
    returns:
        lenRow: an int that determines number of columns in VCF (including metadata)
        samps: a list of sample names from the VCF
    '''
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
                    break
    return lenRow, samps

def read_VCF(vcf, files):
    '''
    open uncompressed VCF and separate columns into individual sample VCFs 
    *note slightly redundant if VCF only contains one sample but allows for editing of VCF wo changing 
    original file
    Args:
        vcf: uncompressed VCF 
        files: dictionary of VCF files to fill in with VCF info
    Outputs:
        none
    '''
    with open(vcf) as v:
        for line in v:
            #for lines not in heading
            if not line.startswith('##'):
                line = line.strip().split()
                #position data for all samples
                position = line[0:9] 
                for f in files:
                    #sample specific data
                    parcel = [line[f]]
                    newline = position+parcel
                    #write position and sample data to file in VCF format
                    files[f].write('\t'.join(newline)+'\n')
            #write position and sample data to file in VCF format
            else:
                #write heading to new VCF file 
                for f in files:
                    files[f].write(line)

def read_VCF_bin(vcf, files):
    '''
    open compressed VCF and separate columns into individual sample VCFs 
    *note: only reads file once regardless of number of samples
    *note slightly redundant if VCF only contains one sample but allows for editing of VCF wo changing 
    original file
    Args:
        vcf: compressed VCF 
        files: dictionary of VCF files to fill in with VCF info
    Outputs:
        none
    '''
    with gzip.open(vcf, 'rt') as v:
        for line in v:
            #for lines not in heading
            if not line.startswith('##'):
                line = line.strip().split()
                #position data for all samples
                position = line[0:9] 
                for f in files:
                    #sample specific data
                    parcel = [line[f]]
                    newline = position+parcel
                    #write position and sample data to file in VCF format 
                    files[f].write('\t'.join(newline)+'\n')
            else:
                for f in files:
                    files[f].write(line)

def check_prev_line(prev, line):
    '''
    when merging lines and masks, make sure all newly added lines are not overlapping with previously added ones
    Args:
        prev: a list containing the most recent added line to all_lines
        line: a list containing the line to be added next to all_lines
    Outputs:
        overlap: a boolean meant to indicate if prev and line overlap
        change: a list containing important information for updating the prev value
    '''
    # NOTE currently not checking overlap to left of prev bc that indicates a bigger error
    overlap = False
    change = None
    newline = None

    '''
    #DEBUG
    logging.debug('masking prev ', prev)
    logging.debug('masking line', line)
    '''

    
    
    prev_s = int(prev[1])
    prev_e = prev_s + int(prev[2])
    line_s = int(line[1])
    line_e = line_s + int(line[2])

    '''
    DEBUG logging.debugS
    logging.debug('prev', prev)
    logging.debug('line', line)
    logging.debug('prev', prev_s, prev_e, 'line', line_s, line_e)
    '''

    if line_s >= prev_s and line_e <= prev_e:
        overlap = True
        logging.info('Full OVERLAP!!!!!')
        logging.debug('prev', prev_s, prev_e, 'line', line_s, line_e)
    elif line_s >= prev_s and line_s < prev_e and line_e >= prev_e:
        overlap = True 
#<<<<<<< fix_mismatch_overlaps
        logging.debug('right overlap vcftodiff', 'prev', prev_s, prev_e, 'line', line_s, line_e)
        
        #logging.debug('prev', prev)
        #logging.debug('line', line)
        #if line_s == prev_e and line[-1]=='1':

        #if line[0] and prev[0] are the same, we can squish these, otherwise, ignore 
        #squish later if necessary 
        if line_s == prev_e:
            logging.debug('start to end', line, prev)
            '''            
        #=======
        logging.debug('right overlap vcftodiff', 'prev', prev, prev_s, prev_e, 'line', line, line_s, line_e)
        #if prev[0] != line[0]:
        #    logging.debug('masking prev ', prev)
        #    logging.debug('masking line', line)
        #logging.debug('prev', prev)
        #logging.debug('line', line)
        logging.debug('prev', prev_s, prev_e, 'line', line_s, line_e, line)
        if line_s == prev_e and line[-1]=='1':
            #logging.debug('SNP')
            #>>>>>>> main
            '''
            overlap = False

        elif line_s < prev_e and line_e>prev_e:
            logging.debug('truly right overlap', 'prev', prev, 'line', line)

            #COME BACK HERE WEDS!!!!!!!!!!! need to figure out how to add two lines 
            if prev[0] == '-' and line[0] != '-':
                logging.debug('masking needed', prev)
                logging.debug('masking needed', line)
                line[1] = str(prev_e)
                line[2] = str(line_e-prev_e)
                logging.debug('newline', line)
                newline = line

            elif prev[0] != '-' and line[0] == '-':
                logging.debug('masking needed complicated', prev)
                logging.debug('masking needed complicated', line)

                
            elif prev[0] == line[0]:
                logging.debug('combine!')
                prev[2] = str(line_e-prev_s)
                logging.debug('prev after ', prev)
                change = prev
            else:
                logging.warning('BAD NEWS', 'line',line, 'prev', prev )
            

            #change = line
        #elif line_s == prev_e and line[-1]!='1':
        #    overlap = True
        #    logging.debug('need to change this!!!', 'prev', prev, 'line', line)
        #    logging.debug('prev')
        else:
            # previously there were arthimatic errors, but this should be fine now
            prev[2] = str(line_e-prev_s)
            logging.debug('prev after ', prev)
            change = prev
    #elif line_s <= prev_s and line_e >= prev_s:
    if newline != None:
        logging.debug('overlap', overlap, 'change', change, 'newlin',newline)
    return overlap, change, newline

#def interpret_overlap()  

def mask_and_write_diff(ld, tb_masks, lines, samps):
    '''
    iterate through masking regions and lines of diff file to mask positions
    Args:
        ld: a dictionary of low depth coverage regions 
        tb_masks: a dictionary of universally masked regions
        lines: a list of diff-formatted lines 
        samps: a list of sample names from the VCF
    Returns:
        all_lines: a list of diff-formatted lines including all of the masked regions
    '''
    
    #if ld != None:
    #    assert len(samps) == 1
    #    masks = condense_mask_regions(ld,tb_masks)
        
    #else:
    #    masks = tb_masks
    if ld != None:
        masks = ld
        masks_key = sorted(ld.keys())
        logging.info("Masking the diff file...")
        logging.debug('masks', masks_key, masks)
    else:
        masks = {}
        masks_key = masks.keys()
    # iterate through all masks and lines one time and combine things as needed
    masks_ind = 0

    #skip header in lines
    lines_ind = 1
    all_lines = [lines[0]]

    prev = None
    while masks_ind < len(masks_key) or lines_ind < len(lines):
        #if both indexes are still going
        if masks_ind < len(masks_key) and lines_ind < len(lines): 
            mask_start = masks_key[masks_ind]
            mask_end =  masks[mask_start]
            line = lines[lines_ind]
            line_start = int(lines[lines_ind][1])
            line_end = int(lines[lines_ind][1])+int(lines[lines_ind][2])

            '''
            DEBUG logging.debugS
            logging.debug('mask_start', mask_start, type(mask_start))
            logging.debug('mask_end', mask_end, type(mask_end))
            logging.debug('line',line)
            logging.debug('line_start', line_start, type(line_start))
            logging.debug('line end', line_end, type(line_end))
            '''

            if line_start >= mask_start and line_end <= mask_end:
                #line and mask fully overlap with line inside
                #logging.debug('full overlap: line inside')
                
                if prev != None:
                    overlap,change,newline = check_prev_line(prev, ['-', mask_start, mask_end-mask_start])
                    if overlap == True and change != None:
                        logging.debug('change', change)
                        #all_sites[change[0]] = change[1]
                    if newline != None:
                        all_lines.append(newline)
                        logging.debug('append!!!!', newline)
                if prev == None or overlap == False:
                    all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])

                masks_ind += 1
                lines_ind += 1

            #need to make sure that if snps overlap they get masked
            elif line_start <= mask_start and line_end >= mask_end:
                #full overlap of line and mask with mask inside

                if prev != None:
                    overlap,change, newline = check_prev_line(prev, line)
                    if overlap == True and change != None:
                        logging.debug('change', change)
                        #all_sites[change[0]] = change[1]
                    if newline != None:
                        all_lines.append(newline)
                        logging.debug('append!!!!', newline)
                if prev == None or overlap == False:
                    
                    all_lines.append(line)
                masks_ind += 1
                lines_ind += 1
        
            elif line_start < mask_start:
                #?
                if line_end <= mask_start:
                    #no overlap, line completely to left
                    #logging.debug(f'line{lines_ind} is below mask{masks_ind}')
                    
                    if prev != None:
                        overlap,change,newline = check_prev_line(prev, line)
                        if overlap == True and change != None:
                            logging.debug('change', change)
                            logging.debug('before',all_lines[-1])
                            #change does this help?
                            
                            all_lines[-1][2] = change[2]
                            logging.debug('after',all_lines[-1])
                        if newline != None:
                            all_lines.append(newline)
                            logging.debug('append!!!!', newline)
                    if prev == None or overlap == False:
                        all_lines.append(line)
                    lines_ind += 1
                
                elif line_end > mask_start:
                    #if line overlaps mask on the left
                    #logging.debug('left overlap')
                    #logging.debug('line',line_start, line_end, 'tb', mask_start, mask_end)
                    
                    if prev != None:
                        overlap,change,newline = check_prev_line(prev, line)
                        if overlap == True and change != None:
                            logging.debug('change', change)
                            #all_sites[change[0]] = change[1]
                        if newline != None:
                            all_lines.append(newline)
                            logging.debug('append!!!!', newline)
                    if prev == None or overlap == False:
                        all_lines.append(['-', str(line_start), str(mask_end-line_end)])
                    masks_ind += 1
                    lines_ind += 1

            elif line_start >= mask_end:
                #if line is completely to the right of mask
                #logging.debug('no overlap', 'line', line_start, line_end, 'mask', mask_start, mask_end)

                if prev != None:
                    #change here
                    overlap,change, newline = check_prev_line(prev, ['-',mask_start,mask_end-mask_start])
                    if overlap == True and change != None:
                        logging.debug('change', change)
                        #all_sites[change[0]] = change[1]
                    if newline != None:
                        all_lines.append(newline)
                        logging.debug('append!!!!', newline)
                if prev == None or overlap == False:
                    all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
                masks_ind += 1
            
            #need to figure out whatn happens if snp is sticking out 
            elif line_start < mask_end and line_end > mask_end:
                #line overlaps mask on the right 
                #logging.debug('right over lap','line',line_start, line_end, 'mask', mask_start, mask_end)

                if prev != None:
                    overlap,change,newline = check_prev_line(prev, line)
                    if overlap == True and change != None:
                        logging.debug('change', change)
                        #all_sites[change[0]] = change[1]
                    if newline != None:
                        all_lines.append(newline)
                        logging.debug('append!!!!', newline)
                if prev == None or overlap == False:
                
                    all_lines.append(['-', str(mask_start), str(line_end-mask_start)])
                masks_ind += 1
                lines_ind += 1

            else:
                logging.debug('other', 'what else could happen?')
                #logging.debug('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                #all_sites[tb_start]


                #logging.debug(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                #logging.debug('ld', ld_start, ld_end)
                #logging.debug('tb', tb_start, tb_end)
            #logging.debug('tb',tb_keys[tb_keys_ind], tb_sites[tb_keys[tb_keys_ind]])
            #logging.debug('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
        
        elif masks_ind >= len(masks_key) and lines_ind < len(lines):
            #after mask list is completely read
            logging.info('no more masks, lines only')

            line = lines[lines_ind]
            line_start = line[1]

            if prev != None:
                    #change here
                    overlap,change,newline = check_prev_line(prev, line)
                    if overlap == True and change != None:
                        #do i need to change this?
                        logging.debug('change!!!!!', change)
                        #all_lines[change[0]] = change[1]
                    if newline != None:
                        all_lines.append(newline)
                        logging.debug('append!!!!', newline)
            if prev == None or overlap == False:
                logging.debug('add as is!!!!!!', 'line', line, 'prev', prev )
                all_lines.append(line)
            lines_ind += 1


        elif lines_ind >= len(lines) and masks_ind < len(masks_key):
            #after lines list is completely read and masks are left 
            logging.info('no more lines, masks only')

            mask_start = masks_key[masks_ind]
            mask_end =  masks[mask_start]

            if prev != None:
                #change here
                overlap,change,newline = check_prev_line(prev, ['-', str(mask_start), str(mask_end-mask_start)])
                if overlap == True and change != None:
                    #do i need to change this?
                    logging.debug('change!!!!!', change)
                    #all_lines[change[0]] = change[1]
                if newline != None:
                    all_lines.append(newline)
                    logging.debug('append!!!!', newline)
            if prev == None or overlap == False:
                #all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])    
                logging.debug('add mask as is!!!!!!', 'mask', mask_start, mask_end, 'prev', prev )
                all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
                

            #all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
            masks_ind += 1 
        
        
        prev = all_lines[-1]
        logging.debug('prev', prev)

    #need to add final prev?  
    logging.debug('final prev', prev)
    return all_lines

def mask2ref(lines, tb_masks):
    logging.debug(tb_masks)
    final = [lines[0]]
    #editing thought: would likely benefit from being a list rather than a dictionary 
    tb_keys = sorted(tb_masks.keys())
    logging.debug('all lines', lines)
    tb_keys_ind = 0
    lines_ind = 1
    cont = 0
    #track previous line in all_sites
    prev = None
    #iterate through tb_keys and ld_keys exactly 1 time, track index of keys as you iterate
    while tb_keys_ind < len(tb_keys) or lines_ind < len(lines):
        #if still iterating through both lists
        if tb_keys_ind < len(tb_keys) and lines_ind < len(lines): 
            tb_start = tb_keys[tb_keys_ind]
            tb_end =  tb_masks[tb_keys[tb_keys_ind]]
            line_start = int(lines[lines_ind][1])
            line_end = line_start+int(lines[lines_ind][2])
            logging.debug('tb start', tb_start, 'tb_end', tb_end, 'line start', line_start, 'line end', line_end)

            if line_end <= tb_start:
                '''
                #line is completely to the left of the universal site
                #if prev != None:
                    #before adding to all_sites, make sure it doesn't overlap prev
                    #overlap, change = check_prev_mask(prev, [ld_start, ld_end])
                        #if overlap is detected AND requires prev to be updated
                        if overlap == True and change != None:
                            #update prev in all_sites
                            all_sites[change[0]] = change[1]
                    if prev == None or overlap == False:
                        #add new site to all_sites
                        all_sites[ld_start] = ld_end
                    #update ld index because site was processed
                    '''
                final.append(lines[lines_ind])
                lines_ind += 1
            
            elif line_start >= tb_end:
                logging.debug('dont add to final yet')
                tb_keys_ind += 1
                #final.append(lines[lines_ind])
            
            elif line_start >= tb_start and line_end <= tb_end:
                logging.debug('full overlap')
                lines_ind += 1
            
            elif line_start < tb_start and line_end > tb_end:
                logging.debug('full overlap OPPOSITE', 'line',lines[lines_ind], 'tb', tb_start, tb_end)
                #first segment of line before TB mask
                new_line = [lines[lines_ind][0], str(line_start), str(tb_start-line_start)]
                #update line after first segment AND TB mask
                lines[lines_ind][1] = str(tb_end)
                lines[lines_ind][2] = str(line_end-tb_end)
                logging.debug('new line', new_line)
                logging.debug('update line', lines[lines_ind])
                final.append(new_line)
                tb_keys_ind += 1

            elif line_start <= tb_start and line_end > tb_start and line_end <= tb_end:
                logging.debug('left overlap')
                new_line = [lines[lines_ind][0], str(line_start), str(tb_start-line_start)]
                logging.debug('old line', lines[lines_ind])
                logging.debug('new line', new_line)
                final.append(new_line)
                lines_ind += 1

            elif line_start >= tb_start and line_start < tb_end and line_end > tb_end:
                logging.debug('right overlap')
                new_line = [lines[lines_ind][0], str(tb_end), str(line_end-tb_end)]
                logging.debug('new line', new_line)
                final.append(new_line)
                tb_keys_ind += 1
                lines_ind += 1



            else:
                logging.debug('overlap?' )
        

        elif tb_keys_ind == len(tb_keys) and lines_ind < len(lines):
            logging.debug(lines[lines_ind])
            final.append(lines[lines_ind])
            lines_ind += 1

        #elif tb_keys_ind < len(tb_keys) and lines_ind == len(lines):
        #    logging.debug(lines[lines_ind])
        #    final.append(lines[lines_ind])
        #    lines_ind += 1

        

            """
                elif ld_end >= tb_start:
                    #this accounts for ld overlapping tb on the left
                    if prev != None:
                        overlap, change = check_prev_mask(prev, [ld_start, tb_end])
                        if overlap == True and change != None:
                            #if new region overlaps, update prev
                            all_sites[change[0]] = change[1]
                    #if there is no overlap, add new region to all_sites
                    if prev == None or overlap == False:
                        all_sites[ld_start] = tb_end
                    #since both regions are added at the same time, update index for both lists
                    tb_keys_ind += 1
                    ld_keys_ind += 1
            
            elif ld_start <= tb_end and ld_end > tb_end:
                #if ld overlaps tb on the right
                if prev != None:
                    #before adding to all_sites, make sure it doesn't overlap prev
                    overlap,change = check_prev_mask(prev, [tb_start, ld_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                #if there is no overlap, add new region to all_sites
                if prev == None or overlap == False:
                    all_sites[tb_start] = ld_end
                #since both regions are added at the same time, update index for both lists
                tb_keys_ind += 1
                ld_keys_ind += 1

            elif ld_start > tb_end:
                #if there is no overlap between ld and tb, and ld is on the right
                if prev != None:
                    #before adding to all_sites, make sure it doesn't overlap prev
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    #if there is no overlap, add new region to all_sites
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1

            elif ld_start >= tb_start and ld_end <= tb_end:
                # if tb and ld fully overlap with ld inside
                if prev != None:
                    #make sure tb doesnt overlap prev 
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != False:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            
            elif ld_start <= tb_start and ld_end >= tb_end:
                #if ld and tb fully overlap with tb inside
                if prev != None:
                    #make sure ld doesnt overlap prev
                    overlap,change = check_prev_mask(prev, [ld_start, ld_end])
                    if overlap == True and change != None:
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[ld_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1
                
            #covered all 6 possible positions of the two regions
            '''
            DEBUG
            logging.debug('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
            '''

        #possible bug: prev overlaps with one list the first time the other list expires
        elif tb_keys_ind >= len(tb_keys) and ld_keys_ind < len(ld_keys):
            #if reach end of tb_masks process ld only
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]
            all_sites[ld_start] = ld_end
            ld_keys_ind += 1

        #possible bug: prev overlaps with one list the first time the other list expires
        elif ld_keys_ind >= len(ld_keys) and tb_keys_ind < len(tb_keys):
            #if reach end of ld masks process tb only
            tb_start = tb_keys[tb_keys_ind]
            tb_end = tb_sites[tb_keys[tb_keys_ind]]
            all_sites[tb_start] = tb_end
            tb_keys_ind += 1 

        #keep all sites sorted
        if len(all_sites) > 0:
            all_sites_keys = sorted(all_sites.keys())
        #track last site for every iteration
        prev = [all_sites_keys[-1],all_sites[all_sites_keys[-1]]]
        """
        #cont += 1
        #if cont == 1000:
        #    break
        logging.debug('end', 'lines', len(lines), 'tb', len(tb_keys))
        logging.debug('tb', tb_keys_ind, 'lines', lines_ind)

    
    #return all_sites
    return final
    


#determine if low-coverage samples exceed 5% of genome length
#note: future iterations of this software may include missing sites in VCF but that is not currently included
#note: future iterations of this software may determine if universal mask sites overlap with low-coverage sites 
# but that is not currently included
def missing_check(lenref, ld):
    '''
    determine what percentage of genome length is low-coverage
    Args:
        lenref: len of reference genome 
        ld: dictionary of low-depth regions 
    Out:
        missing_count/lenref: a float value <1 indicating the percentage of genome 
    '''
    #NOT CURRENTLY NEEDED 
    #how many missing lines ended up in VCF
    #miss_total = 0
    #for line in missing:
    #    miss_total += len(line[3])

    #calculate in main part of script?
    #ld = mask_low_depth(cf,cd)

    #NOT CURRENTLY NECESSARY
    #how many universal mask regions are there
    #mask_count = 0
    #for t in tbmask:
    #    mask_count += int(tbmask[t])- int(t)

    #len of genome after universal masks are excluded
    #eff_lenref = lenref-mask_count


    #how many low depth mask regions are there
    missing_count = 0
    for l in ld:
        missing_count += int(ld[l])-int(l)


    #rules out samples that could never pass quality check no matter what 
    #if missing_count/lenref > cc:
    #    logging.error('fail')
    #    return False, missing_count/lenref

    #return all samples that couldn't fail quality check no matter what 
    #elif (missing_count + miss_total)/eff_lenref <= cc:
    #    return True, missing_count/lenref

    #FOR ALL BORDERLINE SAMPLES (IS THIS NECESSARY?)
    #logging.debug('pass')
    return missing_count/lenref
    

#SCRIPT STARTS HERE

#you can  assume sane people will end their gzipped files with .gz.
#and it erroring if they don't is perfectly acceptable. 
if vcf.endswith('.vcf'):
    binary = False
elif vcf.endswith('.vcf.gz'):
    binary = True
else:
    logging.error("Unsupported file extension for input- use '.vcf' or '.vcf.gz'.")
    logging.debug(f"File name read: {vcf}. Extension detected: {vcf.split('.')[-1]}")
    exit(1)

try:
    if binary:
        lenRow, samps = count_samples_bin(vcf)
    else:
        lenRow, samps = count_samples(vcf)
except IndexError:
    logging.error("Fail to read input vcf columns- check formatting.")
    sys.exit(1)

#be careful w dictionaries!!!
files = make_files(samps, wd)
try:
    if binary:
        read_VCF_bin(vcf, files)
    else:
        read_VCF(vcf, files)
except IndexError:
    logging.error("Fail to read input vcf data- check formatting.")
    sys.exit(1)

if tbmf != None:
    masks = mask_TB(tbmf)
else:
    #default to no masking.
    masks = {}
#this is not parallelized, the more samples in the vcf the longer this will take
#logging.debug(files)
#logging.debug(masks)

final_lines = []
for f in files:
    files[f].close()

    #note if a multisample VCF is submitted to this script, there is no way to mask low-depth
    #find low coverage regions for each sample 
    if cf != None:
        ld = mask_low_depth(cf,cd)
    else:
        ld = None

    sample = os.path.basename(files[f].name)[:-4]
    logging.info('Working on sample', sample)
    filepath = files[f].name
    if args.tag_filter:
        os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}")
        os.system(f"rm {filepath}")
        filepath = f"{filepath}.filt"
    if ld != None:
        #currently quality assessment requires a coverage file, if not coverage is provided the script will fail 
        error = missing_check(len_ref, ld)
    else:
        error = 'nomask'

    #if there is a provided coverage file it will be used to mask low coverage (less than cd) regions 
    #note that only one coverage file can be provided and it will result in an error if the vcf has more samples than coverage files 
    lines = vcf_to_diff(filepath)
    os.remove(filepath)

    all_lines = mask_and_write_diff(ld, masks,lines, samps)
    
    logging.info('MASK2REF')
    filtered_lines = mask2ref(all_lines, masks)
    if args.mask_report:
        with open(f'{wd}{sample}.report.txt','w') as o:
            o.write(f'{sample}.diff\t{error}\t{cd}\n')
    if args.split_diff:
        with open(f'{wd}{sample}.diff','w') as o:
            for line in filtered_lines:
                o.write('\t'.join(line)+'\n')
    final_lines.extend(filtered_lines)

with open(args.output_diff,'w+') as of:
    for line in final_lines:
        print('\t'.join(line),file=of)

logging.info("Finished")