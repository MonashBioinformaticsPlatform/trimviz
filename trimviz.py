#!/usr/bin/env python
#
# TODO: mapping rate vs trim amount; multimapping vs trim
# proper temp files plus cache? (currently its a kind of race condition)
# dependency testing
# length profile of trimmed
# Analysis of soft-clipping from bam
# paired-end (+ paired-end overhang -> suggested clipping)
# test if 'raw' fastq is actually uneven lengths (pretrimmed from seq facility)


__author__ = 'stuartarcher'

import getopt, subprocess, random, re, sys, os, tempfile, time, gzip, pipes, pysam

# An iterator for fastq files.  Takes a open file object.  Returns a dict
class fastq:
    def __init__(self, f):
        self.f = f
    def __iter__(self): return self
    def next(self):
        id = re.sub('\/.+$', '', self.f.next().rstrip()) # remove all after forward slash
        s  = self.f.next().rstrip().upper() # convert to uppercase (<----- maybe change in future? Cutadapt can convert to lowercase instead of trimming, this will be missed)
        self.f.next()
        q = self.f.next().rstrip()
        #return {'id': id.partition(' ')[0], 'seq': s, 'qual': q}
        return ([id.partition(' ')[0], s, q])
    
class dummy_fastq:
    def __init__(self, f):
        self.f = f
    def __iter__(self): return self
    def next(self):
        return ''

    
def main():
    if len(sys.argv) < 2 or not sys.argv[1] in {'FQ','SC'}:
        print('Error: please specify either fastq-fastq mode (FQ) or soft-clipping analysis mode (SC)')
        print_help()
        exit()
    if sys.argv[1] == 'FQ':
        softClipping = False
    else:
        softClipping = True
        
    options, remainder = getopt.getopt(sys.argv[2:], 'u:t:U:T:o:O:b:g:c:n:v:w:a:A:f:g:r:s:Rzeh', ['untrimmed=',
                                                                                   'trimmed=',
                                                                                   'untrimmed_r2=',
                                                                                   'trimmed_r2=',
                                                                                   'out_dir=',
                                                                                   'out_dir_all=',
                                                                                   'bam=',
                                                                                   'genome_fasta=',
                                                                                   'classes=',
                                                                                   'sample_size=', 
                                                                                   'indiv_reads=',
                                                                                   'heatmap_reads=',
                                                                                   'adapt=',
                                                                                   'adapt_file=',
                                                                                   'agg_flank='
                                                                                   'rid_file='
                                                                                   'seed=',
                                                                                   'representative',
                                                                                   'gzipped',
                                                                                   'read2only',
                                                                                   'help'])

    gr_FN = 'indivReadPlots.pdf'
    aggGr_FN = 'aggPlot.pdf'
    tblOut_FN = 'trimVisData.tsv'
    class_opts = ['uncut','5pcut','3pcut','removed','indel','generated_warning']  # not all reads generating warning are included (many are not plottable). indel just for debugging purposes
    
    # set defaults
    orig_FN1 = str('')  # u   'untrimmed='
    proc_FN1 = str('')  # t   'trimmed='
    orig_FN2 = str('')  # U   'untrimmed_r2='
    proc_FN2 = str('')  # T   'trimmed_r2='
    out_DN = str('')    # o|O 'out_dir='
    keepTmp = False     # O                # <------- TODO make use of this parameter
    bam_FN = str('')    # b   'bam='
    gfasta_FN = str('') # g   'genome_fasta='
    coi_raw = ['all']   # c   'classes='
    target_n_pre=50000  # n   'sample_size='
    nvis=20             # v   'indiv_reads='  # if -R is not set, will aim to visualize this many reads from each individual trimming-class
    maxAggN = 200       # w   'heatmap_reads='
    aggFlank=20         # f   'agg_flank=
    madapt=list()       # a   'adapt='...(madapt[0]) # <------- TODO make use of additional adapters (1 only currently)
    adapt_FN = str('')  # A   'adapt_file='  # <------- TODO make use of this parameter
    rid_FN = str('')    # r   'rid_file='
    rseed = 1           # s   'seed='
    balance = True      # R   'representative'
    gzipped = True      # z   'gzipped'
    read2only = False   # e   'read2only'   # <------- TODO paired-end mode
    help = False        # h   'help'
    
    for opt, arg in options:
        if opt in ('-u', '--untrimmed'):
            orig_FN1 = arg
        elif opt in ('-t', '--trimmed'):
            proc_FN1 = arg
        elif opt in ('-U', '--untrimmed_r2'):
            orig_FN2 = arg
        elif opt in ('-T', '--trimmed_r2'):
            proc_FN2 = arg
        elif opt in ('-o', '--out_dir'):
            out_DN = arg
        elif opt in ('-O', '--out_dir_all'):
            out_DN = arg
            keepTmp = True
        elif opt in ('-b', '--bam'):
            bam_FN = arg
        elif opt in ('-g', '--genome_fasta'):
            gfasta_FN = arg # g
        elif opt in ('-c', '--classes'):
            coi_raw = arg.split(',')
        elif opt in ('-n', 'sample_size'): # should be large enough to be able to get a balanced set of read-trim classes
            target_n_pre = int(arg)
        elif opt in ('-v', '--indiv_reads'):  # number of individual reads to visualize 1-by-1
            nvis = int(arg)
        elif opt in ('-w', '--heatmap_reads'): # number of reads to visualize in heatmaps
            maxAggN = int(arg)
        elif opt in ('-f',  '--agg_flank'):
            aggFlank = int(arg)    
        elif opt in ('-a', '--adapt'):
            madapt = arg.split(',')
        elif opt in ('-A', '--adapt_file'):
            adapt_FN = arg     # 
        elif opt in ('-r', '--read_id_file'):
            rid_FN = arg    # t
        elif opt in ('-s', '--seed'):
            rseed = int(arg)
        elif opt in ('-R', '--representative'):
            balance = False
        elif opt in ('-z', '--gzipped'):
            gzipped = True
        elif opt in ('-w', '--read2only'):
            read2only=True
        elif opt in ('-h', '--help'):
            print_help()
            
    # Open a file.  Use gzip based on filename, or 'gzipped' flag
    def gzopen(fname):
        if gzipped or fname.endswith('.gz'):
            return gzip.open(fname, 'rb')
        else:
            return open(fname, 'r')
    
    ###########################################################################################
    #  MAIN part 1:  arguments logic
    ###########################################################################################
    coi = [x for x in coi_raw if x in class_opts]
    if len(coi) < len (coi_raw) and coi_raw != ['all']:
        print ("Warning, values given in -c should include only %s, comma-separated; or just -c 'all'") % (', '.join(class_opts))
    if not cmd_exists('seqtk'):
        print(' *** Could not find seqtk executable. Please install seqtk. Exiting. ***')
        exit()
    if out_DN == '':
        print('Error: output directory required.')
        print_help()
        exit()
    if os.path.exists(out_DN):
        print('Error: output directory already exists. Exiting.')
        exit()
    try:
        os.makedirs(out_DN)
    except:
        print('Error: Could not create output directory %s. Exiting.') % (out_DN)
        exit()   
    
    if not os.path.isfile(orig_FN1):
        if len(orig_FN1) > 0:
            print('Could not find input fastq file: ' + orig_FN1 + '. Exiting.')
        else:
            print orig_FN1
            print ('No untrimmed R1 file given (required). Exiting.')
            print_help()
        exit()
    if softClipping:
        if not os.path.isfile(bam_FN):
            if len(bam_FN) > 0:
                print('Could not find bam file ' + bam_FN + '. Exiting.')
            else:
                print('No bam filename was given (required in SC mode). Exiting.')
                print_help()
            exit()
        if not (proc_FN1=='' and proc_FN2==''):
            print ("Warning, -t and -T options not used in SC mode. Ignoring.")  
    else:   
        if not os.path.isfile(proc_FN1):
            if len(proc_FN1) > 0:
                print('Could not find trimmed R1 fastq file ' + proc_FN1 + '. Exiting.')
            else:
                print ('No trimmed R1 file given, but this is required in FQ mode. Exiting.')
                print_help()
            exit()

    if os.path.isfile(bam_FN) and not os.path.isfile(gfasta_FN):
        if gfasta_FN == '':
            print ('Bam file is given but could not corresponding genome fasta. Exiting.')
        else:
            print ('Bam file is given but could not find corresponding genome fasta file '+ gfasta_FN + '. Exiting.')
        print_help()
        exit()
        
    if not os.path.isfile(adapt_FN) and len(adapt_FN) > 0:
        print('Could not find file: ' + adaptfile + '. Exiting.')
        print_help()
        exit()
        
    random.seed(rseed)

    if nvis > 500:
        nvis = 500
        sys.stderr.write('Warning: Nominated sample size for individual read rendering is larger than 500. Lowering.\n')
    
    paired=False 
    if orig_FN2 != '' or proc_FN2 != '':
        if os.path.isfile(orig_FN2) and (os.path.isfile(proc_FN2) or softClipping):
            paired=True
        else:
            sys.stderr.write('Error: A Read-2 filename was given, but could not locate both Read-2 files (or a bam file).\n')
            quit()
            # <---------- TODO: paried-end mode

    if len(adapt_FN)>0:
        if os.path.exists(adapt_FN) and os.path.isfile(adapt_FN):
            with open (adapt_FN, 'r') as fin:
                for line in fin:
                    madapt.append(line.rstrip())
        else:
            print 'Warning. File: ', adapt_FN, ' not found. Not using adapter file. '  # TODO implement this
    
    
    ###########################################################################################
    #  MAIN part 2:     sample from pre- processed fastq file
    ###########################################################################################     
    if not os.path.exists(out_DN +'/trimVisTmpFiles'):
        os.makedirs(out_DN +'/trimVisTmpFiles') 
    tmpPre1_FN = out_DN + '/trimVisTmpFiles/tmpPre1.fq.gz'
    tmpPost1_FN = out_DN + '/trimVisTmpFiles/tmpPost1.fq.gz'
    tmpPre2_FN = out_DN + '/trimVisTmpFiles/tmpPre2.fq.gz'
    tmpPost2_FN = out_DN + '/trimVisTmpFiles/tmpPost2.fq.gz'
    readIDs1_FN = out_DN + '/trimVisTmpFiles/readIDsPre.lst'
    
    # sample original fq:
    if rid_FN == '':
        print 'target number of reads (orig fastq): %d' % (target_n_pre)
        pcnt_sgn = """%"""
        # split output into a temp fq file (tmpPre1_FN) and a readname list file (readIDs1_FN)
        # remove @, anything after a space in readname; append a tab char to improve fgrep specificity for bam-searching
        # note: sed 's/[ \/].*$//' messes up seqtk search of fastqs with RNs ending in /1 or /2 (but is required to search bams)
        cmd3 = "seqtk sample -s%d %s %d | tee >(awk '1 == NR %s 4' | sed 's/@//'  | sed 's/[ ].*$//' | sed 's/$/\t/' > %s ) | gzip > %s" % (rseed,
                                                                                                                          pipes.quote(orig_FN1),
                                                                                                                          target_n_pre,
                                                                                                                          pcnt_sgn,
                                                                                                                          readIDs1_FN,
                                                                                                                          pipes.quote(tmpPre1_FN))       
    else: # user-specified RIDs
        cmd2 = "cat %s | sed 's/@//'  | sed 's/[ ].*$//' | sed 's/$/\t/' > %s" % (rid_FN, readIDs1_FN) # clean up user-specified rids for seqtk
        with open(out_DN + '/trimVisTmpFiles/tmp_bash2.sh', 'w') as fout:
            print >> fout, cmd2
        sout2 = subprocess.check_output(['bash', out_DN + '/trimVisTmpFiles/tmp_bash2.sh'])
        print (sout2)
        cmd3 = "seqtk subseq %s %s | gzip > %s " % (pipes.quote(orig_FN1), readIDs1_FN, tmpPre1_FN)
    
    with open(out_DN + '/trimVisTmpFiles/tmp_bash3.sh', 'w') as fout:
            print >> fout, cmd3
            
    # retrieve from ORIG FQ using seqtk
    print ('sampling from original fastq using seqtk...')
    sout3 = subprocess.check_output(['bash', out_DN + '/trimVisTmpFiles/tmp_bash3.sh'])
    # this either runs seqtk sample -> readIDs1_FN & tmpPre1_FN; or, if rid_file (-r) is already specified, just seqtk subseq -> tmpPre1_FN
    # either way, files readIDs1_FN & tmpPre1_FN now exist
    print (sout3)
    
    # get rIDs in dict for python and load results into dict 'pre': pre[rID] = [seq, qual]
    pre = dict()
    f_pre = gzopen(tmpPre1_FN)
    for l in fastq(f_pre):
        pre[l[0]] = [l[1], l[2]]                   
    f_pre.close()
    
    both = dict()                                  # this will be used to store pre- and post- trimmed fqs
    rid_class={'uncut':[], '5pcut':[],'3pcut':[], 'removed':[], 'generated_warning':[], 'indel':[]} # will also classify them into the 3 classes and make ID list for each class (-> rid_class)
        
    ####################################################################################
    # MAIN part 2.1: (optional): extract from bam & fasta; save as temp file & dict (fastaD)
    ####################################################################################
 
    if (bam_FN != ''):
        # remove /1 and /2 from end of readnames for bam-searching, replace the terminal tab char
        readIDs1_FN2 = readIDs1_FN + '.2'
        cmd3B = "cat %s | sed 's/[\/].*$/\t/' > %s" % (readIDs1_FN, readIDs1_FN2)
        sub_bamFN = out_DN + '/trimVisTmpFiles/bamEntries1.sam'
        cmd4A = 'samtools view -H %s > %s' % (bam_FN,  sub_bamFN)
        cmd4B = 'samtools view %s | fgrep -f %s | awk \'NF > 8{print} 1\' >> %s' % (bam_FN, readIDs1_FN2, sub_bamFN)  # |  cut -f1-6,9
        with open(out_DN + '/trimVisTmpFiles/tmp_bash6.sh', 'w') as fout:
            print >> fout, cmd3B
            print >> fout, 'sleep 1'
            print >> fout, cmd4A
            print >> fout, 'sleep 1'
            print >> fout, cmd4B
        print ('searching bam file using samtools and fgrep...')
        time.sleep(1)
        sout = subprocess.check_output(['bash', out_DN + '/trimVisTmpFiles/tmp_bash6.sh'])
        print ('getting fasta file...')
        # shove fasta file into giant dictionary
        fastaD=dict()
        currChr=''
        currSeq=''
        with open(gfasta_FN, 'r') as fastaIn:
            for line in fastaIn:
                if line.startswith('>'):
                    if currChr != '':
                        fastaD[currChr]=currSeq
                    currChr=line.rstrip()[1:]
                    currSeq=''
                else:
                    currSeq=currSeq+line.rstrip()
            fastaD[currChr]=currSeq
    
    
    ####################################################################################
    # MAIN part 2.2: (FQ mode only): extract from post-trimmed fastq; merge into dict 'both'
    ####################################################################################  
    
    if not softClipping:     
        # retrieve same reads (in readIDs1_FN) from TRIMMED FQ (proc_FN1) using seqtk
        cmd5 = "seqtk subseq %s %s | gzip > %s " % (proc_FN1, readIDs1_FN, tmpPost1_FN)
        with open(out_DN + '/trimVisTmpFiles/tmp_bash4.sh', 'w') as fout:
            print >> fout, cmd5
        print ('Finding the reads in trimmed fastq file using seqtk...')
        sout4 = subprocess.check_output(['bash', out_DN + '/trimVisTmpFiles/tmp_bash4.sh'])
        print (sout4)
        
        # load results into dict 'post'
        post = dict()
        f_post = gzopen(tmpPost1_FN)
        for l in fastq(f_post):
            post[l[0]] = [l[1], l[2]]  # post[id] = [seq, qual]
        f_post.close()
        
        # find these reads in pre-trimmed fastq:
        # both <- dictionary of lists (both[rID] = [ 0-postSeq 1-postQual 2-preSeq 3-preQual ])
        for id, rdat in pre.items():
            if id in post:
                t = post.pop(id)     # Remove it from post
                both[id] = t + rdat  # both[id] = 0-postSeq 1-postQual 2-preSeq 3-preQual
            else:
                both[id] = ['', ''] + rdat
                #rid_class['removed'].extend([id])
        # check that post is now emptied out
        for r in post:
            print "WARNING: read ---- " + r + "  ---- not found in pre-trimmed data."
            #rid_class['generated_warning'].extend([r]) # actually can't plot this (no 'pre-trimmed' entry)!
            
            
        #########################################################
        #  Guess 5' and 3' trim sites using both seq and qual strings
        #  add as elements [4] and [5] in 'both'
        #########################################################
    
        removeread = list()
        maxReadLen=0
        for r, rdat in both.items(): # rdat: 0-postSeq 1-postQual 2-preSeq 3-preQual
            if len(rdat[0]) > 0: # read present in processed
                if len(rdat) == 4 and rdat[0] in rdat[2]:
                    maxReadLen = max(maxReadLen, len(rdat[2]))
                    seqpos = mfind(rdat[2], rdat[0])  # ref(pre), query(post) mfind returns 0 if no 5' trim
                    qualpos = mfind(rdat[3], rdat[1])
                    fps = [i for i in seqpos if i in qualpos]
                    #if len(fps)==1:  # almost certainly only 1 unique start site for the match (using both Seq and Qual) if length reasonable
                    both[r].append(fps[0] + 1)                # [ both[r][4]=start rel pre    ...so 3' trim len = (len([2])- [5]
                    both[r].append(fps[0] + 1 + len(rdat[0])) # for len=1, start, end on n,n+1 nucleotides
                    if len(fps)> 1:
                        print "Warning - ambiguous trimming (multiple possible trim-points) encountered for read %s" % (r)
                        print len(fps)
                        print rdat
                        rid_class['generated_warning'].extend([r])
                        #both.pop(r)
                else:
                    print "WARNING: read ", r, " does not have all associated data, or the trimmed sequence was not a substring of full sequence."
                    removeread.append(r)
            else:
                both[r].extend([1,1]) # add dummy values for 5'cut and 3' cut site
            for adapt in madapt:
                both[r].append(seqmatch(rdat[2], adapt))  # ref=l[1], query=adapt
                
        for rr in removeread:
            both.pop(rr) # remove incomplete entries

    ###############################################################################################################################################
    # MAIN part 3: Extract from sub-bam file -> baminfo dict. If not soft clipping, use pre-existing 'both' to add trim-lengths to clipped bases at ends
    #    (for NOT soft clipped mode: assumes bam is post-trimmed alignment, not pre-trimmed alignment)
    ###############################################################################################################################################
    
    if (bam_FN != ''):
        bamInfo = dict()
        seqsNotFound=list()
        nmapped={k:0 for k in pre.keys()}
        trimLead = 0   # will remain zero if soft-clipping
        trimTrail = 0  # will remain zero if soft-clipping
        
        bfin = pysam.AlignmentFile(sub_bamFN,'r')
        genome_fa = pysam.Fastafile(gfasta_FN)  #'/references/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa'
        BamAlnsDone=set()
        internal_insertions = set()
        for r in bfin:
            
            # 1) build refBlocks into coords for fethching the ref sequence (non-contiguous) that look like this:
            #   [ (5'trimmed seg*), (5'soft-clipped seg **), (query align seg1), ... (query align seg n), (3' soft-clipped seg**), (3' trimmed seg*) ]
            # * only if softClipping == False (i.e. trimmed file is fastq not bam)  ** only if S in cigar.
            # (trimmed seg) - (SC seg) - (align seg 1) should all form a contiguous block of sequence
            #
            # 2) collect 'trim' coords into both[id2][4] (5' trim length/coord) and both[id2][5] (3' trim COORD in contiguous read bases, rel to start of untrimmed sequence)
            # for fastq-fastq comparison these are 5'trim + 5'softclip; and tot untrimmed len - (3'trim + 3'softclip)
            # for fastq-bam comparison (softClipping==True) these are just 5'SC; and tot read len - 3' SC
            
            id2 = '@'+r.query_name
            if read2only and r.is_read1:
                continue
            if not read2only and r.is_read2:
                continue
            if r.is_secondary: #???
                continue           
            if r.is_unmapped:
                if softClipping:
                    both[id2] = ['', ''] + pre[id2]
                    #rid_class['removed'].extend([id2])
                    both[id2].extend([1,1])
                    for adapt in madapt:
                        both[id2].append(seqmatch(pre[id2][0], adapt))  # ref=pretrim seq, query=adapt
                continue
            if id2 in BamAlnsDone:
                #print 'Warning: read %s has repeated primary alignments. Skipping subsequent one.' % id2
                continue
            BamAlnsDone.add(id2)
            refname = r.reference_name
            if not refname in genome_fa.references:
                seqsNotFound.append(refname)
                continue
            reflen = genome_fa.get_reference_length(refname)
            rvcmp=r.is_reverse
            #if rvcmp:
            #    continue # <-----------------------------------
            cigTuples = r.cigartuples
            refBlocks = r.get_blocks()
            # get_reference_sequence_output = r.get_reference_sequence()
            scLead = 0
            scTrail = 0
            nlpad = 0
            nrpad = 0
            cigTuples2 = [x for x in cigTuples] # copy list
            if cigTuples[0][0] == 4:  # leading soft-clipped. TODO: hard clipped is '5'
                scLead += cigTuples[0][1]
                nlpad = max(0, scLead - refBlocks[0][0]) # if runs off edge of contig:
                scLead -= nlpad
                refBlocks = [(refBlocks[0][0]-scLead, refBlocks[0][0])] +  refBlocks
            if cigTuples[-1][0]==4:  # trailing soft-clipped.
                scTrail += cigTuples[-1][1]
                nrpad = max(0, (refBlocks[-1][1] + scTrail)-reflen) # if runs off edge of contig:
                scTrail -= nrpad
                refBlocks =  refBlocks + [(refBlocks[-1][1], refBlocks[-1][1] + scTrail)]
            if not softClipping:
                # in this case, to get gSeq, need to account for both soft-clipping AND the read trimming that occurred prior to mapping
                if rvcmp:
                    trimTrail = both[id2][4]-1
                    trimLead = len(both[id2][2]) + 1 - both[id2][5]  #len PreSeq minus 3' clipping point  
                else:    
                    trimLead = both[id2][4]-1
                    trimTrail = len(both[id2][2]) + 1 - both[id2][5]  #len PreSeq minus 3' clipping point  
                # if rvcmp:
                #     lhs = refBlocks[0][0] - trimTrail
                #     rhs = refBlocks[-1][1] + trimLead
                # else:
                lhs = refBlocks[0][0] - trimLead
                rhs = refBlocks[-1][1] + trimTrail
                # if non-trimmed runs off edge of contig:
                nlpad += abs(min(0,lhs))
                lhs = max( 0, lhs )
                nrpad += max( 0, rhs - reflen )
                rhs = min(rhs, reflen)
                if trimLead > 0:
                    cigTuples2 =  [(4,refBlocks[0][0]-lhs)] + cigTuples2 # add a leading portion of soft-clipped (we build gseg from gblocks+cigtuples2 so they must be consistent)
                    refBlocks =  [(lhs, refBlocks[0][0])] + refBlocks
                if trimTrail > 0:
                    cigTuples2 =  cigTuples2 + [(4,rhs-refBlocks[-1][1])] # add a trailing portion of soft-clipped (we build gseg from gblocks+cigtuples2 so they must be consistent)
                    refBlocks = refBlocks + [(refBlocks[-1][1], rhs)] # add trimmed bits to outer segments of refBlocks
            lpad = 'X' * nlpad
            rpad = 'X' * nrpad                
            gblocks = [genome_fa.fetch(refname, rb[0], rb[1]) for rb in refBlocks]

            # compensate for any insertions relative to the reference, by inserting X in reference segment:
            gblocks2=list()
            gblocks_pos=0
            for cti in range(len(cigTuples2)):
                if cigTuples2[cti][0] == 1: # "I" / insertion relative to the reference: add X's
                    gblocks2.append('X'* cigTuples2[cti][1])
                    rid_class['indel'].extend([id2])
                elif cigTuples2[cti][0] == 0 or cigTuples2[cti][0] == 4: # "M" / matching segment or "S" soft-clipped: add the gblock
                    gblocks2.append(gblocks[gblocks_pos])
                    if not (len(gblocks2[-1])) == cigTuples2[cti][1]:
                        print ('Warning: genomic block length doesnt match cigar operation (readname: %s )') % ( id2 )
                        rid_class['generated_warning'].extend([ id2 ])
                    gblocks_pos +=1

            gSeg = lpad + ''.join(gblocks2) + rpad
            if r.is_reverse:
                gSeg = srevcomp(gSeg)  # revcomp the alignment, because in viz everything is relative to the read orientation
                joinLocs=[x[0] for x in refBlocks[::-1]]
                refpos=nrpad
                cigTuplesTmp = cigTuples[::-1]
            else:
                joinLocs=[x[1] for x in refBlocks]
                refpos=nlpad
                cigTuplesTmp=cigTuples
                
            bamInfo[id2] = [gSeg, joinLocs] # don't really use joinLocs for anything anymore...
            
            if softClipping: # also fill out 'both' entry for alignment / softclipped portions
                rdat=pre[id2]
                bseq =  r.query_alignment_sequence
                #bqual = r.query_alignment_qualities
                bqual = '.'*len(bseq) # dummy string  <--------------TODO use real quals?
                fps_relRead = scLead
                if rvcmp:
                    bseq = srevcomp(bseq)
                    bqual = bqual[::-1]
                    fps_relRead = scTrail
                newBothEntry = [bseq, bqual] + rdat  # both[id] = 0-postSeq 1-postQual 2-preSeq 3-preQual
                newBothEntry.append(fps_relRead + 1)                # [ both[r][4]=start rel pre    ...so 3' trim len = (len([2])- [5]
                newBothEntry.append(fps_relRead + 1 + len(bseq))    # for len=1, start, end on n,n+1 nucleotides
                # if len(bseq) == len(rdat[0]):
                #     rid_class['uncut'].extend([id2])
                # else:
                #     rid_class['cut'].extend([id2])
                for adapt in madapt:
                    newBothEntry.append(seqmatch(rdat[0], adapt))  # ref=pretrim seq, query=adapt
                both[id2] = newBothEntry
                
        if len(seqsNotFound) > 0:
            print "Warning;  reference chromosome/contig names in bam file were not in fasta file. These are: " + ', '.join([ x for x in set(seqsNotFound)]) 
       
    ##########################################################################
    #   MAIN   Part 4 classify reads into 5pcut, 3pcut, uncut and removed
    ##########################################################################
    for r, rdat in both.items():
        if len(rdat[0]) == 0:
            rid_class['removed'].extend([r])
        else:
            if rdat[4] > 1: #5' trim (1-based)
                rid_class['5pcut'].extend([r])
            if len(rdat[2]) - len(rdat[0]) > rdat[4]-1: # if the length diff is more than what can be accounted for by 5' trim (1-based)
                rid_class['3pcut'].extend([r])
            if len(rdat[2]) == len(rdat[0]):
                rid_class['uncut'].extend([r])
        
        
    ##################################################################################
    #   MAIN    part 5: prepare long-form output for 1-by-1 read vis R
    ##################################################################################            

   # sub-select balanced proportions of uncut, cut, removed for plotting
    lookupCls = dict()
    if (balance and coi_raw != ['all']):
        toPlot=[]
        for cls in coi:
            ids = rid_class[cls]
            for r in ids:
                #if r in lookupCls:
                #    print "warning: read %s is in multiple classes: %s and %s " % (r, lookupCls[r], cls)
                lookupCls[r]=cls
            rs = random.sample(ids, k=min( len(ids) , nvis ))
            toPlot = toPlot + rs
    else:
        toPlot = random.sample(both.keys(), nvis)
    print ("Trim-classification results from the random sample of %d reads") % (target_n_pre)
    trimClassTbl = dict()
    for cls in rid_class.keys():
        ncls = len(rid_class[cls])
        if balance and ncls < nvis:
            print ("There were %d reads classed as %s. (Warning: this was less than requested for plotting (-v). Increase sample size using -n)") % (ncls , cls)
        else:
            print ("There were %d reads classed as %s.") % (ncls , cls)
        trimClassTbl[cls] = ncls

    
    #with tempfile.NamedTemporaryFile(delete = False) as tempf:
    tempfname = out_DN + '/trimVisTmpFiles/trimviz_readData.tsv'
    with open (tempfname,'w') as fout:
        line = ['read', 'position', 'seq', 'qual', 'fp_cutoff', 'tp_cutoff']
        line.extend(madapt)
        if not bam_FN == '':
            line.append('genomic_seq')
        print >> fout, '\t'.join(line)
        j=0
        for r in toPlot:
            j+=1
            rdat=both[r]
            for i in range(0, len(rdat[2])):
                line = [ str(x) for x in [   r, (i+1), rdat[2][i], ord(rdat[3][i])-34, rdat[4], rdat[5]   ]  ]
                for ad in range(0, len(madapt)):
                    if 6+ad >= len(rdat):
                        print "it is rdat!"
                        print ad
                        print r
                        print j
                        print rdat
                    if i >= len (rdat[6+ad]):
                        print "it is string index!"
                    line.append(str(ord(rdat[6+ad][i])-34))  # convert ascii+33 -> zero-based num
                if bam_FN != '':
                    if r in bamInfo:
                        if len(bamInfo[r][0]) > i:
                            line.extend(bamInfo[r][0][i]) # add genome
                        else:
                            line.extend(['N'])
                    else:
                        line.extend(['N'])
                print >> fout, '\t'.join(line)
    #print "FILE1 READY:"
    #print tempfname
  
    ##################################################################################
    #   MAIN    part 6: prepare aggregate / trim-anchored read matrices
    ##################################################################################            
      
    #     runsheet={'3pcut':out_DN + '/trimVisTmpFiles/seq3psites.txt', '5pcut':out_DN + '/trimVisTmpFiles/seq5psites.txt', 'indel':out_DN + '/trimVisTmpFiles/indel_sites.txt'}
    runsheet={'3pcut':out_DN + '/trimVisTmpFiles/seq3psites.txt'} #, '5pcut':out_DN + '/trimVisTmpFiles/seq5psites.txt'}
    difflens = 0
    for cls, clsfile in runsheet.items():
        with open(clsfile, 'w') as fout:
            bs = range(aggFlank*2+1)
            if bam_FN == '':
                print >> fout, '\t'.join(['readID', 'tpCutPos'] + ['s'+str(i+1) for i in bs] + ['q'+str(i+1) for i in bs])
            else:
                print >> fout, '\t'.join(['readID', 'tpCutPos'] + ['s'+str(i+1) for i in bs] + ['q'+str(i+1) for i in bs] + ['g'+str(i+1) for i in bs])
            for r in rid_class[cls]:
                    rdat=both[r]                                # both[id] = 0-postSeq 1-postQual 2-preSeq 3-preQual ("post" is the "-t" fastq file; "pre" is the "-u" fastq file)
                    ps=padstr(rdat[2], rdat[5]-1, aggFlank, r)  # pre-seq; 3' trim site (0-based so should not exceed len(rdat[2]) ); flanking seq = 20. returns list
                    pq=padstr(rdat[3], rdat[5]-1, aggFlank, r)  # pad quals with 'N's which are not legit qual characters (highest is 'J')
                    if bam_FN != '':                      
                        if r in bamInfo:
                            gSeg=bamInfo[r][0]
                        else:
                            gSeg='N'*len(rdat[2])
                        gs=padstr(gSeg, rdat[5]-1, aggFlank, r) 
                        if softClipping:
                            if not len (gSeg) == len (rdat[2]): #if SC, bam align len should match -u file
                                difflens += 1
                                rid_class['generated_warning'].extend([r])
                        elif not len (gSeg) == len (rdat[2]):    # gSeg has already been adjusted in length to expand and match the -u read, if there was trimming between -u and -t
                            difflens += 1
                        print >> fout, '\t'.join([r, str(rdat[5])] + ps + [str(ord(x)) for x in pq] + gs)
                    else:
                        print >> fout, '\t'.join([r, str(rdat[5])] + ps + [str(ord(x)) for x in pq])  # N's -> 78. Highest legit q-val is 'J' (-> 74)
    
    if difflens > 0:
        if softClipping:
            print ( " ---- " )
            print ( "Warning: %d reads showed read length not equal between -u fastq and bam file. Was there an intervening read-trimming step? If not, it is not just soft-clipping being visualized." ) % difflens
            print ( "This may have unexpected effects. In SC mode it is strongly recommended to use only the fastq file that was directly input to the aligner as the -u argument." )
            print ( " ---- " )
        else:
            print ( " ---- " )
            print ("Warning: %d reads showed read length not equal between -u fastq and bam file, and not accounted for by trimming between -u and -t. Was there an intervening read-trimming step?" ) % difflens
            print ( "This may have unexpected effects. In FQ mode it is strongly recommended to use only the fastq file that was directly input to the aligner as the -t argument.")
            print ( " ---- " )
 
        
    ###########################
    # << CALL R SCRIPT HERE >>#
    ###########################
    cmd6 = 'Rscript ' +'./graph_ts.R '+ out_DN + ' ' + str(maxAggN)  #os.curdir

    #cmd6 = 'Rscript ' +'./graph_ts.R '+ tempfname + ' ' + out_DN + '/' + gr_FN + ' ' + out_DN +'/' + aggGr_FN + ' ' + str(maxAggN)  #os.curdir
    print (cmd6)
    rout = subprocess.check_output(cmd6, shell=True)
    
    print "R stdout:", rout
    
    ###########################
    # << MAKE HTML REPORT   >>#
    ###########################
    cmd7 = 'python ./' #os.curdir

    #cmd6 = 'Rscript ' +'./graph_ts.R '+ tempfname + ' ' + out_DN + '/' + gr_FN + ' ' + out_DN +'/' + aggGr_FN + ' ' + str(maxAggN)  #os.curdir
    print (cmd6)
    rout = subprocess.check_output(cmd6, shell=True)
    
    mode = "Fastq-fastq comparison mode, with alignment"
    if softClipping:
        mode = "Fastq-bam soft-clipping mode"
    elif bam_FN == '':
        mode = "Fastq-fastq comparison mode"
    report = makeReport(mode, out_DN, trimClassTbl, len(pre), orig_FN1, proc_FN1, orig_FN2, proc_FN2, bam_FN, gfasta_FN,
                        PDF_FNs=['TVheatmap_S.pdf','TVheatmap_Q.pdf','indiv_reads.pdf'], TPCs=['3pcut'])
    with open (out_DN+'/trimvis_report.html', 'w') as fout:
        print >> fout, report
   
# <<<<<<<<<<<<<<<<< END MAIN >>>>>>>>>>>>>>>>>>>

def cmd_exists(cmd):
    return any(
        os.access(os.path.join(path, cmd), os.X_OK) 
        for path in os.environ["PATH"].split(os.pathsep)
    )

def mfind (text, query):
    i=0
    idxs=list()
    while i < len(text):
        i = text.find(query, i)
        if i == -1:
            break
        idxs.append(i)
        i += 1 #finds overlapping too
    return(idxs)

def padstr(s, offset, flnk, r):
    if offset > len (s):
        print "warning, offset exceeds len s. ^%s^ " % (r)
        print str(offset) + "   " + s
    leftpad=max(0, flnk-offset)
    rightpad=max(0, offset + flnk - len(s))
    leftstr= 'N'*leftpad + s[max(0,offset-flnk):offset]
    rightstr=s[offset:(offset+flnk-rightpad)] + 'N'*rightpad
    totstr=leftstr + 'X' + rightstr
    return([x for x in totstr])
 
def srevcomp (seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X', '|':'|',
                  'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n', 'x':'x'}
    return"".join(complement.get(base, base) for base in reversed(seq))

def seqmatch (ref, query):
    aln_code = [0] * len(ref) # prepare list of n zeros
    #home-made adapter matcher
    for qp in range(len(query)-4):
        queryss = query[qp:qp+5]
        matches = [m.start() for m in re.finditer('(?='+queryss+')', ref)]
        for m in matches:
            m_extend = 5
            while ref[m : m + m_extend + 1] == query[qp : qp + m_extend + 1] and (m + m_extend) < len(ref) and (qp + m_extend) < len(query) :
                m_extend += 1  #extend match to the right, 5->1 6->2 ...
            m_extend = min(126, m_extend)   # limit to 128 to encode as ascii
            for nucpos in range(m, m + m_extend):
                # assign each position a score of m_extend (if m_extend is larger than pre-existing score)
                aln_code[nucpos] = max((aln_code[nucpos], m_extend)) 
    asascii33 = [chr(x+34) for x in aln_code]
    return( "".join(asascii33))

def flagDecoder (flag):
    bitlist=['template having multiple segments in sequencing',
             'each segment properly aligned according to the aligner',
             'segment unmapped',
             'next segment in the template unmapped',
             'SEQ being reverse complemented',
             'SEQ of the next segment in the template being reversed',
             'the first segment in the template',
             'the last segment in the template',
             'secondary alignment',
             'not passing quality controls',
             'PCR or optical duplicate']
    return([bitlist[i] for i in range(len(bitlist)) if flag & 2**(i)])

    
def print_help ():
    print '''
    trimviz takes a random sample of trimmed reads from a fastq file,
    looks up the same reads in an untrimmed fastq file and visualises the
    trimmed reads with respect to surrounding base call quality values and
    adapter sequence. In soft-clipping mode, trimviz will instead
    visualize the soft-clipping of reads by an aligner.
    
    Usage:
    ./trimsvis.py FQ -o/-O output_dir -u untrimmed.fq.gz -t trimmed.fq.gz [ -b align.bam -g reference.fa ]
    ./trimsvis.py SC -o/-O output_dir -u unaligned.fq.gz -b align.bam -g reference.fa
    
    trimsvis.py FQ        Fastq-fastq comparison. Bam file and genome fasta file can be optionally given to view the mapping outcomes for trimmed reads.
    trimsvis.py SC        Treat soft clipping as the trimming of interest. Bam and genome fasta file are required, with only one fastq file.
    
    options:
    -o/--out_dir          Directory for output. If it already exists, an error will be generated.
    -O/--out_dir_fat      Directory for output (retain temporary files; choose this option to keep the sub-sampled fastq files)
    -u/--untrimmed        Untrimmed (original) fastq file name. In SC mode, this should be the fastq file that was directly input into the aligner.
    -t/--trimmed          Trimmed fastq file name
    -U/--untrimmed_r2     Read2 fastq file corresponding to -u file (not implemented yet)
    -T/--trimmed_r2       Read2 fastq file corresponding to -t file (not implemented yet)
    -b/--bam              Bam file (optional in FQ mode; required in SC mode). Only R1 alignments will be extracted (or only R2 if -e is set).
    -g/--genome_fasta     Fasta file of genome sequence (required if using .bam alignment)
    -c/--classes          ['all'] Which trim-classes to analyse. Either 'all' or one/several of 'uncut','5pcut','3pcut','removed','generated_warning','indel'
    -a/--adapt:           ['AAAAAATGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTTCAGAGTTCTACAGTCCGACGATC'] comma-separated adapter sequences to highlight
    -A/--adaptfile        Text file containing adapter sequences
    -n/--sample_size      [50000] internal parameter: max reads to subsample in file (should be >> -m and -v, especially if only a small proportion are trimmed)
    -v/--nvis             [20] number of reads in each category (or in total if -R is set) to use for detailed individual plots
    -w/--heatmap_reads    [200] number of reads to plot in heatmaps
    -f/--agg_flank        [20] number of flanking nucleotides around trim point to plot in heatmaps 
    -r/--rid_file         File of read-ids to select, instead of using random sampling
    -s/--rseed            [1] random seed for sampling
    
    flags:
    -R/--representative   Ignore read-classes and take a representative sample. (This often results in untrimmed reads dominating the visualization)
    -z/--gzipped          Assume fastq files are gzipped (default behaviour is to guess via .gz file extension)
    -e/--read2only        If set, only alignments with the 'read-2' flag will be extracted from the .bam file (default: only read-1 is extracted)
    -h/--help:            Print help
     
    Requires:
    Rscript
    seqtk
    zcat
    python libs:
    getopt, subprocess, random, re, sys, os, tempfile, time, gzip, pipes, pysam
    R libs:
    ggplot2,dplyr,ape,reshape2,gridExtra
    '''

def makeReport(mode, out_DN, trimClassTbl, lenpre, orig_FN1, proc_FN1, orig_FN2, proc_FN2, bam_FN, gfasta_FN, PDF_FNs=[], TPCs=[]):      
    report=list()                                                                            
    report.append ('''<!DOCTYPE html>
    <style>
    table.mytable-marg{
        border-collapse: collapse;
    }
    table.mytable-marg td, table.mytable-marg th{
      border: 1px solid #ccc;
      text-align: left;
      padding-right: 18px; 
      padding-left: 10px;
      font-weight: normal
    }
    </style>
    ''')
    
    report.append( '<h1> Trimviz trimming summary: %s </h1> <br><hr/><br>' % mode )
    
    VAR1={'Untrimmed R1 fastq file':orig_FN1,
          'Trimmed R1 fastq file':proc_FN1,
          'Untrimmed R2 fastq file':orig_FN2,
          'Trimmed R2 fastq file':proc_FN2,
          'bam file':bam_FN,
          'genome fasta file':gfasta_FN}
    
    report.append ('<br>'.join( (k + ': ' + v) for (k, v) in VAR1.items() if not v=='')  )

    htmlTbl = '''
    <h3> Summary of trimming frequency </h3>
    <table class="mytable-marg" ><tr><th><b>Trim class</b></th><th><b>Number of reads</b></th></tr>
    '''
    for (cls, ncls) in trimClassTbl.items():
        htmlTbl += '<tr><th> %s </th><th> %d </th>' % (cls , ncls)
    htmlTbl += '<tr><th><b>Tot unique</b></th><th><b>%d</b></th></tr>' % (lenpre)     
    htmlTbl += '</table><br><hr/>'
    report.append (htmlTbl)
    
    for c in TPCs:
        htmlblock = '<p> %s trim position profile: <img src="%s_hist.png" title="trimming profile"' % (c,c)
        htmlblock += 'alt="" style="display: block; margin: auto;" /></p><br>'
        report.append( htmlblock )
        
    #<!-- VAR7 example: ['3p','5p']; (repeat html block len(VAR7) times) --> 
    
    
    pdfs = {'TVheatmap_S.pdf':'Heatmap, reads clustered by untrimmed read sequence',
            'TVheatmap_Q.pdf':'Heatmap, reads clustered by quality patterns',
            'TVheatmap_G.pdf':'Heatmaps, clustered by local reference sequence',
            'indiv_reads.pdf':'Individual read trimming, grouped by trimming class'}
 
    for (fn) in PDF_FNs:
        htmlblock = '<hr/><details> <summary> %s </summary> <br>' % (pdfs[fn])
        htmlblock += '<embed src= %s type="application/pdf" width="100%%" height="600px" />' % (fn)
        htmlblock += '<br> </details> <p> <br></p>'
        report.append( htmlblock )
    
    report.append('<br>')
    return(''.join(report))


if __name__ == "__main__":
    main()