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
        s  = self.f.next().rstrip()
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
    options, remainder = getopt.getopt(sys.argv[1:], 'o:p:O:P:B:c:n:v:s:a:A:x:d:g:t:r:bzhf:e:m:w', ['original=',
                                                                                   'processed=',
                                                                                   'originalR2=',
                                                                                   'processedR2=',
                                                                                   'bam_file=',
                                                                                   'classes=',   # cut,uncut,removed
                                                                                   'sample_size=', 
                                                                                   'nvis=',
                                                                                   'seed=',
                                                                                   'adapt=',
                                                                                   'adaptfile=',
                                                                                   'plotfiles_prefix=',
                                                                                   'data_out_file=',
                                                                                   'genome_fasta=',
                                                                                   'genome_gtf=',
                                                                                   'rid_file='
                                                                                   'balance_classes',
                                                                                   'gzipped',
                                                                                   'help',
                                                                                   'aggflank=',
                                                                                   'aggPlot=',
                                                                                   'maxNaggplot=',
                                                                                   'read2only'])
    # set defaults
    orig_FN1 = str('')  # o
    proc_FN1 = str('')  # p
    orig_FN2 = str('')  # O
    proc_FN2 = str('')  # P
    bam_FN = str('')    # B
    coi=list()          # c
    target_n_pre=50000  # n 
    nvis=20             # v
    rseed = 1           # s
    madapt=list()       # a   ...(madapt[0])
    adapt_FN = str('')  # A
    gr_FN = str('')     # x  
    tblOut_FN = str('') # d
    gfasta_FN = str('') # g
    gtf_FN = str('')    # t
    rid_FN = str('')    # r
    balance = True      # b
    gzipped = True      # z
    nfirstpass = 4000   # ?   ... TODO add to args
    class_opts = ['uncut','5pcut','3pcut','removed']  # codify as 0,1 and 2
    coi=['all']
    aggFlank=20         # f
    aggGr_FN = 'aggPlot.pdf' # e
    maxAggN = 200       # m
    read2only=False     #w
    
    for opt, arg in options:
        if opt in ('-o', '--original'):
            orig_FN1 = arg
        elif opt in ('-p', '--processed'):
            proc_FN1 = arg
        elif opt in ('-O', '--originalR2'):
            orig_FN2 = arg
        elif opt in ('-P', '--processedR2'):
            proc_FN2 = arg
        elif opt in ('-c', '--classes'):
            coi_raw = arg.split(',')
            coi = [x for x in coi_raw if x in class_opts]
            if len(coi) < len (coi_raw) and coi_raw != ['all']:
                print ("Warning, values given in -c should include only 'uncut', 'cut' and/or 'removed', comma-separated; or 'all'")
        elif opt in ('-n', 'sample_size'):
            target_n_pre = int(arg)
        elif opt in ('-v', '--nvis'):
            nvis = int(arg)
        elif opt in ('-s', '--seed'):
            rseed = int(arg)
        elif opt in ('-a', '--adapt'):
            madapt = arg.split(',')
        elif opt in ('-A', '--adaptfile'):
            adapt_FN = arg
        elif opt in ('-x', '--plotfiles_prefix'):
            gr_FN = arg
        elif opt in ('-d', '--data_out_file'):
            tblOut_FN = arg
        elif opt in ('-B', '--bam_file'):
            bam_FN = arg
        elif opt in ('-g', '--genome_fasta'):
            gfasta_FN = arg # g
        elif opt in ('-t', '--genome_gtf'):
            gtf_FN = str('')    # t
        elif opt in ('-r', '--rid_file'):
            rid_FN = arg    # t
        elif opt in ('-b', '--balance_classes_for_plot'):
            balance = True
        elif opt in ('-z', '--gzipped'):
            gzipped = True
        elif opt in ('-h', '--help'):
            print_help()
        elif opt in ('-f', '--aggflank'):
            aggFlank = int(arg)
        elif opt in ('-e', '--aggPlot'):
            aggGr_FN=arg
        elif opt in ('-m', '--maxNaggplot'):
            maxAggN = int(arg)
        elif opt in ('-w', '--read2only'):
            read2only=True    
    # Open a file.  Use gzip based on filename, or 'gzipped' flag
    def gzopen(fname):
        if gzipped or fname.endswith('gz'):
            return gzip.open(fname, 'rb')
        else:
            return open(fname, 'r')
    
    ###########################################################################################
    #  MAIN part 1:  arguments logic
    ###########################################################################################
    
    softClipping=False
    if not os.path.isfile(proc_FN1):
        if (len(proc_FN1) > 0):
            print('Could not find processed R1 file ' + proc_FN1 + '. Exiting.')
            print_help()
            exit()
        if not os.path.isfile(bam_FN ):
            print('No processed R1 file given; and no valid bam filename was given in lieu. Exiting.')
            print_help()
            exit()
        else:
            print ('Bam file given but not processed fastq file. Will treat soft-clipping in bam file as the trimming to analyze.')
            softClipping=True
    if not os.path.isfile(orig_FN1): 
        print('Could not find file: ' + orig_FN1 + '. Exiting.')
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
    if not os.path.exists('trimVisTmpFiles'):
        os.makedirs('trimVisTmpFiles') 
    tmpPre1_FN = os.curdir+'/trimVisTmpFiles/tmpPre1.fq.gz'
    tmpPost1_FN = os.curdir+'/trimVisTmpFiles/tmpPost1.fq.gz'
    tmpPre2_FN = os.curdir+'/trimVisTmpFiles/tmpPre2.fq.gz'
    tmpPost2_FN = os.curdir+'/trimVisTmpFiles/tmpPost2.fq.gz'
    readIDs1_FN = os.curdir+'/trimVisTmpFiles/readIDsPre.lst'
    
    # sample original fq:
    if rid_FN == '':
        print 'target number of reads (orig fastq): %d' % (target_n_pre)
        pcnt_sgn = """%"""
        # split output into a temp fq file (tmpPre1_FN) and a readname list file (readIDs1_FN)
        # remove @, anything after a space in readname; append a tab char to improve fgrep specificity for bam-searching 
        cmd3 = "seqtk sample -s%d %s %d | tee >(awk '1 == NR %s 4' | sed 's/@//'  | sed 's/[ \/].*$//' | sed 's/$/\t/' > %s ) | gzip > %s" % (rseed,
                                                                                                                          pipes.quote(orig_FN1),
                                                                                                                          target_n_pre,
                                                                                                                          pcnt_sgn,
                                                                                                                          readIDs1_FN,
                                                                                                                          pipes.quote(tmpPre1_FN))       
    else: # user-specified RIDs
        cmd2 = "cat %s | sed 's/@//'  | sed 's/[ \/].*$//' | sed 's/$/\t/' > %s" % (rid_FN, readIDs1_FN) # clean up user-specified rids for seqtk
        with open('./trimVisTmpFiles/tmp_bash2.sh', 'w') as fout:
            print >> fout, cmd2
        sout2 = subprocess.check_output(['bash', './trimVisTmpFiles/tmp_bash2.sh'])
        print (sout2)
        cmd3 = "seqtk subseq %s %s | gzip > %s " % (pipes.quote(orig_FN1), readIDs1_FN, tmpPre1_FN)
    
    with open('./trimVisTmpFiles/tmp_bash3.sh', 'w') as fout:
            print >> fout, cmd3
            
    # retrieve from ORIG FQ using seqtk
    print ('sampling from original fastq using seqtk...')
    sout3 = subprocess.check_output(['bash', './trimVisTmpFiles/tmp_bash3.sh'])
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
    rid_class={'uncut':[], '5pcut':[],'3pcut':[], 'removed':[]} # will also classify them into the 3 classes and make ID list for each class (-> rid_class)
        
    ####################################################################################
    # MAIN part 2.1: (optional): extract from bam & fasta; save as temp file & dict (fastaD)
    ####################################################################################
 
    if (bam_FN != ''):
        sub_bamFN = './trimVisTmpFiles/bamEntries1.sam'
        cmd4A = 'samtools view -H %s > %s' % (bam_FN,  sub_bamFN)
        cmd4B = 'samtools view %s | fgrep -f %s | awk \'NF > 8{print} 1\' >> %s' % (bam_FN, readIDs1_FN, sub_bamFN)  # |  cut -f1-6,9
        with open('./trimVisTmpFiles/tmp_bash6.sh', 'w') as fout:
            print >> fout, cmd4A
            print >> fout, 'sleep 1'
            print >> fout, cmd4B
        print ('searching bam file using samtools and fgrep...')
        time.sleep(1)
        sout = subprocess.check_output(['bash', './trimVisTmpFiles/tmp_bash6.sh'])
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
    # MAIN part 2.2: (optional): extract from post-trimmed fastq; merge into dict 'both'
    ####################################################################################  
    
    if not softClipping:     
        # retrieve same reads (in readIDs1_FN) from TRIMMED FQ (proc_FN1) using seqtk
        cmd5 = "seqtk subseq %s %s | gzip > %s " % (proc_FN1, readIDs1_FN, tmpPost1_FN)
        with open('./trimVisTmpFiles/tmp_bash4.sh', 'w') as fout:
            print >> fout, cmd5
        print ('Finding the reads in trimmed fastq file using seqtk...')
        sout4 = subprocess.check_output(['bash', './trimVisTmpFiles/tmp_bash4.sh'])
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
            
        #########################################################
        #  Guess 5' and 3' trim sites using both seq and qual strings
        #  add as elements [4] and [5] in 'both'
        #########################################################
    
        removeread = list()
        maxReadLen=0
        for r, rdat in both.items():
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
                        print "Warning - ambiguous trimming encountered for read %s - removing" % (r)
                        print len(fps)
                        print rdat
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
            if cigTuples[0][0] == 4:  # leading soft-clipped. TODO: hard clipped is '5'
                scLead += cigTuples[0][1]
                nlpad = max(0, scLead - refBlocks[0][0]) # if runs off edge of contig:
                scLead -= nlpad
                refBlocks = [(refBlocks[0][0]-scLead, refBlocks[0][0])] +  refBlocks   
            if cigTuples[-1][0]==4:
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
                    refBlocks =  [(lhs, refBlocks[0][0])] + refBlocks
                if trimTrail > 0:
                    refBlocks = refBlocks + [(refBlocks[-1][1], rhs)] # add trimmed bits to outer segments of refBlocks
            lpad = 'X' * nlpad
            rpad = 'X' * nrpad                
            gblocks = [genome_fa.fetch(refname, rb[0], rb[1]) for rb in refBlocks]
            gSeg = lpad + ''.join(gblocks) + rpad
            if r.is_reverse:
                gSeg = srevcomp(gSeg)
                joinLocs=[x[0] for x in refBlocks[::-1]]
            else:
                joinLocs=[x[1] for x in refBlocks]
            bamInfo[id2] = [gSeg, joinLocs]
            
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
    if (balance and coi != ['all']):
        toPlot=[]
        for cls in coi:
            ids = rid_class[cls]
            for r in ids:
                #if r in lookupCls:
                #    print "warning: read %s is in multiple classes: %s and %s " % (r, lookupCls[r], cls)
                lookupCls[r]=cls
            print ("There were %d reads classed as %s") % (len(rid_class[cls]) , cls)
            rs = random.sample(ids, k=min( len(ids) , nvis ))
            toPlot = toPlot + rs
    else:
        toPlot = random.sample(keys(both), nvis)
    

    
    #with tempfile.NamedTemporaryFile(delete = False) as tempf:
    tempfname = './trimVisTmpFiles/trimviz_readData.tsv'
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
    print "FILE1 READY:"
    print tempfname
  
    ##################################################################################
    #   MAIN    part 6: prepare aggregate / trim-anchored read matrices
    ##################################################################################            
      
    runsheet={'3pcut':'./trimVisTmpFiles/seq3psites.txt', '5pcut':'./trimVisTmpFiles/seq5psites.txt'}
    for cls, clsfile in runsheet.items():
        i=0
        with open(clsfile, 'w') as fout:
            bs = range(aggFlank*2+1)
            if bam_FN == '':
                print >> fout, '\t'.join(['readID', 'tpCutPos'] + ['s'+str(i+1) for i in bs] + ['q'+str(i+1) for i in bs])
            else:
                print >> fout, '\t'.join(['readID', 'tpCutPos'] + ['s'+str(i+1) for i in bs] + ['q'+str(i+1) for i in bs] + ['g'+str(i+1) for i in bs])
            for r in rid_class[cls]:
                if i < 200: #if r in toPlot:
                    rdat=both[r]
                    ps=padstr(rdat[2], rdat[5]-1, aggFlank, r)  # pre-seq; 3' trim site (0-based so should not exceed len(rdat[2]) ); flanking seq = 20. returns list
                    pq=padstr(rdat[3], rdat[5]-1, aggFlank, r)
                    if bam_FN != '':                      
                        if r in bamInfo:
                            gSeq=bamInfo[r][0]
                        else:
                            gSeq='N'*len(rdat[2])
                        gs=padstr(gSeq, rdat[5]-1, aggFlank, r)
                        print >> fout, '\t'.join([r, str(rdat[5])] + ps + [str(ord(x)) for x in pq] + gs)
                    else:
                        print >> fout, '\t'.join([r, str(rdat[5])] + ps + [str(ord(x)) for x in pq])
                    
 
        
    ###########################
    # << CALL R SCRIPT HERE >>#
    ###########################
    
    cmd6 = 'Rscript ' +'./graph_ts.R '+ tempfname + ' ' + gr_FN + ' ' + aggGr_FN + ' ' + str(maxAggN)  #os.curdir
    print (cmd6)
    rout = subprocess.check_output(cmd6, shell=True)
    
    print "R stdout:", rout
        
   
# <<<<<<<<<<<<<<<<< END MAIN >>>>>>>>>>>>>>>>>>>     
 
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
    adapter sequence. If a bam file is provided AND the processed fastq
    file (-p) is ommitted, trimviz will instead visualize the soft-clipping
    of reads by the aligner.
    
    Usage:
    ./trimsvis.py -o orig.fq.gz -x plotfiles_prefix [ -p processed.fq.gz | -b align.bam -g reference.fa ]
    options:
    -a/--adapt:           ['AAAAAATGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTTCAGAGTTCTACAGTCCGACGATC'] comma-separated adapter sequences to highlight
    -A/--adaptfile        text file containing adapter sequences
    -h/--help:            print help
    -g/--genome_fasta     fasta file of genome sequence (required if using .bam alignment)
    -m/--maxNaggplot      [200]
    -v/--nvis             [20]
    -r/--rid_file         file of read-ids to select, instead of using random sampling
    -x/--plotfiles_prefix prefix of output pdf file (individual read visualization)
    -e/--aggPlot          prefix of output pdf file (clustered read heatmaps)
    -n/sample_size        [50000] internal parameter: max reads to subsample in file (should be >> -m and -v, especially if only a small proportion are trimmed)
     
    Requires:
    Rscript
    seqtk
    zcat
    python libs:
    getopt, subprocess, random, re, sys, os, tempfile, time, gzip, pipes, pysam
    R libs:
    ggplot2,dplyr,ape,reshape2,gridExtra
    '''

if __name__ == "__main__":
    main()
    
