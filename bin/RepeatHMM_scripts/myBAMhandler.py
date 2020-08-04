
import re
import os
import sys
import string
import math
import copy

import numpy
import time
import resource
import peakutils
import argparse

import logging

from scipy.stats import norm

import heapq

from Bio import pairwise2
import Bio

from . import getAlignment
from . import myHMM
from . import myGaussianMixtureModel
from . import printHMMmatrix
#from . import myRepeatReAlignment

#from .myheader import *
from . import myheader

def myReadTxtFile(filename):
    f = open(filename, 'r')

    data = f.readlines()
    while len(data) > 0 and string.strip(data[-1]) == "":
        data = data[:-1]
    f.close()

    return data


def myWriteTxtFile(mlist, filename):
    f = open(filename, 'w')

    for ml in mlist:
        if ml[-1] == '\n':
            f.write(ml)
        else:
            f.write(ml + '\n')

    f.close()


def get_gLoc(repeatname, commonOptions):
    repeatname = repeatname.lower()

    # 0       1         2        3      4       5      6
    #mdict[repname] = [chr, start_pos, end_pos, reppat, strand, range, others]

    gLoc = commonOptions['gLoc']
    nece_field_num = 5
    if gLoc.has_key(repeatname):
        curloc = gLoc[repeatname]
    else:
        curloc = ['', '', '', '',   '', '', '']

    newinfo = commonOptions['specifiedRepeatInfo']
    infospt = newinfo.split('/')
    if len(infospt) < nece_field_num:
        logging.error("Error: wrong input for the gene of interes: %s whose length is %d less than the expected length %d" % (
            newinfo, len(infospt), nece_field_num))
        sys.exit(101)
    for i in range(len(infospt)):
        curnew = string.strip(infospt[i])
        if not curnew == '':
            curloc[i] = curnew
    errorstr = ''
    for i in range(len(curloc)):
        if isinstance(curloc[i], str) and string.strip(curloc[i]) == '' and i < nece_field_num:
            errorstr += ('Error no information for %d\n' % i)
    if not errorstr == '':
        logging.error(errorstr)
        if commonOptions['outlog'] <= myheader.M_WARNING:
            print(errorstr, curloc, repeatname)
        sys.exit(102)

    return get_Loc1(curloc, commonOptions)


def get_Loc1(curloc, commonOptions):
    curloc[1] = int(curloc[1])
    curloc[2] = int(curloc[2])
    retOptions = {}
    retOptions['repeat_start_end'] = [curloc[1], curloc[2]]
    retOptions['chr'] = curloc[0]
    retOptions['repPat'] = curloc[3]
    retOptions['forw_rerv'] = curloc[4]
    retOptions['mgloc'] = curloc

    if commonOptions['SeqTech'] == "Illumina":
        rep_up_down_size = 300
    elif commonOptions['SeqTech'] == "Pacbio":
        rep_up_down_size = 1000
    elif commonOptions['SeqTech'] == "Nanopore":
        rep_up_down_size = 1000
    else:
        rep_up_down_size = 1000
    retOptions['gene_start_end'] = [curloc[1] - rep_up_down_size, curloc[2] + rep_up_down_size]
    if retOptions['gene_start_end'][0] < 1:
        retOptions['gene_start_end'] = 1

    return retOptions


def insert_n_for_flanking(oldstr, rep_predata, rep_sufdata, repeatFlankLength, repPat, otherinfo):
    #min_flank_len = 16;

    lengthmore = 5 * len(repPat)
    if lengthmore > 100:
        lengthmore = 100
    query_pos = [[0, repeatFlankLength + lengthmore], [-repeatFlankLength - lengthmore, ]]
    query_pre_suf = [oldstr[query_pos[0][0]:query_pos[0][1]], oldstr[query_pos[1][0]:]]

    if Bio.__version__ == '1.66':
        match, mismatch, gap, gap_extension = 2, -1, -1, -10  # from Bio import pairwise2;
    else:
        match, mismatch, gap, gap_extension = 1, -1, -1, -1
    twopairs = [[query_pre_suf[0], rep_predata], [query_pre_suf[1], rep_sufdata]]
    detectpos = []
    detail2 = []
    for curp_ind in range(len(twopairs)):
        curp = twopairs[curp_ind]
        if len(curp[1]) < 5:
            detectpos.append([None, None])
            continue

        if curp_ind == 0:
            curp[1] = curp[1][:-1]
        else:
            curp[1] = curp[1][1:]

        la = pairwise2.align.localms(curp[0], curp[1], match, mismatch, gap, gap_extension)
        seq1, seq2, score, abegin, aend = '', '', 0, -1, -1
        for a in la:
            p, q, s, b, e = a
            if score < s:
                seq1, seq2, score, abegin, aend = p, q, s, b, e
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print('Fatal not equal', seq1, curp[0])
            print('Fatal not equal', seq2, curp[1])
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print(seq1, curp[0])
            print(seq2, curp[1])

        align_len = aend - abegin
        align_len = len(curp[0])
        if align_len == 0 and len(curp[1]) > 0:
            align_len = len(curp[1])
        elif align_len > len(curp[1]) and len(curp[1]) > myheader.min_flank_len:
            align_len = len(curp[1])
        if align_len < myheader.min_flank_len:
            align_len = myheader.min_flank_len

        if myheader.cur_M_STAT <= myheader.M_INFO:
            print(seq1, '<', curp[0])
            print(seq2, '<', curp[1])
            print(abegin, aend, len(seq1), len(seq2))

        while abegin < aend and (seq2[abegin] == '-' or (not seq2[abegin] == seq1[abegin])):
            abegin += 1
        while abegin < aend and (seq2[aend - 1] == '-' or (not seq2[aend - 1] == seq1[aend - 1])):
            aend -= 1

        numsame = 0
        numnot = 0
        numgap = 0
        query_start = abegin
        query_end = abegin
        for ai in range(abegin, aend):
            if not seq1[ai] == '-':
                query_end += 1
                if seq1[ai] == seq2[ai]:
                    numsame += 1
                elif seq2[ai] == '-':
                    numgap += 1
                else:
                    numnot += 1
            else:
                if seq2[ai] == '-':
                    if myheader.cur_M_STAT <= myheader.M_ERROR:
                        print('Fatal!!! both gap', seq1, seq2)
                else:
                    numgap += 1
        identical_fraction = 0
        if align_len > 0:
            identical_fraction = numsame / float(align_len)
        if align_len > myheader.min_flank_len and identical_fraction > 0.75:
            detectpos.append([query_start, query_end])
        else:
            if myheader.cur_M_STAT <= myheader.M_INFO:
                print('xxxxxxxxxx', oldstr)
            detectpos.append([None, None])
            #detail2.append(('Identical(%d) of %d bps(per=%.2f) for two sequences with len(%d, %d)' % (numsame, align_len, identical_fraction, len(seq1), len(seq2))))
            detail2.append(('%d/%d=%.2f for len(%d, %d)' %
                            (numsame, align_len, identical_fraction, len(curp[0]), len(curp[1]))))
            if align_len == 0:
                if myheader.cur_M_STAT <= myheader.M_ERROR:
                    print('Error!!! zero alignment ', curp, aend, abegin)
                    logging.error('Error!!! zero alignment ' + str(curp) +
                                  ' ' + str(aend) + ' ' + str(abegin))
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print('info', seq1, detectpos[-1])
            print('    ', seq2, numsame, len(curp[1]), align_len,)
            if align_len > 0:
                print(round(numsame / float(align_len), 2))
            else:
                print('')

    if myheader.cur_M_STAT <= myheader.M_INFO:
        print('seq', oldstr)
    if (detectpos[1][0] is not None and detectpos[0][1] is None):
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print(detectpos)
    else:
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print('NO flank', detectpos, otherinfo, detail2)
    num_n = 10
    if detectpos[1][0] is not None:
        curpos = -len(query_pre_suf[1]) + detectpos[1][0]
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print('suf', oldstr[curpos:], detectpos[1])
        oldstr = oldstr[:curpos] + ('N' * (num_n * len(repPat))) + oldstr[curpos:]
    if detectpos[0][1] is not None:
        curpos = detectpos[0][1]
        if myheader.cur_M_STAT <= myheader.M_INFO:
            print('pre', oldstr[:(curpos + 1)], detectpos[0])
        oldstr = oldstr[:(curpos)] + ('N' * (num_n * len(repPat))) + oldstr[(curpos):]
    if myheader.cur_M_STAT <= myheader.M_INFO:
        print('n-seq', oldstr)
    return oldstr


def getUnsymAlignAndHMM(repPat, forw_rerv, repeatFlankLength, hmmoptions, queryrep, commonOptions, otherinfo):
    rep_predata, rep_sufdata = commonOptions['rep_flanking_data']

    #queryrep = queryrep.replace('N', ''); queryrep = queryrep.replace('n', '')
    queryrep = queryrep.replace('n', 'N')

    if commonOptions['isGapCorrection'] == 0:
        unewstr = queryrep
    else:
        match, mismatch, gap_in_perf, cur_gap_in_read, gap_before_after = commonOptions['MatchInfo']
        mult_gap_in_perf = 3
        if len(queryrep[repeatFlankLength:(len(queryrep) - repeatFlankLength)]) < 150:
            if cur_gap_in_read < gap_in_perf * mult_gap_in_perf:
                cur_gap_in_read = gap_in_perf * mult_gap_in_perf

        if (len(queryrep) - 2 * repeatFlankLength) > 10:
            unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatFlankLength:(len(
                queryrep) - repeatFlankLength)], forw_rerv, match, mismatch, gap_in_perf, cur_gap_in_read, gap_before_after, mismatchstr=hmmoptions[8][0], mismnum=hmmoptions[8][2]);
        else:
            unewstr = queryrep
    unewstr = unewstr.replace('-', '')

    if repeatFlankLength > 0 and commonOptions['isGapCorrection'] > 0:
        unewstr = queryrep[:repeatFlankLength] + unewstr + \
            queryrep[(len(queryrep) - repeatFlankLength):]

    unewstr = insert_n_for_flanking(unewstr, rep_predata, rep_sufdata,
                                    repeatFlankLength, repPat, otherinfo)
    # print 'after_n_insert', unewstr

    if len(unewstr) > 1:
        newstr, pre0, predstats = myHMM.hmmpred(
            unewstr, repPat, forw_rerv, hmmoptions, commonOptions)
    else:
        newstr, pre0, predstats = '', 0, ''

    return [newstr, pre0, predstats]


def getGene(repeatName, chr, gene_start_end, unique_file_id, analysis_file_id, hgfn, flanking_len, specifiedOptions):  # ='hg38.fa'
    alignfolder = specifiedOptions['align']  # 'align/'
    if not os.path.isdir(alignfolder):
        os.system('mkdir ' + alignfolder)

    fastafile = alignfolder + repeatName + unique_file_id + '.fasta'

    #flanking_len = 30;
    #if flanking_len<7: flanking_len = 7
    if flanking_len < myheader.min_flank_len:
        flanking_len = myheader.min_flank_len
    interest_region_pos = [[gene_start_end[0] - flanking_len, gene_start_end[0] - 1], [
        gene_start_end[0], gene_start_end[1]], [gene_start_end[1] + 1, gene_start_end[1] + flanking_len]]
    my_interest_region = []
    for cur_pos in interest_region_pos:
        # print cur_pos, gene_start_end, chr, hg_reference_and_index, hgfn, chr, fastafile
        if cur_pos[0] < 0:
            cur_pos[0] = 0
        if cur_pos[1] < 0:
            cur_pos[1] = 0

        #get_alg_cmd = 'samtools faidx '+hg_reference_and_index+'/'+hgfn+' '+ chr+':'+str(cur_pos[0])+'-'+str(cur_pos[1])+' > '+fastafile
        get_alg_cmd = 'samtools faidx ' + hgfn + ' ' + chr + ':' + \
            str(cur_pos[0]) + '-' + str(cur_pos[1]) + ' > ' + fastafile
        # print get_alg_cmd

        os.system(get_alg_cmd)
        if not os.path.isfile(fastafile):
            logging.error('Cannot produce ' + fastafile + ' for ' + repeatName)
            sys.exit(1)
        fastadata = myReadTxtFile(fastafile)
        mfadata = ''
        for li in fastadata:
            if li[0] == '>':
                continue
            if li[-1] == '\n':
                li = li[:-1]
            mfadata = mfadata + li

        os.system('rm ' + fastafile)
        my_interest_region.append(mfadata.upper())

    return my_interest_region


def getHMMOptions(repeatFlankLength, repPat, forw_rerv, commonOptions):
    hmmoptions = []
    hmmoptions = myHMM.getTransition_start_emission_prob(repPat, commonOptions)
    return hmmoptions


def getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions):
    logging.info(moreOptions['chr'] + ' ' + str(moreOptions['repeat_start_end']))
    chr = moreOptions['chr']
    repeatName = moreOptions['repeatName']
    gene_start_end = moreOptions['gene_start_end']
    repeat_start_end = moreOptions['repeat_start_end']
    repPat = moreOptions['repPat']
    forw_rerv = moreOptions['forw_rerv']

    bamfile = specifiedOptions['bamfile']
    unique_file_id = specifiedOptions['unique_file_id']
    analysis_file_id = specifiedOptions['analysis_file_id']

    isGapCorrection = commonOptions['isGapCorrection']
    repeatFlankLength = commonOptions['repeatFlankLength']
    MinSup = commonOptions['MinSup']

    len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
    logging.info("len_repPat=" + str(len_repPat))
    # print repeatName,
    alignfolder = specifiedOptions['align']  # 'align/'
    if not os.path.isdir(alignfolder):
        os.system('mkdir ' + alignfolder)

    ref_repeat = (repeat_start_end[1] - repeat_start_end[0] + 1) / float(len_repPat)  # 3.0

    alignfile = alignfolder + repeatName + unique_file_id + '.alignment.txt'
    get_alg_cmd = 'samtools view ' + bamfile + ' ' + chr + ':' + \
        str(gene_start_end[0]) + '-' + str(gene_start_end[1]) + ' > ' + alignfile
    if 'thread' not in specifiedOptions:
        logging.info('Running ' + get_alg_cmd)
    os.system(get_alg_cmd)
    if os.path.getsize(alignfile) == 0:
        if commonOptions['outlog'] <= myheader.M_WARNING:
            logging.info(get_alg_cmd + '\n')
            logging.info('The file %s have zero size in the function of getRepeatForGivenGene.\nTry without chr' % alignfile)
            #print ('The file %s have zero size\nTry without chr' % alignfile)
        get_alg_cmd = 'samtools view ' + bamfile + ' ' + \
            chr[3:] + ':' + str(gene_start_end[0]) + '-' + \
            str(gene_start_end[1]) + ' > ' + alignfile
        if commonOptions['outlog'] <= myheader.M_INFO and ('thread' not in specifiedOptions):
            logging.info('Running ' + get_alg_cmd)
        os.system(get_alg_cmd)
    if commonOptions['outlog'] <= myheader.M_INFO:
        logging.info('Produced ' + alignfile + ' done!')

    if (not os.path.isfile(alignfile)) or os.path.getsize(alignfile) == 0:
        if commonOptions['outlog'] <= myheader.M_FATAL:
            logging.error('!!!!Cannot produce ' + alignfile + ' for ' + repeatName)
            # sys.exit(1)
            os.system('rm ' + alignfile)
            return None
    aligndata = myReadTxtFile(alignfile)
    os.system('rm ' + alignfile)

    repregion_len_threhold = len_repPat  # 3;

    predata, mfadata, sufdata = getGene(repeatName, chr, gene_start_end, unique_file_id,
                                        analysis_file_id, commonOptions['hgfile'], repeatFlankLength, specifiedOptions)
    rep_predata, rep_mfadata, rep_sufdata = getGene(
        repeatName, chr, repeat_start_end, unique_file_id, analysis_file_id, commonOptions['hgfile'], repeatFlankLength, specifiedOptions)

    commonOptions['rep_flanking_data'] = rep_predata, rep_sufdata

    wrongalign = 0

    hmmoptions = getHMMOptions(repeatFlankLength, repPat, forw_rerv, commonOptions)

    repeats = []
    repeats_dict = {}
    ids = []
    for line in aligndata:
        lsp = line.split('\t')
        readid = lsp[0]
        cchr = lsp[2]

        if not (chr == cchr or (len(chr) > 3 and chr[3:] == cchr) or (len(cchr) > 3 and cchr[3:] == chr)):
            continue

        pos = int(lsp[3])
        aligninfo = lsp[5]
        aainfo = lsp[9]

        if pos > repeat_start_end[0] - repeatFlankLength:
            wrongalign += 1
            # continue;
            #logging.error('The start pos in ref Genome is greater than the start position of repeats' + str(pos) +' ' + str(repeat_start_end[0]));
        if not (cchr == chr or cchr == chr[3:]):
            logging.error('Not same ' + cchr + ' ' + chr)
            continue

        numreg = re.compile('\d+')
        numinfo = numreg.findall(aligninfo)

        mdireg = re.compile('[MIDNSHPX=]{1}')
        mdiinfo = mdireg.findall(aligninfo)

        if not len(numinfo) == len(mdiinfo):
            logging.error('Num is equal to mid' + str(len(numinfo)) + ' ' + str(len(mdiinfo)))
            continue

        queryind = 0
        hpadd = 0
        queryrep = ''
        longer = False
        query_start_ind = None
        query_end_ind = None

        for n1ind in range(len(numinfo)):
            n1 = int(numinfo[n1ind])
            mdi = mdiinfo[n1ind]

            for n1i in range(n1):
                qrepadd = False
                if mdi == 'M' or mdi == '=' or mdi == 'X':
                    pos = pos + 1
                    queryind = queryind + 1
                    qrepadd = True
                elif mdi == 'I':
                    qrepadd = True
                    queryind = queryind + 1
                elif mdi == 'D':
                    pos = pos + 1
                elif mdi == 'S':
                    queryind = queryind + 1
                    qrepadd = True
                elif mdi == 'H':
                    if qrepadd:
                        hpadd += 1  # pass
                elif mdi == 'P':
                    if qrepadd:
                        hpadd += 1  # pass
                else:
                    logging.warning('Warning unknown CIGAR element ' + str(n1) + ' ' + mdi)
                if qrepadd:
                    if pos - 1 >= repeat_start_end[0] - repeatFlankLength and pos - 1 <= repeat_start_end[1] + repeatFlankLength:
                        queryrep = queryrep + aainfo[queryind - 1]
                if pos - 1 < repeat_start_end[0] - repeatFlankLength:
                    query_start_ind = queryind - 1
                if pos - 1 >= repeat_start_end[1] and pos - 1 < repeat_start_end[1] + repeatFlankLength:
                    query_end_ind = queryind

            if pos - 1 > repeat_start_end[1] + repeatFlankLength:
                longer = True

        if readid not in repeats_dict:
            repeats_dict[readid] = [query_start_ind, query_end_ind, aainfo]
        else:
            if query_start_ind is not None:
                if repeats_dict[readid][0] is None or repeats_dict[readid][0] > query_start_ind:
                    repeats_dict[readid][0] = query_start_ind
            if query_end_ind is not None:
                if repeats_dict[readid][1] is None or repeats_dict[readid][1] < query_end_ind:
                    repeats_dict[readid][1] = query_end_ind
            if len(repeats_dict[readid][2]) < len(aainfo):
                repeats_dict[readid][2] = aainfo

        if len(queryrep) >= repregion_len_threhold:
            repeats.append([longer, queryrep, lsp[0]])
            ids.append(readid)

    handleint = True
    if handleint:
        repeats = []
        ids = []
        repeatskeys = repeats_dict.keys()
        for rk in repeatskeys:
            if repeats_dict[rk][0] is None or repeats_dict[rk][1] is None:
                repeats.append([False, str(repeats_dict[rk][0]) +
                                '-to-' + str(repeats_dict[rk][1]), rk])
                ids.append(rk)
            else:
                if repeats_dict[rk][1] - repeats_dict[rk][0] > 0:
                    repeats.append([True, repeats_dict[rk][2]
                                    [repeats_dict[rk][0]:(repeats_dict[rk][1] + 1)], rk])
                    ids.append(rk)
                else:
                    if commonOptions['outlog'] <= myheader.M_WARNING:
                        print('Warning!!! negative ', rk, repeats_dict[rk][:2])
                    repeats.append([False, str(repeats_dict[rk][0]) +
                                    '-to-' + str(repeats_dict[rk][1]), rk])
                    ids.append(rk)

    rptrue = []
    rpfalse = []
    orignial = []
    for currep_ind in range(len(repeats)):
        currep = repeats[currep_ind]
        newstr = currep[1]

        pre0 = 0
        predstats = ''
        if len(newstr) < commonOptions['MaxRep'] * len_repPat:
            if currep[0]:
                    # print 'BAMhandler', repeat_start_end, chr
                newstr, pre0, predstats = getUnsymAlignAndHMM(
                    repPat, forw_rerv, repeatFlankLength, hmmoptions, currep[1], commonOptions, ids[currep_ind])
            else:
                if 'thread' not in specifiedOptions:
                    logging.warning('The sequence is partial: ' + str(len(newstr)) + ' ' + chr + ' ' + repeatName + ' ' + repPat + ' ' + str(
                        currep[0]) + ' reads name:' + currep[2] + " " + str(commonOptions['MaxRep']) + " " + str(commonOptions['MaxRep'] * len_repPat))
                    if handleint:
                        logging.warning(str(repeats_dict[currep[2]][:2]))
        else:
            logging.warning('The sequence is too long: ' + str(len(newstr)) + ' ' + chr + ' ' + repeatName + ' ' + repPat + ' ' + str(
                currep[0]) + ' reads name:' + currep[2] + " " + str(commonOptions['MaxRep']) + " " + str(commonOptions['MaxRep'] * len_repPat))
            if handleint:
                logging.warning(str(repeats_dict[currep[2]][:2]))
        orignial.append([currep[1], pre0, predstats])
        currep[1] = newstr
        if currep[0]:
            rptrue.append(len(currep[1]) / float(len_repPat))  # 3.0);
        else:
            rpfalse.append(len(currep[1]) / float(len_repPat))  # 3.0);

    rptrue.sort()
    rpfalse.sort()
    trstr = 'true ' + str(len(rptrue)) + ' ['
    for rpt in rptrue:
        trstr = trstr + ('%.0f,' % rpt)
    trstr = trstr[:-1] + ']'
    logging.debug(trstr)

    p2, allocr = myGaussianMixtureModel.get2Peaks(rptrue, MinSup, commonoptions=commonOptions)

    if len(rpfalse) > 0:
        flstr = 'fals ' + str(len(rpfalse)) + ' ['
        for rpf in rpfalse:
            flstr = flstr + ('%.0f,' % rpf)
        flstr = flstr[:-1] + ']'
        logging.debug(flstr)

    logging.info('ref_repeat ' + ('%.0f' % ref_repeat) + '\t' + repPat + '\t' + forw_rerv)

    for currep_ind in range(len(repeats)):
        currep = repeats[currep_ind]

        aaprinindex = -1
        if not (currep[0]):
            aaprinindex = 300

        logging.debug('\t' + str(currep[0]) + ' o:' + str(len(orignial[currep_ind]
                                                              [0])) + '\t' + orignial[currep_ind][0][:aaprinindex]);
        prestr = ''
        # print currep_ind, orignial[currep_ind][1], orignial[currep_ind]
        for i in range(orignial[currep_ind][1]):
            prestr = prestr + ' '
        logging.debug('\t' + str(currep[0]) + ' p:' + str(len(currep[1])
                                                          ) + '\t' + prestr + (currep[1][:aaprinindex]))

    return [repeatName, ref_repeat, p2, allocr, len(rptrue), len(rpfalse) + wrongalign]


def fixsize2(p2, ind):
    if len(p2[ind]) > 1:
        pass
    elif len(p2[ind]) == 1:
        p2[ind].append(p2[ind][0])
    else:
        p2[ind].append(0)
        p2[ind].append(0)


def addSumForAGene(p2, myret, myretdetail, mstr, ind=2):
    fixsize2(p2, ind)
    myretdetail[mstr] = (('%10s' % mstr) + ' ' + str(p2))
    # print mstr, p2
    #logging.info(mstr+' '+str(p2)+'\n')
    sys.stdout.flush()
    myret[mstr] = p2[ind]


def getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions={}):
    if "repeatName" in moreOptions:
        repeatName = moreOptions['repeatName'].lower()  # repeatName.lower()
    else:
        repeatName = commonOptions['repeatName'].lower()
        moreOptions['repeatName'] = repeatName

    retoptions = get_gLoc(repeatName, commonOptions)
    mgloc = retoptions['mgloc']
    if commonOptions['outlog'] <= myheader.M_INFO:
        print('mgloc', mgloc)

    moreOptions.update(retoptions)
    myHMM.produce_for_repPat(commonOptions, moreOptions)

    if specifiedOptions["SepbamfileTemp"] is not None:
        specifiedOptions["bamfile"] = (specifiedOptions["SepbamfileTemp"] % moreOptions['chr'][3:])

    myret = {}
    myretdetail = {}
    if (commonOptions['SplitAndReAlign'] in [0, 2]) or myheader.testall:
        start_time = time.time()
        if commonOptions['outlog'] <= myheader.M_INFO and ('thread' not in specifiedOptions or specifiedOptions['thread']<2):
            print('p2bamhmm start')
        p2bamhmm = getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions)
        memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        if p2bamhmm is None:
            if  (not specifiedOptions.has_key('thread')) or specifiedOptions['thread']<2:
               print('ERROR None detection', moreOptions['repeatName'], mgloc)
               logging.error('ERROR None detection: ' +
                          str(moreOptions['repeatName']) + ' ' + str(mgloc))
        
        addSumForAGene(p2bamhmm, myret, myretdetail, 'p2bamhmm', 2)
        end_time = time.time()
        if commonOptions['outlog'] <= myheader.M_WARNING and ('thread' not in specifiedOptions or specifiedOptions['thread']<2):
            print('p2bamhmm end---running time%.0f mem%d' % (end_time - start_time, memres))
            sys.stdout.flush()
    if ((commonOptions['SplitAndReAlign'] in [1, 2]) or myheader.testall) and (commonOptions['SeqTech'] not in ["Illumina"]):
        from . import myRepeatReAlignment
        start_time = time.time()
        if commonOptions['outlog'] <= myheader.M_INFO and ('thread' not in specifiedOptions or specifiedOptions['thread']<2):
            print('p2sp start')
        moreOptions['fafqfile'] = specifiedOptions["bamfile"]
        moreOptions['fafqtype'] = 'bam'
        p2sp = myRepeatReAlignment.getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
        memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        if p2sp is None:
            if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
               print('ERROR None detection (sp)', moreOptions['repeatName'], mgloc)
               logging.error('ERROR None detection (sp): ' +
                          str(moreOptions['repeatName']) + ' ' + str(mgloc))
        
        addSumForAGene(p2sp, myret, myretdetail, 'p2sp', 2)
        end_time = time.time()
        if commonOptions['outlog'] <= myheader.M_WARNING and ('thread' not in specifiedOptions or specifiedOptions['thread']<2):
            print('p2sp end---running time%.0f mem%d' % (end_time - start_time, memres))
            sys.stdout.flush()

    return [myret, myretdetail]


def getRepeat(commonOptions, specifiedOptions):
    summary = {}

    gLoc = commonOptions['gLoc']
    moreOptions = {}

    glkeys = gLoc.keys()
    glkeys.sort()
    for glk in glkeys:
        moreOptions['repeatName'] = glk
        summary[glk] = getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions)

    return summary
