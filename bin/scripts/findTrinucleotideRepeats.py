
import re;
import os;
import sys;
import string;
import math;
import copy

import numpy;
import peakutils;
import argparse;

import logging

from scipy.stats import norm

import heapq



import getAlignment
import myHMM


from myheader import *

#LOG_FILENAME = 'TrinRepDis.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

def myReadTxtFile(filename):
        f = open(filename,'r')

        data = f.readlines();
        while string.strip(data[-1])=="": data = data[:-1];
        f.close();

        return data;

def myWriteTxtFile(mlist, filename):
        f = open(filename, 'w')

        for ml in mlist:
                if ml[-1]=='\n': f.write(ml);
                else: f.write(ml+'\n');

        f.close();

def getDiseseGeneInRefGenomeLocation_hg19():
	gLoc = {};  plusminus = '-'; plusminus = '+';
	#CSF1PO  5q33.1   TCTA +12       chr5:149355735 - 149556053    chr5:149455883 - 149455959 AGAT(top strand commonly used
	locus_name='CSF1PO'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr5', '149355735','149556053', '149455883','149455959', 'TCTA', plusminus+'12', '']

	#TH01    11p15.5  AATG +9        chr11:2092277 - 2292522       chr11:2192322 - 2192345    AATG bottom strand (commonly used)	
	locus_name='TH01'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr11', '2092277','2292522', '2192322','2192345', 'AATG', plusminus+'9', '']

	#TPOX    2p25.3   AATG +11       chr2:1393368 - 1593481        chr2:1493421 - 1493456     AATG top strand (commonly used
	locus_name='TPOX'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr2', '1393368','1593481', '1493421','1493456', 'AATG', plusminus+'11', '']

	#D3S1358 3p21.31  TCTA +16       chr3:45482205 - 45682335      chr3:45582231 - 45582294   TCTA bottom strand
	locus_name='D3S1358'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr3', '45482205','45682335', '45582231','45582294', 'TCTA', plusminus+'16', '']

	#D5S818  5q23.2   TCTA +11       chr5:123011125 - 123211402    chr5:123111247 - 123111290 AGAT top strand
	locus_name='D5S818'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr5', '123011125','123211402', '123111247','123111290', 'TCTA', plusminus+'11', '']

	#D7S820  7q21.11  TATC +12       chr7:83689381 - 83889718      chr7:83789530 - 83789593   GATA top strand
	locus_name='D7S820'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr7', '83689381','83889718', '83789530','83789593', 'TATC', plusminus+'12', '']

	#D8S1179 8q24.13  TATC +12       chr8:125807064 - 126007405    chr8:125907105 - 125907160 TATC top strand
	locus_name='D8S1179'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr8', '125807064','126007405', '125907105','125907160', 'TATC', plusminus+'12', '']
	
	#D13S317 13q31.1  TATC +12       chr13:82622033 - 82822314     chr13:82722160 - 82722203  GATA bottom strand (commonly used
	locus_name='D13S317'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr13', '82622033','82822314', '82722160','82722203', 'TATC', plusminus+'12', '']
	
	#D16S539 16q24.1  GATA +11       chr16:86286034 - 86486428     chr16:86386308 - 86386351  GATA top strand
	locus_name='D16S539'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr16', '86286034','86486428', '86386308','86386351', 'GATA', plusminus+'11', '']
	
	#D18S51  18q21.33 GAAA +19       chr18:60848678 - 61049364     chr18:60948901 - 60949006  GAAA top strand
	locus_name='D18S51'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr18', '60848678','61049364', '60948901','60949006', 'GAAA', plusminus+'19', '']

	#D21S11  21q21.1  TCTA +27       chr21:20454263 - 20654483     chr21:20554291 - 20554441  TCTG top strand
	locus_name='D21S11'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr21', '20454263','20654483', '20554291','20554441', 'TCTA', plusminus+'27', '']
	
	#D2S441  2p14     TCTA +12       chr2:68139016 - 68339157      chr2:68239079 - 68239126   TCTA top strand
	locus_name='D2S441'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr2', '68139016','68339157', '68239079','68239126', 'TCTA', plusminus+'12', '']

	#D10S1248 10q26.3 GGAA +13       chr10:130992374 - 131192796   chr10:131092508 - 131092559
	locus_name='D10S1248'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr10', '130992374','131192796', '131092508','131092559', 'GGAA', plusminus+'13', '']

	#D22S1045 22q12.3 ATT  +17       chr22:37436285 - 37636570     chr22:37536327 - 37536380  ATT top GenBank strand
	locus_name='D22S1045'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr22', '37436285','37636570', '37536327','37536380', 'ATT', plusminus+'17', '']


	locus_name='D8S1179'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr8','125907063','125907405','125907104','125907159','TATC','+13','',]
	locus_name='D21S11'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr21','20554262','20554483','20554290','20554441','TCTA','+38','',]
	locus_name='D7S820'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr7','83789380','83789718','83789528','83789593','CTAT','+16','',]
	locus_name='CSF1PO'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr5','149455734','149456053','149455884','149455960','CTAT','+19','',]
	locus_name='D3S1358'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr3','45582204','45582335','45582228','45582295','TATC','+16','',]
	locus_name='TH01'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr11','2192276','2192522','2192315','2192346','TGAA','+7','',]
	locus_name='D13S317'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr13','82722032','82722314','82722159','82722223','TATC','+16','',]
	locus_name='D16S539'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr16','86386033','86386428','86386307','86386351','GATA','+11','',]
	locus_name='D19S433'; locus_name = locus_name.lower(); ######
	gLoc[locus_name] = ['chr19','30416989','30417261','30417140','30417206','TCCT','+16','',]
	locus_name='TPOX'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr2','1493367','1493481','1493422','1493456','TGAA','+8','',]
	locus_name='D18S51'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr18','60948677','60949364','60948899','60949006','AGAA','+26','',]
	locus_name='D5S818'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr5','123111124','123111402','123111246','123111293','TCTA','+11','',]
	#locus_name='FGA'; locus_name = locus_name.lower();
	#gLoc[locus_name] = ['FGA','155508847','155509043','155508872','155508975','AAAG','+25','',]
	locus_name='D1S1656'; locus_name = locus_name.lower(); ######
	gLoc[locus_name] = ['chr1','230905195','230905509','230905362','230905429','CTAT','+16','',]
	locus_name='D2S441'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr2','68239015','68239157','68239078','68239127','TCTA','+12','',]
	locus_name='D10S1248'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr10','131092373','131092796','131092507','131092559','GGAA','+13','',]
	locus_name='D22S1045'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr22','37536284','37536570','37536326','37536382','ATT','+18','',]
	locus_name='D20S1082'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr20','53865825','53866083','53865937','53865979','ATA','+14','',]
	locus_name='D6S474'; locus_name = locus_name.lower();  #******
	gLoc[locus_name] = ['chr6','112878949','112879283','112879151','112879231','TAGA','+20','',]
	locus_name='D1S1677'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr1','163559699','163560041','163559815','163559892','TTCC','+19','',]
	locus_name='D11S4463'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr11','130872238','130872721','130872403','130872462','TCTA','+14','',]
	locus_name='D4S2364'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr4','93517335','93517632','93517367','93517421','ATTC','+13','',]
	locus_name='D9S1122'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr9','79688593','79688849','79688733','79688791','TAGA','+14','',]
	locus_name='D2S1776'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr2','169645171','169645775','169645402','169645454','AGAT','+13','',]
	locus_name='D10S1435'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr10','2243219','2243552','2243331','2243392','TATC','+15','',]
	locus_name='D3S3053'; locus_name = locus_name.lower();  #******
	gLoc[locus_name] = ['chr3','171750701','171751224','171750963','171751000','AGAT','+9','',]
	locus_name='D5S2500'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr5','58697040','58697348','58697269','58697313','CTAT','+11','',]
	locus_name='D1S1627'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr1','106963664','106963777','106963713','106963752','ATT','+13','',]
	locus_name='D3S4529'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr3','85852473','85852736','85852632','85852684','AGAT','+13','',]
	locus_name='D17S974'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr17','10518665','10518972','10518722','10518791','AGAT','+17','',]
	locus_name='D6S1017'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr6','41677173','41677509','41677267','41677308','TGGA','+10','',]
	locus_name='D4S2408'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr4','31304230','31304513','31304419','31304456','ATCT','+9','',]
	locus_name='D9S2157'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr9','136035466','136035859','136035668','136035699','ATA','+10','',]
	locus_name='D17S1301'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr17','72680785','72681109','72680993','72681041','AGAT','+12','',]
	locus_name='D18S853'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr18','3990523','3990853','3990628','3990665','ATA','+12','',]
	locus_name='D20S482'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr20','4506247','4506466','4506337','4506396','AGAT','+14','',]
	locus_name='D14S1434'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr14','95308122','95308685','95308401','95308443','TCTA','+10','',]
	#locus_name='DYS390'; locus_name = locus_name.lower(); ######
	#gLoc[locus_name] = ['chrY','17274873','17275219','17274932','17275050','AGAT','+29','',]
	#locus_name='DYS19'; locus_name = locus_name.lower(); ######
	#gLoc[locus_name] = ['chrY','9521933','9522128','9521986','9522052','TATC','+16','',]
	#locus_name='DYS393'; locus_name = locus_name.lower();
	#gLoc[locus_name] = ['chrY','88861983','88862275','88862009','88862064','ATAG','+13','',]
	#locus_name='DYS391'; locus_name = locus_name.lower();
	#gLoc[locus_name] = ['chrY','14102739','14103068','14102780','14102841','TATC','+15','',]
	#locus_name='DYS392'; locus_name = locus_name.lower();
	#gLoc[locus_name] = ['chrY','22633713','22634011','22633871','22633912','AAT','+13','',]
	locus_name='D5S2505'; locus_name = locus_name.lower();  #******
	gLoc[locus_name] = ['chr5','5816996','5817453','5817134','5817338','AGAT','+53','',]
	locus_name='D8S1115'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr8','42536554','42536839','42536588','42536615','AAT','+9','',]
	locus_name='D10S2325'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr10','12792945','12793245','12793050','12793126','ATAAG','+15','',]
	locus_name='D11S554'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr11','44930233','44930472','44930285','44930381','AAAG','+23','',]
	locus_name='D13S308'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr13','26457002','26457149','26457032','26457111','GAT','+26','',]
	locus_name='D14S306'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr14','38328082','38328451','38328290','38328392','TCTA','+25','',]
	locus_name='D18S535'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr18','38148774','38149329','38148825','38148894','TAGA','+17','',]
	locus_name='D18S849'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr18','55564577','55564888','55564619','55564758','TCTA','+35','',]
	locus_name='D19S253'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr19','15728135','15728364','15728294','15728345','ATCT','+12','',]
	locus_name='D1S103'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr1','230836735','230836986','230836894','230836933','TG','+19','',]
	locus_name='D20S85'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr20','38051608','38051920','38051753','38051830','TTTC','+20','',]
	locus_name='D21S2055'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr21','41191368','41191662','41191434','41191586','CTAT','+39','',]
	locus_name='D2S1360'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr2','17491921','17492295','17491984','17492069','AGAT','+21','',]
	locus_name='D2S410'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr2','116240906','116241220','116241001','116241056','TCTA','+13','',]
	locus_name='D2S436'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr2','107243132','107243450','107243293','107243381','TCTA','+22','',]
	locus_name='D3S1359'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr3','49948548','49948757','49948611','49948679','AGAT','+16','',]
	locus_name='D3S1744'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr3','147092494','147092853','147092538','147092646','ATAG','+27','',]
	locus_name='D4S2366'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr4','6484677','6484947','6484828','6484901','GATA','+18','',]
	locus_name='D5S815'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr5','90990446','90990802','90990491','90990590','TATC','+24','',]
	locus_name='D6S1043'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr6','92449866','92450188','92449940','92450005','CTAT','+16','',]
	locus_name='D6S502'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr6','125907063','125907405','125907104','125907159','TATC','+13','',]
	locus_name='D7S1517'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr7','123497650','123498050','123497698','123497787','CTTT','+22','',]
	locus_name='D8S1132'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr8','107328820','107329072','107328919','107329003','TCTA','+21','',]
	locus_name='D8S344'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr8','131076949','131077362','131077139','131077180','TG','+20','',]
	locus_name='D8S347'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr8','129340452','129340831','129340529','129340619','TCTA','+22','',]
	locus_name='D8S639'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr8','16770880','16771227','16771051','16771172','ATAG','+30','',]
	#locus_name='DXS7132'; locus_name = locus_name.lower(); #******
	#gLoc[locus_name] = ['chrX','64655267','64655632','64655524','64655595','GATA','+17','',]
	#locus_name='DXS7423'; locus_name = locus_name.lower(); #******
	#gLoc[locus_name] = ['chrX','149710902','149711089','149710970','149711041','TGGA','+17','',]
	#locus_name='DXS8378'; locus_name = locus_name.lower(); #******
	#gLoc[locus_name] = ['chrX','9370149','9370429','9370301','9370341','ATAG','+10','',]
	locus_name='HUMTH01'; locus_name = locus_name.lower(); #******
	gLoc[locus_name] = ['chr11','2192276','2192522','2192315','2192346','TGAA','+7','',]

	#SE33                    6q14            AAAG   chr6:88986839-88987087   TCTT   chr6:88986849-88987076 '26'
	locus_name='SE33'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr6', '88986839', '88987087', '88986849', '88987076', 'TCTT', '+26', '']
	#Penta_D                 21q22.3         AAAGA  chr21:45055834-45056398  AAAGA chr21:45056073-45056153  '13'
	locus_name='Penta_D'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr21', '45055834', '45056398', '45056073', '45056153','AAAGA',  '+13', '']
	#Penta_E                 15q26.2         AAAGA chr15:97374212-97374565 TTTTC chr15:97374242-97374269    '5'
	locus_name='Penta_E'; locus_name = locus_name.lower();
        gLoc[locus_name] = ['chr15', '97374212', '97374565', '97374242', '97374269', 'TTTTC', '+5', '']
	#F13B                    1q31-q32.1      AAAT  chr1:197007783-197007951 TAAA chr1:197007832-197007869 '10'
	locus_name='F13B'; locus_name = locus_name.lower();
        gLoc[locus_name] = ['chr1', '197007783', '197007951', '197007832', '197007869', 'TAAA', '+10', '']

	glockeys = gLoc.keys(); glockeys.sort()
	for glk in glockeys:
		genstr = int(gLoc[glk][1]); genend = int(gLoc[glk][2])
		repstr = int(gLoc[glk][3]); repend = int(gLoc[glk][4])
		if repstr-genstr<0 or genend-repend<0:
			print ('Wrong gene and repeat information: gene = [%d, %d], and repeat = [%d, %d]' % genstr,genend,repstr,repend)

	return gLoc;

def getDiseseGeneInRefGenomeLocation(hg):
	if hg=='' or hg=='hg38': return getDiseseGeneInRefGenomeLocation_hg38();
	elif hg=='hg19': return getDiseseGeneInRefGenomeLocation_hg19();
	else:
		logging.error('Error: not supported yet')

#from wiki: https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
#RefSeq chr:start-end would be used.	
def getDiseseGeneInRefGenomeLocation_hg38():
        gLoc = {}
        #HTT at chr4:3074510-3243960
        gLoc['htt'] = ['chr4','3074510','3243960','3074877','3074939', 'CAG', '+19', '6-35:36-250']                # (CAG)* +1CAACAG 
        #ATN1 at         chr12:   6924463 -  6942321   6936729   6936773
        gLoc['atn1'] = ['chr12', '6924463', '6942321','6936717', '6936773', 'CAG', '+10', '6-35:49-88']   # 2(CAGCAA) (CAG)*
        #ar: (Q)* <+ 5non-Q + 7 Q>
        gLoc['ar'] = ['chrX', '67544032', '67730619', '67545318', '67545386', 'CAG', '+22', '9-36:38-62']         # (CAG)* +1CAA
        gLoc['atxn1'] = ['chr6', '16299112', '16761490', '16327636', '16327722', 'CAG', '-14', '6-35:49-88']      # (CAG)* +2(ATGCTG) 11(CTG)
        gLoc['atxn2'] = ['chr12', '111452214', '111599676', '111598951', '111599019', 'CAG', '-13', '14-32:33-77'] # 10Q  (CAG)*
        #atxn3: chr:start-end for RefSeq is chr14:92058552-92106621. chr:start-end for UCSC is chr14:92038652-92106610
        #ATXN3 at chr14:92058552-92106621 
        gLoc['atxn3'] = ['chr14', '92038652', '92106610', '92071011', '92071052', 'CAG', '-14', '12-40:55-86']       # (CAG)* + CTGTTGCTGCTTTTGCTGCTG
        #CACNA1A at chr19:13206443-13506460 
        gLoc['cacna1a'] = ['chr19', '13206443', '13506460', '13207859', '13207897', 'CAG', '-8', '4-18:21-30']     #y
        #ATXN7 at chr3:63864557-64003460 
        #could not be found using UCSC genome browser with "simple repeat" on
        gLoc['atxn7'] = ['chr3','63864557', '64003460', '63912686', '63912715', 'CAG', '+7', '7-17:38-120']        # (CAG)* + 3CAG
        #TBP at chr6:           170554333 -  170572870
        gLoc['tbp'] = ['chr6', '170554333', '170572870', '170561899', '170562021', 'CAG', '+19', '25-42:47-63']    # 20*3 (CAG)* +1CAACAG
        #FMR1 at chrX:147911951-147951127 
        gLoc['fmr1'] = ['chrX', '147911951', '147951127', '147912051', '147912110', 'CGG', '+10', '6-53:230+/55-200']     # (CGG)* + 1AGG + 9(CGG)
        #aff2: 148500612--148500639: 6 CCG + 3 non-CCG
        #gLoc['aff2'] = ['chrX', '148500619', '149000663', '148500639', '148500692', 'CCG', '+15', '6-35:200+']   # (CCG)* + 1CTG + 2CCG
        gLoc['aff2'] = ['chrX', '148500619', '149000663', '148500606', '148500692', 'CCG', '+15', '6-35:200+']   
        #chr9:69,037,262-69,037,374    FXN at chr9:69035563-69100178
        #could not be found using UCSC genome browser with "simple repeat" on
        gLoc['fxn'] = ['chr9', '69035563', '69100178', '69037287', '69037305', 'GAA', '+5', '7-34:100+']         #xxxy
        #RefSeq for dmpk: chr19:45769709-45782557 
        gLoc['dmpk'] = ['chr19', '45769709', '45782557', '45770205', '45770264', 'CTG', '-20', '5-37:50+']      #y
        #gLoc['sca8'] = ['chr13', '70681345', '70713561', '70139384', '70139428', 'CTG', '15']      #XXXXy
        #gLoc['sca8'] = ['chr13', '70107213', '70139552', '70139351', '70139428', 'CTG', '15']      #Xy  1TTA + 10 CTA + (CTG)*
        #ATXN8OS at       chr13 :  70107213 -  70139753
        gLoc['atxn8os'] = ['chr13', '70107213', '70139753', '70139351', '70139428', 'CTG', '+15', '16-37:110-250']   # 1TTA + 10 CTA + (CTG)*
       
	#https://en.wikipedia.org/wiki/Spinocerebellar_ataxia 
        #PPP2R2B at chr5:           146589505 -  147081520 RefSeq
        gLoc['ppp2r2b'] = ['chr5', '146589505', '147081520', '146878729', '146878759', 'CAG', '-10', '7-28:66-78'] #may contain errors;

        locus_name='D8S1179'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','124894821','124895163','124894862','124894917','TATC','+13','',]
        locus_name='D21S11'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr21','19181944','19182165','19181972','19182123','TCTA','+38','',]
        locus_name='D7S820'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr7','84160064','84160402','84160212','84160277','CTAT','+16','',]
        locus_name='CSF1PO'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','150076171','150076490','150076321','150076397','CTAT','+19','',]
        locus_name='D3S1358'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr3','45540712','45540843','45540736','45540803','TATC','+16','',]
        locus_name='TH01'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','2171046','2171292','2171085','2171116','TGAA','+7','',]
        locus_name='D13S317'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr13','82147897','82148179','82148024','82148088','TATC','+16','',]
        locus_name='D16S539'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr16','86352427','86352822','86352701','86352745','GATA','+11','',]
        locus_name='D19S433'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr19','29926082','29926354','29926233','29926299','TCCT','+16','',]
        locus_name='D18S51'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','63281444','63282131','63281666','63281773','AGAA','+26','',]
        locus_name='D5S818'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','123775430','123775708','123775552','123775599','TCTA','+11','',]
        #locus_name='FGA'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chr4','154587695','154587891','154587720','154587823','AAAG','+25','',]
        locus_name='D1S1656'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','230769449','230769763','230769616','230769683','CTAT','+16','',]
        locus_name='D2S441'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','68011883','68012025','68011946','68011995','TCTA','+12','',]
        locus_name='D10S1248'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','129294109','129294532','129294243','129294295','GGAA','+13','',]
        locus_name='D22S1045'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr22','37140244','37140530','37140286','37140342','ATT','+18','',]
        locus_name='D20S1082'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr20','55249286','55249544','55249398','55249440','ATA','+14','',]
        locus_name='D6S474'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr6','112557747','112558081','112557949','112558029','TAGA','+20','',]
        locus_name='D1S1677'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','163589909','163590251','163590025','163590102','TTCC','+19','',]
        locus_name='D11S4463'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','131002343','131002826','131002508','131002567','TCTA','+14','',]
        locus_name='D4S2364'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr4','92596184','92596481','92596216','92596270','ATTC','+13','',]
        locus_name='D9S1122'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr9','77073677','77073933','77073817','77073875','TAGA','+14','',]
        locus_name='D2S1776'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','168788661','168789265','168788892','168788944','AGAT','+13','',]
        locus_name='D10S1435'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','2201025','2201358','2201137','2201198','TATC','+15','',]
        locus_name='D3S3053'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr3','172032911','172033434','172033173','172033210','AGAT','+9','',]
        locus_name='D5S2500'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','59401214','59401522','59401443','59401487','CTAT','+11','',]
        locus_name='D1S1627'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','106421042','106421155','106421091','106421130','ATT','+13','',]
        locus_name='D3S4529'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr3','85803323','85803586','85803482','85803534','AGAT','+13','',]
        locus_name='D17S974'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr17','10615348','10615655','10615405','10615474','AGAT','+17','',]
        locus_name='D6S1017'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr6','41709435','41709771','41709529','41709570','TGGA','+10','',]
        locus_name='D4S2408'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr4','31302608','31302891','31302797','31302834','ATCT','+9','',]
        locus_name='D9S2157'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr9','133160079','133160472','133160281','133160312','ATA','+10','',]
        locus_name='D17S1301'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr17','74684646','74684970','74684854','74684902','AGAT','+12','',]
        locus_name='D18S853'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','3990523','3990853','3990628','3990665','ATA','+12','',]
        locus_name='D20S482'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr20','4525601','4525820','4525691','4525750','AGAT','+14','',]
        locus_name='D14S1434'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr14','94841785','94842348','94842064','94842106','TCTA','+10','',]
        #locus_name='DYS390'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','15162993','15163339','15163052','15163170','AGAT','+29','',]
        #locus_name='DYS19'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','9684324','9684519','9684377','9684443','TATC','+16','',]
        #locus_name='DYS393'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrX','89606984','89607276','89607010','89607065','ATAG','+13','',]
        #locus_name='DYS391'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','11982033','11982362','11982074','11982135','TATC','+15','',]
        #locus_name='DYS392'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','20471827','20472125','20471985','20472026','AAT','+13','',]
        locus_name='D5S2505'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','5816883','5817340','5817021','5817225','AGAT','+53','',]
        locus_name='D8S1115'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','42681411','42681696','42681445','42681472','AAT','+9','',]
        locus_name='D10S2325'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','12750946','12751246','12751051','12751127','ATAAG','+15','',]
        locus_name='D11S554'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','44908682','44908921','44908734','44908830','AAAG','+23','',]
        locus_name='D14S306'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr14','37858877','37859246','37859085','37859187','TCTA','+25','',]
        locus_name='D18S535'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','40568810','40569365','40568861','40568930','TAGA','+17','',]
        locus_name='D18S849'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','57897345','57897656','57897387','57897526','TCTA','+35','',]
        locus_name='D19S253'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr19','15617324','15617553','15617483','15617534','ATCT','+12','',]
        locus_name='D1S103'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','230700989','230701240','230701148','230701187','TG','+19','',]
        locus_name='D20S85'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr20','39422965','39423277','39423110','39423187','TTTC','+20','',]
        locus_name='D21S2055'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr21','39819441','39819735','39819507','39819659','CTAT','+39','',]
        locus_name='D2S1360'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','17310654','17311028','17310717','17310802','AGAT','+21','',]
        locus_name='D2S410'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','115483330','115483644','115483425','115483480','TCTA','+13','',]
        locus_name='D2S436'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','106626676','106626994','106626837','106626925','TCTA','+22','',]
        locus_name='D3S1744'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr3','147374707','147375066','147374751','147374859','ATAG','+27','',]
        locus_name='D4S2366'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr4','6482950','6483220','6483101','6483174','GATA','+18','',]
        locus_name='D5S815'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','91694629','91694985','91694674','91694773','TATC','+24','',]
        locus_name='D6S1043'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr6','91740148','91740470','91740222','91740287','CTAT','+16','',]
        locus_name='D6S502'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','124894821','124895163','124894862','124894917','TATC','+13','',]
        locus_name='D7S1517'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr7','123857596','123857996','123857644','123857733','CTTT','+22','',]
        locus_name='D8S1132'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','106316592','106316844','106316691','106316775','TCTA','+21','',]
        locus_name='D8S344'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','130064703','130065116','130064893','130064934','TG','+20','',]
        locus_name='D8S347'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','128328206','128328585','128328283','128328373','TCTA','+22','',]
        locus_name='D8S639'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','16913371','16913718','16913542','16913663','ATAG','+30','',]
        #locus_name='DXS7132'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrX','65435387','65435752','65435644','65435715','GATA','+17','',]
        #locus_name='DXS7423'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrX','150542453','150542640','150542521','150542592','TGGA','+17','',]
        #locus_name='DXS8378'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrX','9402109','9402389','9402261','9402301','ATAG','+10','',]
        locus_name='HUMTH01'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','2171046','2171292','2171085','2171116','TGAA','+7','',]
        locus_name='CSF1PO'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','150076171','150076490','150076321','150076397','CTAT','+19','',]
        #locus_name='FGA'; locus_name = locus_name.lower();   
        #gLoc[locus_name] = ['chr4','154587695','154587891','154587720','154587823','AAAG','+25','',]
        locus_name='TH01'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','2171046','2171292','2171085','2171116','TGAA','+7','',]
        locus_name='D1S1656'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','230769449','230769763','230769616','230769683','CTAT','+16','',]
        locus_name='D2S441'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','68011883','68012025','68011946','68011995','TCTA','+12','',]
        locus_name='D3S1358'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr3','45540712','45540843','45540736','45540803','TATC','+16','',]
        locus_name='D5S818'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr5','123775430','123775708','123775552','123775599','TCTA','+11','',]
        locus_name='D6S1043'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr6','91740148','91740470','91740222','91740287','CTAT','+16','',]
        locus_name='D7S820'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr7','84160064','84160402','84160212','84160277','CTAT','+16','',]
        locus_name='D8S1179'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr8','124894821','124895163','124894862','124894917','TATC','+13','',]
        locus_name='D10S1248'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','129294109','129294532','129294243','129294295','GGAA','+13','',]
        locus_name='D13S317'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr13','82147897','82148179','82148024','82148088','TATC','+16','',]
        locus_name='D14S1434'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr14','94841785','94842348','94842064','94842106','TCTA','+10','',]
        locus_name='D16S539'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr16','86352427','86352822','86352701','86352745','GATA','+11','',]
        locus_name='D18S51'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','63281444','63282131','63281666','63281773','AGAA','+26','',]
        locus_name='D19S433'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr19','29926082','29926354','29926233','29926299','TCCT','+16','',]
        locus_name='D21S11'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr21','19181944','19182165','19181972','19182123','TCTA','+38','',]
        locus_name='D22S1045'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr22','37140244','37140530','37140286','37140342','ATT','+18','',]
        #locus_name='DYS19'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','9684324','9684519','9684377','9684443','TATC','+16','',]
        #locus_name='DYS388'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','12635276','12635761','12635603','12635639','AAT','+12','',]
        #locus_name='DYS390'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','15162993','15163339','15163052','15163170','AGAT','+29','',]
        #locus_name='DYS391'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','11982033','11982362','11982074','11982135','TATC','+15','',]
        #locus_name='DYS392'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','20471827','20472125','20471985','20472026','AAT','+13','',]
        #locus_name='DYS393'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrX','89606984','89607276','89607010','89607065','ATAG','+13','',]
        #locus_name='Y-GATA-A4'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','12403374','12403649','12403460','12403566','GATA','+26','',]
        #locus_name='Y-GATA-A7.1'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','18888689','18889221','18888801','18888995','TATC','+49','',]
        #locus_name='Y-GATA-A7.2'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','18888689','18889221','18888801','18888995','TATC','+49','',]
        #locus_name='Y-GATA-A10'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','16606889','16607107','16607006','16607074','TATC','+17','',]
        #locus_name='Y-GATA-H4'; locus_name = locus_name.lower();
        #gLoc[locus_name] = ['chrY','16631292','16631877','16631623','16631759','CTAT','+33','',]

        locus_name='D2S1338'; locus_name = locus_name.lower(); ########
        gLoc[locus_name] = ['chr2','218014645','218014994','218014853','218014950','AGGA','+17','',]
        locus_name='TPOX'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr2','1489595','1489709','1489650','1489684','TGAA','+8','',]
        locus_name='D12S391'; locus_name = locus_name.lower(); ########
        gLoc[locus_name] = ['chr12','12296939','12297292','12297013','12297094','ATAG','+13','',]
        locus_name='D8S320'; locus_name = locus_name.lower();
        gLoc[locus_name] = ['chr8','88818683','88819178','88818857','88819107','TTTC','+9','',]

        return gLoc;

def get_gLoc(dis_type, gLoc):
        dis_type = dis_type.lower();
        return gLoc[dis_type];

def getValues(mdict, minreads):
        valdict = {}
        mdkeys = mdict.keys();
        for mdk in mdkeys:
                if not valdict.has_key(mdict[mdk]):
                        valdict[mdict[mdk]] = 0;
                valdict[mdict[mdk]] += 1;
        valkeys = valdict.keys(); valkeys.sort();
        more = []
        #for i in range(1, len(valkeys)*2/3+1):
        #print 'test', valkeys, minreads
        for i in range(1, len(valkeys)+1):
                #if valkeys[-i]<4: continue;
                #
                #revised on Dec 12, 2016
                #change from above to the blow
                if valkeys[-i]<minreads: continue;
                #print -i, valkeys[-i], valdict[valkeys[-i]]
                if valdict[valkeys[-i]]>=2: more.append(valkeys[-i])
        #print 'test', more
        return more;

def findSameValues(mdict, v, mdkeys):
        start = -1; end = -1;
        for mdk_ind in range(len(mdkeys)-1):
                if mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]] and mdict[mdkeys[mdk_ind]]==v:
                        if start==-1: start = mdk_ind
                        end = mdk_ind+1
                if (not mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]]) and (not start==-1): break;
        return [start, end]

def reviseDict(mdict, v):
        has_revised = False;
        mdkeys = mdict.keys();  mdkeys.sort()
        while True:
                [start, end] = findSameValues(mdict, v, mdkeys)
		#print 'revisev', v, start, end, mdict
                if start==-1 or end==-1: break;
                else:
                        i=start; j = end;
                        while True:
                                #print i, j, mdkeys[i], mdkeys[j], mdict[mdkeys[i]], mdict[mdkeys[j]]
                                if i>=j:
                                        print 'Error i is larger than j', i, j, start, end;
                                        break;
                             
                                has_revised = True;  
                                i += 1;
                                mdict[mdkeys[i]] += 1;
                                if i==j: break;
                                j -= 1;
                                if i==j: break;
                                mdict[mdkeys[j]] += 1;
                                ##print 'ij', i, j, mdkeys[i], mdkeys[j], mdict[mdkeys[i]], mdict[mdkeys[j]]

        return [mdict, has_revised];

def reviseDictAccordingV(mdict, minreads):
	while True:
           more = getValues(mdict, minreads)
           #print 'more', more, minreads
           more.sort();
           has_revised = False;
           for m in more:
                   mdict, has1_revised = reviseDict(mdict, m);
                   if not has_revised: has_revised = has1_revised
           if not has_revised: break;
        return mdict


def selectFromTwoX(x_index2a, x_index2b, lendict, x_index1):
	if x_index2a==x_index2b: x_index2 = x_index2b
	else:
		neighbor3 = [0,0]; point3 = [x_index2a, x_index2b]
		for pi in range(len(point3)):
			for xi in range(point3[pi]-1, point3[pi]+2):
				if lendict.has_key(xi): neighbor3[pi] += lendict[xi]
		
		#print x_index2a, neighbor3[0], x_index2b, neighbor3[1]
		
		if neighbor3[0]*x_index2a < neighbor3[1]*x_index2b:
			x_index2 = x_index2b
		elif neighbor3[0]*x_index2a > neighbor3[1]*x_index2b:
			x_index2 = x_index2a
		elif x_index2a*2 < x_index2b and neighbor3[0]*0.9 < neighbor3[1]:
			x_index2 = x_index2b
		elif x_index2b*2 < x_index2a and neighbor3[1]*0.9 < neighbor3[0]:
			x_index2 = x_index2a
		else:
			if neighbor3[0]>neighbor3[1]: x_index2 = x_index2a
			elif neighbor3[0]<neighbor3[1]: x_index2 = x_index2b
			else:
				if x_index2a>x_index2b: x_index2 = x_index2a;
				elif x_index2a<x_index2b: x_index2 = x_index2b
				else:
					if math.fabs(x_index1-x_index2a)>math.fabs(x_index1-x_index2b): x_index2 = x_index2a;
					else: x_index2 = x_index2b
	return x_index2

def selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x):
	if y_index2a==y_index2b: y_index2 = y_index2b
	else:
		x2 = selectFromTwoX(x[y_index2a], x[y_index2b], lendict, x[y_index1])
		if x2==x[y_index2a]: y_index2 = y_index2a
		else: y_index2 = y_index2b
	return y_index2
	

def selectOne(iy, curyvalue, lendict, y_index1, indexes, x):
	y_index2a = indexes[iy.index(curyvalue)];
	iy.reverse();
	y_index2b = indexes[len(iy) - iy.index(curyvalue) - 1]

	if y_index2a==y_index1: return y_index2b
	elif y_index2b==y_index1: return y_index2a
	
	return selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x)


#def getPeaks2(x, y, lendict, mm, mdebug):
#revised frmo above to below on Dec 12 2016.
def getPeaks2(x, y, lendict, mm, mdebug, minreads, pnearby_size=3):
	indexes = peakutils.indexes(numpy.array(y), thres=0.01) - mm
	peak2 = [];
	pnearby = []; 
	
	if mdebug:
		for i in range(len(x)): print ('%d=%d;' % (x[i], y[i])),
		print 'mindexes',
		for i in indexes: print x[i],
		print ''
		print '#x =',
		for i in range(len(x)): print ('%3d' % x[i]),
		print ''
		print 'y  =',
		for i in range(len(x)): print ('%3d' % y[i]),
		print ''

	if len(indexes)>1:
		iy = []
		for i in indexes: iy.append(y[i]);
		ylargest = heapq.nlargest(2, iy);
		y_index1 = indexes[iy.index(ylargest[0])];
		
		p_y_index2 = selectOne(copy.deepcopy(iy), ylargest[1], lendict, y_index1, indexes, x)
		if len(iy)>2:
			ylargest = heapq.nlargest(3, iy);
			p_y_index3 = selectOne(iy, ylargest[2], lendict, y_index1, indexes, x)
			y_index2 = selectFromTwo(p_y_index2, p_y_index3, lendict, y_index1, x)
		else: y_index2 = p_y_index2


		#ylargest = heapq.nlargest(2, iy);
		#y_index1 = indexes[iy.index(ylargest[0])];
		#
		#y_index2a = indexes[iy.index(ylargest[1])];
		#iy.reverse();
		#y_index2b = indexes[len(iy) - iy.index(ylargest[1]) - 1]
		#
		#if y_index2a==y_index2b: y_index2 = y_index2b
		#else:
		#	neighbor3 = [0,0]; point3 = [y_index2a, y_index2b]
		#	for pi in range(len(point3)):
		#		for xi in range(point3[pi]-1, point3[pi]+2):
		#			if lendict.has_key(xi):	neighbor3[pi] += lendict[xi]
		#	if neighbor3[0]>neighbor3[1]: y_index2 = y_index2a
		#	elif neighbor3[0]<neighbor3[1]: y_index2 = y_index2b
		#	else:
		#		if y_index2a>y_index2b: y_index2 = y_index2a;
		#		elif y_index2a<y_index2b: y_index2 = y_index2b
		#		else:
		#			if math.fabs(y_index1-y_index2a)>math.fabs(y_index1-y_index2b): y_index2 = y_index2a;
		#			else: y_index2 = y_index2b

		peak2.append(x[y_index1]); 
		for ni in range(x[y_index1]-pnearby_size+1, x[y_index1]+pnearby_size):
			pnearby.append(ni);
		#if y_index1-1>0: pnearby.append(x[y_index1-1]);
		#if y_index1+1<len(x): pnearby.append(x[y_index1+1]);
		peak2.append(x[y_index2])
		for ni in range(x[y_index2]-pnearby_size+1, x[y_index2]+pnearby_size):
			pnearby.append(ni);
		#if y_index2-1>0: pnearby.append(x[y_index2-1]); 
		#if y_index2+1<len(x): pnearby.append(x[y_index2+1]);
        else:
		if len(indexes)==1:
			logging.debug("Only find one peak");
			peak2.append(x[indexes[0]])
			for ni in range(x[indexes[0]]-pnearby_size+1, x[indexes[0]]+pnearby_size):
				pnearby.append(ni);
			#pnearby.append(x[indexes[0]]);
			#if indexes[0]-1>0: pnearby.append(x[indexes[0]-1]); 
			#if indexes[0]+1<len(x): pnearby.append(x[indexes[0]+1]);
		else:
			logging.debug("Cannot find a peak");

	#print 'for testing'
	#print x
	#print y;
	#print peak2

	if len(peak2)>1 and y[peak2.index(peak2[1])] > [y[peak2.index(peak2[0])]]:
		peak2[0], peak2[1] = peak2[1], peak2[0]
	
	p2len = len(peak2)
	twotails = [];
	if (x[0] not in pnearby): twotails.append(x[0]);
	if (x[len(x)-1] not in pnearby): twotails.append(x[len(x)-1]);
        #if (x[0] not in peak2) and (x[1] not in peak2):
        #       if lendict[x[0]]==lendict[x[1]]: twotails.append(x[0]);
        #       else: twotails.append(x[0]);twotails.append(x[1]);
        #if (x[len(x)-2] not in peak2) and (x[len(x)-1] not in peak2):
        #       if lendict[x[len(x)-2]]==lendict[x[len(x)-1]]:  twotails.append(x[len(x)-1]);
        #       else: twotails.append(x[len(x)-2]); twotails.append(x[len(x)-1]);
	
	#print twotails
	tailschoosen = []
	for ti in twotails:
		tinearby = [ti-1, ti+1]; shouldchoose = True;
		if lendict[ti]<minreads: shouldchoose=False
		else:
			for tnb in tinearby:
				if lendict.has_key(tnb) and lendict[tnb]>lendict[ti]: shouldchoose=False;
		tailschoosen.append(shouldchoose);
		if not shouldchoose: continue;
		
		for i in range(p2len):
			if selectFromTwoX(ti, peak2[i], lendict, -1)==ti:
			#if lendict[ti] > lendict[peak2[i]]: 
				peak2.insert(i, ti); break;

	if len(peak2)==0:
		for curtt_ind in range(len(twotails)):
			if not tailschoosen[curtt_ind]: continue;
			curtt = twotails[curtt_ind]
			if curtt>4 and lendict[curtt] >= minreads: peak2.append(curtt);
	#print peak2
	if len(peak2)>2: peak2 = peak2[:2];

	return peak2	

def getLargerOverSmallRatio(p1, p2):
	mp = [p1,p2]; mp.sort();
	
	mratio = 2**(-mp[1]/float(mp[0])) #math.exp(-mp[1]/float(mp[0]));
	mratio = math.exp(-mp[1]/float(mp[0]));
	if mratio<0.1: mratio=0.1;
	
	return mratio

def checkSmallSupport(peak2, lendict, minreads):
	newpeak2 = []
	for p in peak2:
		if lendict[p]<minreads: 
			print 'checkSmallSupport', p, lendict[p], minreads, ' for peak2 ', peak2
		else: newpeak2.append(p)
	return newpeak2

def get2Peaks(lengd, SeqDepth):
#def get2Peaks(lengd):
	#prnucleotideRepeats.pyint lengd

	minreads = int(SeqDepth)/20;
	if minreads<0: minreads = 2
	if minreads>5: minreads = 5

	pvalue = 0.05; ratio = 0.4;
	smalleratio = 0.5; #0.55 #0.4; #0.55 # 0.5 #0.6;

	mdebug = True; mdebug = False;
	
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;

	ldkeys = lendict.keys(); ldkeys.sort();
	
	allocr = '';
	for ldk in ldkeys:
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))
	if mdebug: print allocr
	logging.info(allocr)
	#print 'allocr:', allocr

	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: return [[ldkeys[0]], allocr]
	elif len(ldkeys)==2 or len(ldkeys)==3:
		maxk = ldkeys[0]
		if maxk==0: maxk = ldkeys[1]
		for ldk in ldkeys:
			if lendict[ldk]>=lendict[maxk] and ldk>0: maxk = ldk
		return [[maxk], allocr]
		

	len3dict = {}
	for i in range(len(ldkeys)):
                curlk = ldkeys[i]
                if i==0:
                        if curlk==0: len3dict[curlk] = 0
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
                elif i==len(ldkeys)-1:
                        if curlk-1==0: len3dict[curlk] = lendict[curlk]
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk-1): len3dict[curlk] += lendict[curlk-1]
                                else: len3dict[curlk-1] = 0
                else:
                        if curlk-1==0:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk-1): len3dict[curlk] += lendict[curlk-1]
                                else: len3dict[curlk-1] = 0
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
	'''
	for i in range(len(ldkeys)):
                if i==0:
	                if ldkeys[i]==0: len3dict[ldkeys[i]] = 0; #lendict[ldkeys[i]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
                elif i==len(ldkeys)-1:
                        if ldkeys[i-1]==0: len3dict[ldkeys[i]] = lendict[ldkeys[i]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i-1]]
                else:
                        if ldkeys[i-1]==0: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i-1]] + lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
	'''
	
	len3dict = reviseDictAccordingV(len3dict, minreads)
	
	x = []; yo = [];
	for ldk in ldkeys:
		x.append(ldk); 
		yo.append(len3dict[ldk]);

	peak2 = getPeaks2(x, yo, lendict, 0, mdebug, minreads)

	'''
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;

	lendict = reviseDictAccordingV(lendict)
	#print lendict

	x = []; yo = []; allocr = ''; 
	yo = [0]; 
	ldkeys = lendict.keys(); ldkeys.sort();
	for ldk in ldkeys:
		x.append(ldk); yo.append(lendict[ldk]);
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))

	#print lendict	
	logging.info(allocr)
	
	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: return [[yo[1]], allocr]

	yo.append(0); 
	
	pea1k2 = getPeaks2(x, yo[1:-1], lendict, 0, mdebug)	
	
	y2a = []
	for i in range(len(yo)-1): y2a.append(yo[i]+yo[i+1])
	pea2k2 = getPeaks2(x, y2a, lendict, 1, mdebug)

	pea1knearby = []; 
	for p1 in pea1k2:
		for ni in range(p1-2, p1+3):
			pea1knearby.append(ni);
		#pea1knearby.append(p1);
		#p1ind = x.index(p1);
		#if p1ind-1>0: pea1knearby.append(x[p1ind-1]);
		#if p1ind+1<len(x): pea1knearby.append(x[p1ind+1]);

	pea2k2new = []
	for p2 in pea2k2:
		p2ind = x.index(p2);
		nb3 = [p2]
		if p2ind-1>0: nb3.append(x[p2ind-1]);
		if p2ind+1<len(x): nb3.append(x[p2ind+1]);
		curmx = nb3[-1]  #p2;
		for curi in nb3: 
			if lendict[curi]>lendict[curmx]: curmx = curi
		pea2k2new.append(curmx)

	notinpea1k = []
	for p2 in pea2k2new:
		if p2 not in pea1knearby: notinpea1k.append(p2);

	#print pea1k2, notinpea1k, pea2k2new 

	#if len(pea1k2)==0: pea1k2 = pea2k2new
	#else:
	for p2 in notinpea1k:
			has_insert = False;
			for i in range(len(pea1k2)):
				x_index_find = selectFromTwoX(p2, pea1k2[i], lendict, -1)
				if x_index_find == p2:
					pea1k2.insert(i, p2); has_insert = True; break;
				#if lendict[p2] > lendict[pea1k2[i]]:
				#	 pea1k2.insert(i, p2); has_insert = True; break;
			if not has_insert: pea1k2.append(p2)

	peak2 = pea1k2
	'''
	if mdebug: print peak2

	samevalues = []
	for curp_ind in range(len(peak2)):
                curp = peak2[curp_ind]
                pin = ldkeys.index(curp)
                nb3 = [curp]
                if pin-1>0: nb3.append(ldkeys[pin-1])
                if pin+1<len(ldkeys): nb3.append(ldkeys[pin+1])
                curmx = nb3[-1]
                for curi in nb3:
                        #if lendict[curi]>lendict[curmx]: curmx = curi
                        #
                        #revised from the line above to below:
                        if (lendict[curi]*curi>lendict[curmx]*curmx): curmx = curi
                        if (lendict[curi]*curi==lendict[curmx]*curmx) and curi>curmx: curmx = curi
                peak2[curp_ind] = curmx

                samevalues.append([curmx, nb3])

	peak2 = checkSmallSupport(peak2, lendict, minreads)

	if len(peak2)>2: peak2 = peak2[:2]

	if mdebug:
		for i in range(len(x)): print ('%d=%d;' % (x[i], lendict[x[i]])),
		print 'mPeak', peak2
		print '#xF=',
		for i in range(len(x)): print ('%3d' % x[i]),
		print ''
		print 'yF =',
		for i in range(len(x)): print ('%3d' % lendict[x[i]]),
		print ''

	ally = [lendict[ldkeys[0]], lendict[ldkeys[-1]]]
	if len(ldkeys)>1: last2 = -2;
	else: last2 = -1;
	for ik in range(ldkeys[1], ldkeys[last2]+1):
		if (ik in peak2) and len(ldkeys)>1: continue;
		if lendict.has_key(ik): ally.append(lendict[ik]);
		else: ally.append(0)

	#mmean = numpy.mean(numpy.array(yo[1:-1]));
	#mstd = numpy.std(numpy.array(yo[1:-1]));
	mmean = numpy.mean(numpy.array(ally));
	mstd = numpy.std(numpy.array(ally));
	if mdebug:
		print ('<%.2f,%.2f/%d>' % (mmean, mstd, len(yo[1:-1]))),
		for pf in peak2:
			print (' %d=%.6f:%.3f; ' % (pf, 1-norm.cdf((lendict[pf]-mmean)/mstd), (lendict[pf]-mmean)/mstd)), 

	#if len(peak2)<2:
	#	if lendict[ldkeys[-1]]>minreads and (len(peak2)==0 or (len(peak2)==1 and math.fabs(peak2[0]-ldkeys[-1])>1 and ((not lendict.has_key(ldkeys[-1]-1)) or (lendict.has_key(ldkeys[-1]-1) and lendict[ldkeys[-1]]>lendict[ldkeys[-1]-1]) ) )):
	#		peak2.append(ldkeys[-1])

	#print 'mtest d3s1359: ', peak2;
	
	if len(peak2)==2:  # pvalue = 0.05; ratio = 0.55
		mstr = 'In Peak:'
		for pf in peak2:
			mstr += ('for %d, pvalue=%.6f; ' % (pf, 1-norm.cdf((lendict[pf]-mmean)/mstd)))
		if not lendict.has_key(peak2[0]): mstr += (' No %d' % peak2[0]);
		elif not lendict.has_key(peak2[1]): mstr += (' No %d' % peak2[1]);
		else:  mstr += (' ratio=%.3f(%d/%d)' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]))
		#logging.info(mstr)

		if mdebug: 
			print ('%.3f=%d/%d: ' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]) ), 
		#if lendict[peak2[1]]/float(lendict[peak2[0]])<ratio or 1-norm.cdf((lendict[peak2[1]]-mmean)/mstd)>pvalue:
		#if 1-norm.cdf((lendict[peak2[1]]-mmean)/mstd)>pvalue:
		sum3_0 = lendict[peak2[0]]; sum3_1 = lendict[peak2[1]];
		if lendict.has_key(peak2[0]-1) and lendict.has_key(peak2[1]-1):
			sum3_0 += lendict[peak2[0]-1]
			sum3_1 += lendict[peak2[1]-1]
		if lendict.has_key(peak2[0]+1) and lendict.has_key(peak2[1]+1):
			sum3_0 += lendict[peak2[0]+1]
			sum3_1 += lendict[peak2[1]+1]
		
		mstr += (' >>> %.3f(%d/%d)' % (sum3_1*peak2[1]/float(sum3_0*peak2[0]), sum3_1*peak2[1], sum3_0*peak2[0]))
		if mdebug: 
			print mstr
			print sum3_1, peak2[1], '=', sum3_1*peak2[1], '/', sum3_0, peak2[0], '=', sum3_0*peak2[0], '===', sum3_1*peak2[1]/float(sum3_0*peak2[0])
		logging.info(mstr)
		
		#print 'mtest d3s1359: ', mstr
		#print 'mtest d3s1359: ', sum3_1, sum3_0, sum3_1/float(sum3_0), getLargerOverSmallRatio(peak2[0], peak2[1])

		#if sum3_1*peak2[1]/float(sum3_0*peak2[0])<0.1:
 		#
		#revised on 13 dec 2016
		#
		if sum3_1/float(sum3_0)<getLargerOverSmallRatio(peak2[0], peak2[1]):
			peak2 = peak2[:1]

	peak2.sort();
	if len(peak2)>1:
		#revised on 13 dec 2016
		#if len3dict[peak2[0]]<len3dict[peak2[1]]*0.6:
		if len3dict[peak2[0]]<len3dict[peak2[1]]*smalleratio or lendict[peak2[0]]<lendict[peak2[1]]*smalleratio:
		#if len3dict[peak2[0]]<len3dict[peak2[1]]*smalleratio and lendict[peak2[0]]<lendict[peak2[1]]*smalleratio:
                        peak2 = [peak2[1]] #, peak2[1]]


	#
	#add on Dec 12, 2016
	#
	if len(peak2)>1 and peak2[0]==peak2[1]: peak2=peak2[1:]
	
	peaksteststr = 'Peak2='+str(len(peak2))+' info:';
	for pst in peak2: peaksteststr += ' '+str(pst)
	logging.info(peaksteststr)
	
	if len(peak2)==1:
		curPeak1 = peak2[0]
		lenP1dict = {}
		for px in range(curPeak1-6, curPeak1+7):
			if lendict.has_key(px):
				lenP1dict[px] = lendict[px];
			else: pass #lenP1dict[px] = 0
		ld1keys = lenP1dict.keys(); ld1keys.sort();
		forbelowtest = False; #True;
		if forbelowtest:
			print 'test', int(minreads/2.0+0.5), 
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			print allocr; 
		lenP1dict = reviseDictAccordingV(lenP1dict, int(minreads/2.0+0.5))
		if forbelowtest:
			print 'test', int(minreads/2.0+0.5),
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			print allocr
		
		ld1keys = lenP1dict.keys(); ld1keys.sort();
		x = []; yo = []
		#x = [ld1keys[0]-1]; yo = [0]; lenP1dict[ld1keys[0]-1] = 0;
		for ldk in ld1keys:
			x.append(ldk);
			yo.append(lenP1dict[ldk]);
		#x.append(ld1keys[-1]+1); yo.append(0); lenP1dict[ld1keys[-1]+1] = 0
		peak2 = getPeaks2(x, yo, lenP1dict, 0, mdebug, int(minreads/2.0+0.5), 0)
		peak2 = checkSmallSupport(peak2, lendict, minreads)

		peaksteststr = 'Peak2 for P1='+str(len(peak2))+' info:';
		for pst in peak2: peaksteststr += ' '+str(pst)
		logging.info(peaksteststr)
	
		if curPeak1 not in peak2:
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			
			p1pstr = ' '
			for p1p in peak2: p1pstr += str(p1p)+','
			logging.info("Warning!!!! previously detected peak "+str(curPeak1)+" not in this detection:"+p1pstr[1:-1]+" for "+allocr)
			
			newpeak2 = []; 
			for curdetp in peak2:
				issame = False;
				for i in range(len(samevalues)):
					if (curPeak1 in samevalues[i][1]) and (curdetp in samevalues[i][1]) and (lendict[curPeak1]==lendict[curdetp]):
						issame = True; newpeak2.append(samevalues[i][0])
				if not issame: newpeak2.append(curdetp)
			peak2 = newpeak2
			p1pstr = ' '
			for p1p in peak2: p1pstr += str(p1p)+','
			logging.info("New peaks "+p1pstr[1:-1])
				
	if len(peak2)>2: peak2 = peak2[:2]
	if len(peak2)<2:
		if lendict[ldkeys[-1]]>minreads and (len(peak2)==0 or (len(peak2)==1 and math.fabs(peak2[0]-ldkeys[-1])>1 and ((not lendict.has_key(ldkeys[-1]-1)) or (lendict.has_key(ldkeys[-1]-1) and lendict[ldkeys[-1]]>lendict[ldkeys[-1]-1]) ) )):
			peak2.append(ldkeys[-1])
	
	peak2.sort();
	#copy from above
	if len(peak2)==2:
		print peak2, lendict #, lenP1dict
		sum3_0 = lendict[peak2[0]]; sum3_1 = lendict[peak2[1]];
		if lendict.has_key(peak2[0]-1) and lendict.has_key(peak2[1]-1):
			sum3_0 += lendict[peak2[0]-1]; sum3_1 += lendict[peak2[1]-1]
		if lendict.has_key(peak2[0]+1) and lendict.has_key(peak2[1]+1):
			sum3_0 += lendict[peak2[0]+1]; sum3_1 += lendict[peak2[1]+1]
		mstr = ('P1 <1=%d, 0=%d> >>> %.3f(%d/%d)' % (peak2[1], peak2[0], sum3_1*peak2[1]/float(sum3_0*peak2[0]), sum3_1*peak2[1], sum3_0*peak2[0]))
		logging.info(mstr)
		#if sum3_1*peak2[1]/float(sum3_0*peak2[0])<0.1:
		#
		#revised on 13 dec 2016
		#
		if sum3_1/float(sum3_0)<getLargerOverSmallRatio(peak2[0], peak2[1]):
			peak2 = peak2[:1]

	#copy from above
	peak2.sort();
	if len(peak2)>1:
		#if len3dict[peak2[0]]<len3dict[peak2[1]]*0.6 and lendict[peak2[0]]<lendict[peak2[1]]*0.6:
		if len3dict[peak2[0]]<len3dict[peak2[1]]*smalleratio or lendict[peak2[0]]<lendict[peak2[1]]*smalleratio:
		#if len3dict[peak2[0]]<len3dict[peak2[1]]*smalleratio and lendict[peak2[0]]<lendict[peak2[1]]*smalleratio:
			peak2 = [peak2[1]] #, peak2[1]]
	#
	#end the revised on 12 dec 2016
	#

	if mdebug: print 'mPeak', peak2

	peak2 = checkSmallSupport(peak2, lendict, minreads)

	return [peak2, allocr[:-1]];

#def getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, queryrep, match=10, mismatch=-9, gap_in_perf=-2, gap_in_read=-13, gap_before_after = -1):
def getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, queryrep):
	unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv, match, mismatch, gap_in_perf, gap_in_read, gap_before_after);
	if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
	newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, hmmoptions, repeatbeforeafter)
	
	return [newstr, pre0, predstats]
	

def getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id, hgfn): #='hg38.fa'):
        alignfolder = 'align/'
	if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        #unique_file_id = simulation_file_id + analysis_file_id

        fastafile = alignfolder + repeatgene + unique_file_id +'.fasta'
        #get_alg_cmd = 'samtools faidx ./hg38/Homo_sapiens_assembly38.fasta '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile

        get_alg_cmd = 'samtools faidx '+hg_reference_and_index+'/'+hgfn+' '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile
        print get_alg_cmd

        os.system(get_alg_cmd)
        if not os.path.isfile(fastafile):
                logging.error('Cannot produce '+fastafile+' for '+repeatgene)
                sys.exit(1)
        fastadata = myReadTxtFile(fastafile)
        mfadata = ''
        for li in fastadata:
            if li[0]=='>': continue;
            if li[-1]=='\n': li = li[:-1]
            mfadata = mfadata + li

	os.system('rm '+fastafile)

	return mfadata.upper()

def getHMMOptions(repeatbeforeafter, repPat, forw_rerv):
        bp = myHMM.getBasePair()
        if forw_rerv[0]=='-': repPat = myHMM.getComplementary3(bp, repPat)
        
        hmmoptions = []
        if repeatbeforeafter>0:
                hmmoptions = myHMM.getTransition_start_emission_prob(repPat)
        else:
                hmmoptions = myHMM.getTransition_start_emission_prob_without0(repPat)
        return hmmoptions

def getRepeatForGivenGene2(commonOptions, specifiedOptions, moreOptions):
#def getRepeatForGivenGene2(chr, repeatgene, gene_start_end, repeat_orig_start_end, bamfile, repPat, forw_rerv, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        chr = moreOptions['chr']
        repeatgene = moreOptions['repeatgene']
        gene_start_end = moreOptions['gene_start_end']
        repeat_orig_start_end = moreOptions['repeat_orig_start_end']
        repPat = string.strip(moreOptions['repPat'])
        forw_rerv = moreOptions['forw_rerv']

        bamfile = specifiedOptions['bamfile']
        unique_file_id = specifiedOptions['unique_file_id']
        analysis_file_id = specifiedOptions['analysis_file_id']

        isRemInDel = commonOptions['isRemInDel']
        isupdown = commonOptions['isupdown']
        isExtend = commonOptions['isExtend']
        SeqDepth = commonOptions['SeqDepth']

        #print repeatgene,
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/float(len(repPat)) #3.0
        alignfile = alignfolder + repeatgene + unique_file_id +'.alignment.txt'
        get_alg_cmd = 'samtools view '+bamfile+' ' + chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
        logging.info('Running '+get_alg_cmd)
        os.system(get_alg_cmd);
        if os.path.getsize(alignfile)==0:
           logging.info('The file %s have zero size\nTry without chr' % alignfile)
           print ('The file %s have zero size\nTry without chr' % alignfile)
           get_alg_cmd = 'samtools view '+bamfile+' ' + chr[3:]+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
           logging.info('Running '+get_alg_cmd)
           os.system(get_alg_cmd);
        logging.info('Produced ' + alignfile + ' done!');

        if not os.path.isfile(alignfile):
                logging.error('Cannot produce '+alignfile+' for '+repeatgene)
                sys.exit(1)
        aligndata = myReadTxtFile(alignfile)
        os.system('rm '+alignfile)

        mfadata = getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id, commonOptions['hgfile'])

        covermorebeforeafter = 30
        repregion_len_threhold = len(repPat); #3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]
        #repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
        #if isExtend>0 and repeat_start_end[0]<1: repeat_start_end[0]=1

        simplebeforeafter = len(repPat); #3;

        check = False; #True;
        wrongalign = 0;
        if check: repeat_beforeafter = [];

        hmmoptions = getHMMOptions(repeatbeforeafter, repPat, forw_rerv)

        repeats = [];
        for line in aligndata:
                if check:
                   beforematch = {}; aftermatch = {}; inmatch = 0; allmatch = 0; neighbeforeafter = 200;
                   beforematch[50] = 0; aftermatch[50] = 0
                   beforematch[100] = 0; aftermatch[100] = 0
                   beforematch[150] = 0; aftermatch[150] = 0
                   beforematch[200] = 0; aftermatch[200] = 0
                   beforematch[250] = 0; aftermatch[250] = 0
                   beforematch[300] = 0; aftermatch[300] = 0

                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]
                #qualifyinfo = lsp[12]
                #print 'qualifyinfo', qualifyinfo[:100]

                if pos > repeat_start_end[0] - covermorebeforeafter:
                        wrongalign += 1;
                        continue;
                        #logging.error('The start pos in ref Genome is greater than the start position of repeats' + str(pos) +' ' + str(repeat_start_end[0]));
                if not (cchr==chr or cchr==chr[3:]):  logging.error('Not same ' + cchr +' ' + chr); continue;

                numreg = re.compile('\d+')
                numinfo = numreg.findall(aligninfo)

                mdireg = re.compile('[MIDNSHPX=]{1}')
                mdiinfo = mdireg.findall(aligninfo)

                if not len(numinfo)==len(mdiinfo):
                        logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo))); continue;

                queryind = 0;
                queryrep = '';
                longer = False;
                if check:
                   totalbefore = 0; totalafter = 0;  matchinfo = '';
                   neighmatch = 0; neighref=''; neightest = ''

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                if check:
                                   if totalbefore<repeat_start_end[0]-pos:
                                      totalbefore=repeat_start_end[0]-pos;
                                   if totalafter<pos-repeat_start_end[1]:
                                      totalafter=pos-repeat_start_end[1]

                                qrepadd = False;
                                if mdi=='M':
                                        if check:
                                           faind = pos - (repeat_start_end[0] -  neighbeforeafter)
                                           if (faind>=0 and pos-repeat_start_end[0]<0) or \
                                              (pos-repeat_start_end[1]>=0 and pos-(repeat_start_end[1]+neighbeforeafter)<=0):
                                                if mfadata[faind] == aainfo[queryind]: neighmatch += 1;
                                                neighref = neighref + mfadata[faind]
                                                neightest = neightest + aainfo[queryind]

                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;

                                        if check:
                                           allmatch += 1;
                                           if pos-1 < repeat_start_end[0]:
                                                bef = repeat_start_end[0] - pos + 1
                                                bmkeys = beforematch.keys(); bmkeys.sort();
                                                for bmk in bmkeys:
                                                    if bmk>=bef: beforematch[bmk] += 1;
                                           elif pos-1 > repeat_start_end[1]:
                                                aft = pos-1 - repeat_start_end[1]
                                                afkeys = aftermatch.keys(); afkeys.sort();
                                                for afk in afkeys:
                                                    if afk >= aft: aftermatch[afk] += 1;
                                           else: inmatch += 1;

                                elif mdi =='I':
                                        qrepadd = True;
                                        queryind = queryind + 1;
                                elif mdi == 'D':
                                        pos = pos + 1;
                                elif mdi == 'S':
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi == 'H':
                                        pass;
                                elif mdi == 'P':
                                        pass;
                                else:
                                        logging.warning('Warning unknow CIGAR element ' + str(n1) + ' ' + mdi)
                                if qrepadd:
                                        #if pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter:
                                        if pos-1 >= repeat_start_end[0]-simplebeforeafter and pos-1 <= repeat_start_end[1]+simplebeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                                if check and pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter: matchinfo += mdi
                        #if pos-1 > repeat_start_end[1] + covermorebeforeafter: longer = True;
                        if pos-1 > repeat_start_end[1]+simplebeforeafter: longer = True;

                if check and len(queryrep)>=repregion_len_threhold:
                   bmkeys = beforematch.keys(); bmkeys.sort(); befstr = '';
                   for bmk in bmkeys:
                       befstr += (str(bmk)+'='+str(beforematch[bmk])+';');
                   afkeys = aftermatch.keys(); afkeys.sort(); aftstr = '';
                   for afk in afkeys:
                       aftstr += (str(afk)+'='+str(aftermatch[afk])+';');

                   print repeatgene, repeat_start_end, len(queryrep), queryrep, gene_start_end

                   newstr = '';
                   print repPat, ':', 'MI', matchinfo
                   if len(queryrep)<1000:
                        newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, queryrep)

                        #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv);
                        #if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
                        #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)

                        #newstr, pre0, predstats = myHMM.hmmpred(queryrep, repPat, forw_rerv, repeatbeforeafter)
                   print ('Gene %s(%s:%d-%d) Match for repeat(%s%s):INmatch=%d/%d(%d) test_rep=%d; beforematch:%s(%d) aftermatch:%s(%d)' % (repeatgene, chr, repeat_start_end[0], repeat_start_end[1], repPat, forw_rerv[0], inmatch, len(queryrep), neighmatch, len(newstr)/float(len(repPat)), befstr, totalbefore, aftstr, totalafter)), longer
                   #print neighref;
                   #print neightest;

                   if longer: repeat_beforeafter.append([len(newstr)/float(len(repPat))-simplebeforeafter*2/len(repPat), totalbefore, totalafter, neighmatch])

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        print len(repeats), repeat_start_end, repeat_start_end[1]-repeat_start_end[0],gene_start_end, gene_start_end[1]-gene_start_end[0] 

        if check:
           for rbf in repeat_beforeafter:
                print ('%6d %6d %6d %6d' % (rbf[0], rbf[1], rbf[2], rbf[3]))

        rptrue = []; rpfalse = []; orignial = [];
        for currep in repeats:
                #print chr, repeatgene, repPat, len(currep[1]),
                newstr = currep[1]

                pre0 = 0; predstats=''
                if len(newstr)<1000:
                        #newstr = getAlignment.correctSeq(repPat, currep[1], forw_rerv);
                        if isRemInDel>0: #isAlign>0:
                                pass #newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, currep[1])

                                #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, currep[1][repeatbeforeafter:(len(currep[1])-repeatbeforeafter)], forw_rerv);
                                #if repeatbeforeafter>0: unewstr =  currep[1][:repeatbeforeafter]+ unewstr + currep[1][(len(currep[1])-repeatbeforeafter):]
                                #
                                ##newstr, pre0, predstats = myHMM.hmmpred(currep[1], repPat, forw_rerv, repeatbeforeafter)
                                #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)


                        else:
                                pass #newstr, pre0, predstats = myHMM.hmmpred(newstr, repPat, forw_rerv, repeatbeforeafter)
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2])
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr
                #if len(currep[1])==0: continue;
                if currep[0]: #repPat, forw_rerv
                        #rptrue.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rptrue.append(len(currep[1])/3.0);
                        rptrue.append(len(currep[1])/float(len(repPat))-simplebeforeafter*2/len(repPat));
                else:
                        #rpfalse.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rpfalse.append(len(currep[1])/3.0);
                        rpfalse.append(len(currep[1])/float(len(repPat))-simplebeforeafter*2/len(repPat));

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
        logging.debug(trstr)

        #print repeatgene, repPat, rptrue
        p2, allocr = get2Peaks(rptrue, SeqDepth)

        if len(rpfalse)>0:
                flstr = 'fals ' + str(len(rpfalse)) + ' ['
                for rpf in rpfalse:
                        flstr = flstr + ('%.0f,' % rpf)
                flstr = flstr[:-1] + ']'
                logging.debug(flstr);

        logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

        '''
        for currep_ind in range(len(repeats)):
                currep = repeats[currep_ind]

                aaprinindex = -1;
                if not (currep[0]): aaprinindex = 300

                logging.debug('\t'+str(currep[0]) + ' o:' + str(len(orignial[currep_ind][0]))  +'\t'+ orignial[currep_ind][0][:aaprinindex]);
                prestr = '';
                for i in range(orignial[currep_ind][1]): prestr += ' ';
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+prestr+ orignial[currep_ind][2]);
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+orignial[currep_ind][2]);
                logging.debug('\t'+str(currep[0]) + ' p:' + str(len(currep[1])) +'\t' + prestr+ (currep[1][:aaprinindex]))
        '''
        #p2, allocr = get2Peaks(rptrue)

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions):
#def getRepeatForGivenGene(chr, repeatgene, gene_start_end, repeat_orig_start_end, bamfile, repPat, forw_rerv, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        chr = moreOptions['chr']
        repeatgene = moreOptions['repeatgene']
        gene_start_end = moreOptions['gene_start_end']
        repeat_orig_start_end = moreOptions['repeat_orig_start_end']
        repPat = moreOptions['repPat']
        forw_rerv = moreOptions['forw_rerv']

        bamfile = specifiedOptions['bamfile']
        unique_file_id = specifiedOptions['unique_file_id']
        analysis_file_id = specifiedOptions['analysis_file_id']

        isRemInDel = commonOptions['isRemInDel']
        isupdown = commonOptions['isupdown']
        isExtend = commonOptions['isExtend']
        SeqDepth = commonOptions['SeqDepth']

        #print repeatgene,
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/float(len(repPat)) #3.0
        #repeat_start_end[1] += 1;

        ##unique_file_id = simulation_file_id + analysis_file_id
        alignfile = alignfolder + repeatgene + unique_file_id +'.alignment.txt'
        get_alg_cmd = 'samtools view '+bamfile+' ' + chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
        logging.info('Running '+get_alg_cmd)
        os.system(get_alg_cmd);
        if os.path.getsize(alignfile)==0:
           logging.info('The file %s have zero size\nTry without chr' % alignfile)
           print ('The file %s have zero size\nTry without chr' % alignfile)
           get_alg_cmd = 'samtools view '+bamfile+' ' + chr[3:]+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
           logging.info('Running '+get_alg_cmd)
           os.system(get_alg_cmd);
        logging.info('Produced ' + alignfile + ' done!');

        if not os.path.isfile(alignfile):
                logging.error('Cannot produce '+alignfile+' for '+repeatgene)
                sys.exit(1)
        aligndata = myReadTxtFile(alignfile)
        os.system('rm '+alignfile)

        mfadata = getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id, commonOptions['hgfile'])

        #fastafile = alignfolder + repeatgene+'.1.fasta'
        #get_alg_cmd = 'samtools faidx ./hg38/Homo_sapiens_assembly38.fasta '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile
        #os.system(get_alg_cmd)
        #if not os.path.isfile(fastafile):
        #        logging.error('Cannot produce '+fastafile+' for '+repeatgene)
        #        sys.exit(1)
        #fastadata = myReadTxtFile(fastafile)
        #mfadata = ''
        #for li in fastadata:
        #    if li[0]=='>': continue;
        #    if li[-1]=='\n': li = li[:-1]
        #    mfadata = mfadata + li

        covermorebeforeafter = 30
        repregion_len_threhold = len(repPat) #3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]
        repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
        if isExtend>0 and repeat_start_end[0]<1: repeat_start_end[0]=1

        check = False; #True;
        wrongalign = 0;
        if check: repeat_beforeafter = [];

        hmmoptions = getHMMOptions(repeatbeforeafter, repPat, forw_rerv)

        repeats = [];
        for line in aligndata:
                if check:
                   beforematch = {}; aftermatch = {}; inmatch = 0; allmatch = 0; neighbeforeafter = 200;
                   beforematch[50] = 0; aftermatch[50] = 0
                   beforematch[100] = 0; aftermatch[100] = 0
                   beforematch[150] = 0; aftermatch[150] = 0
                   beforematch[200] = 0; aftermatch[200] = 0
                   beforematch[250] = 0; aftermatch[250] = 0
                   beforematch[300] = 0; aftermatch[300] = 0

                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]
                #qualifyinfo = lsp[12]
                #print 'qualifyinfo', qualifyinfo[:100]

                if pos > repeat_start_end[0] - covermorebeforeafter:
                        wrongalign += 1;
                        continue;
                        #logging.error('The start pos in ref Genome is greater than the start position of repeats' + str(pos) +' ' + str(repeat_start_end[0]));
                if not (cchr==chr or cchr==chr[3:]):  logging.error('Not same ' + cchr +' ' + chr); continue;

                numreg = re.compile('\d+')
                numinfo = numreg.findall(aligninfo)

                mdireg = re.compile('[MIDNSHPX=]{1}')
                mdiinfo = mdireg.findall(aligninfo)

                if not len(numinfo)==len(mdiinfo):
                        logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo))); continue;

                queryind = 0;
                queryrep = '';
                longer = False;
                if check:
                   totalbefore = 0; totalafter = 0;  matchinfo = '';
                   neighmatch = 0; neighref=''; neightest = ''

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                if check:
                                   if totalbefore<repeat_start_end[0]-pos:
                                      totalbefore=repeat_start_end[0]-pos;
                                   if totalafter<pos-repeat_start_end[1]:
                                      totalafter=pos-repeat_start_end[1]

                                qrepadd = False;
                                if mdi=='M':
                                        if check:
                                           faind = pos - (repeat_start_end[0] -  neighbeforeafter)
                                           if (faind>=0 and pos-repeat_start_end[0]<0) or \
                                              (pos-repeat_start_end[1]>=0 and pos-(repeat_start_end[1]+neighbeforeafter)<=0):
                                                if mfadata[faind] == aainfo[queryind]: neighmatch += 1;
                                                neighref = neighref + mfadata[faind]
                                                neightest = neightest + aainfo[queryind]

                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;

                                        if check:
                                           allmatch += 1;
                                           if pos-1 < repeat_start_end[0]:
                                                bef = repeat_start_end[0] - pos + 1
                                                bmkeys = beforematch.keys(); bmkeys.sort();
                                                for bmk in bmkeys:
                                                    if bmk>=bef: beforematch[bmk] += 1;
                                           elif pos-1 > repeat_start_end[1]:
                                                aft = pos-1 - repeat_start_end[1]
                                                afkeys = aftermatch.keys(); afkeys.sort();
                                                for afk in afkeys:
                                                    if afk >= aft: aftermatch[afk] += 1;
                                           else: inmatch += 1;

                                elif mdi =='I':
                                        qrepadd = True;
                                        queryind = queryind + 1;
                                elif mdi == 'D':
                                        pos = pos + 1;
                                elif mdi == 'S':
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi == 'H':
                                        pass;
                                elif mdi == 'P':
                                        pass;
                                else:
                                        logging.warning('Warning unknow CIGAR element ' + str(n1) + ' ' + mdi)
                                if qrepadd:
                                        if pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                                if check and pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter: matchinfo += mdi
                        if pos-1 > repeat_start_end[1] + covermorebeforeafter: longer = True;

                if check and len(queryrep)>=repregion_len_threhold:
                   bmkeys = beforematch.keys(); bmkeys.sort(); befstr = '';
                   for bmk in bmkeys:
                       befstr += (str(bmk)+'='+str(beforematch[bmk])+';');
                   afkeys = aftermatch.keys(); afkeys.sort(); aftstr = '';
                   for afk in afkeys:
                       aftstr += (str(afk)+'='+str(aftermatch[afk])+';');

                   print repeatgene, repeat_start_end, len(queryrep), queryrep, gene_start_end

                   newstr = '';
                   print repPat, ':', 'MI', matchinfo
                   if len(queryrep)<1000:
                        newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, queryrep)

                        #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv);
                        #if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
                        #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)

                        #newstr, pre0, predstats = myHMM.hmmpred(queryrep, repPat, forw_rerv, repeatbeforeafter)
                   print ('Gene %s(%s:%d-%d) Match for repeat(%s%s):INmatch=%d/%d(%d) test_rep=%d; beforematch:%s(%d) aftermatch:%s(%d)' % (repeatgene, chr, repeat_start_end[0], repeat_start_end[1], repPat, forw_rerv[0], inmatch, len(queryrep), neighmatch, len(newstr)/float(len(repPat)), befstr, totalbefore, aftstr, totalafter)), longer
                   #print neighref;
                   #print neightest;

                   if longer: repeat_beforeafter.append([len(newstr)/float(len(repPat)), totalbefore, totalafter, neighmatch])

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        if check:
           for rbf in repeat_beforeafter:
                print ('%6d %6d %6d %6d' % (rbf[0], rbf[1], rbf[2], rbf[3]))

        rptrue = []; rpfalse = []; orignial = [];
        for currep in repeats:
                #print chr, repeatgene, repPat, len(currep[1]),
                newstr = currep[1]

                pre0 = 0; predstats=''
                if len(newstr)<1000:
                        #newstr = getAlignment.correctSeq(repPat, currep[1], forw_rerv);
                        if isRemInDel>0: #isAlign>0:
                                newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, currep[1])

                                #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, currep[1][repeatbeforeafter:(len(currep[1])-repeatbeforeafter)], forw_rerv);
                                #if repeatbeforeafter>0: unewstr =  currep[1][:repeatbeforeafter]+ unewstr + currep[1][(len(currep[1])-repeatbeforeafter):]
                                #
                                ##newstr, pre0, predstats = myHMM.hmmpred(currep[1], repPat, forw_rerv, repeatbeforeafter)
                                #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)


                        else:
                                newstr, pre0, predstats = myHMM.hmmpred(newstr, repPat, forw_rerv, hmmoptions, repeatbeforeafter)
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2])
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr
                #if len(currep[1])==0: continue;
                if currep[0]: #repPat, forw_rerv
                        #rptrue.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rptrue.append(len(currep[1])/3.0);
                        rptrue.append(len(currep[1])/float(len(repPat))) #3.0);
                else:
                        #rpfalse.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rpfalse.append(len(currep[1])/3.0);
                        rpfalse.append(len(currep[1])/float(len(repPat))) #3.0);

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
        logging.debug(trstr)

        #print repeatgene, repPat, rptrue
        p2, allocr = get2Peaks(rptrue, SeqDepth)

        if len(rpfalse)>0:
                flstr = 'fals ' + str(len(rpfalse)) + ' ['
                for rpf in rpfalse:
                        flstr = flstr + ('%.0f,' % rpf)
                flstr = flstr[:-1] + ']'
                logging.debug(flstr);

        logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

        for currep_ind in range(len(repeats)):
                currep = repeats[currep_ind]

                aaprinindex = -1;
                if not (currep[0]): aaprinindex = 300

                logging.debug('\t'+str(currep[0]) + ' o:' + str(len(orignial[currep_ind][0]))  +'\t'+ orignial[currep_ind][0][:aaprinindex]);
                prestr = '';
                for i in range(orignial[currep_ind][1]): prestr += ' ';
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+prestr+ orignial[currep_ind][2]);
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+orignial[currep_ind][2]);
                logging.debug('\t'+str(currep[0]) + ' p:' + str(len(currep[1])) +'\t' + prestr+ (currep[1][:aaprinindex]))

        #p2, allocr = get2Peaks(rptrue)

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions={}):
#def getRepeatForKnownGene(gLoc, repeatgene, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        if moreOptions.has_key('repeatgene'): 
           repeatgene = moreOptions['repeatgene'].lower(); #repeatgene.lower()
        else: 
           repeatgene = commonOptions['repeatgene'].lower()
           moreOptions['repeatgene'] = repeatgene
        
        gLoc = commonOptions['gLoc'];
        mgloc = get_gLoc(repeatgene, gLoc);
        print mgloc

        gene_start_end = [int(mgloc[1]), int(mgloc[2])]
        repeat_start_end = [int(mgloc[3]), int(mgloc[4])]
        #      getRepeatForGivenGene(chr, repeatgene, gene_start_end, repeat_orig_start_end, bamfile, repPat, forw_rerv, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id)
        
        moreOptions['chr'] = mgloc[0]
        moreOptions['gene_start_end'] = gene_start_end
        moreOptions['repeat_orig_start_end'] = repeat_start_end
        moreOptions['repPat'] = mgloc[5]
        moreOptions['forw_rerv']= mgloc[6]
        moreOptions['mgloc'] = mgloc

        if not specifiedOptions["SepbamfileTemp"]==None: 
           specifiedOptions["bamfile"] = (specifiedOptions["SepbamfileTemp"] % moreOptions['chr'][3:])

        res = getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions) #mgloc[0], repeatgene, gene_start_end, repeat_start_end, bamfile, mgloc[5], mgloc[6],isAlign, isupdown, isExtend, unique_file_id, analysis_file_id)
        res.append(mgloc[7])
        return res;

def getRepeat(commonOptions, specifiedOptions):
#def getRepeat(gLoc, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        summary = []

        gLoc = commonOptions['gLoc']
        moreOptions = {}

        glkeys = gLoc.keys(); glkeys.sort()
        for glk in glkeys:
                moreOptions['repeatgene'] = glk
                summary.append(getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions)); #getRepeatForKnownGene(gLoc, glk, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id));

        return summary;

