
import re;
import os;
import sys;
import string;
import math;
import copy

import numpy;
import time
import resource
import peakutils;
import argparse;

import logging

from scipy.stats import norm

import heapq



import getAlignment
import myHMM
import myPeakDetection
import printHMMmatrix
import myRepeatReAlignment

from myheader import *

def myReadTxtFile(filename):
        f = open(filename,'r')

        data = f.readlines();
        while len(data)>0 and string.strip(data[-1])=="": data = data[:-1];
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

	locus_name='htt'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr4', '3076237', '3245687', '3076604', '3076666', 'CAG', '+19']
	locus_name='atn1'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr12', '7033626', '7051484', '7045880', '7045936', 'CAG', '+10']
	locus_name='atxn1'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr6', '16299343', '16761721', '16327867', '16327953', 'CAG', '-14']
	locus_name='atxn2'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr12', '111890018', '112037480', '112036755', '112036823', 'CAG', '-13']
	locus_name='atxn3'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr14', '92504996', '92572954', '92537355', '92537396', 'CAG', '-14']
	locus_name='cacna1a'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr19', '13317257', '13617274', '13318673', '13318711', 'CAG', '-8']
	locus_name='atxn7'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr3', '63850233', '63989136', '63898362', '63898391', 'CAG', '+7']
	locus_name='tbp'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr6', '170863421', '170881958', '170870987', '170871109', 'CAG', '+19']
	locus_name='fxn'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr9', '71650479', '71715094', '71652203', '71652221', 'GAA', '+5']
	locus_name='dmpk'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr19', '46272967', '46285815', '46273463', '46273522', 'CTG', '-20']
	locus_name='atxn8os'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr13', '70681345', '70713885', '70713483', '70713560', 'CTG', '+15']
	locus_name='ppp2r2b'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chr5', '145969068', '146461083', '146258292', '146258322', 'CAG', '-10']
	locus_name='ar'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chrX', '66763874', '66950461', '66765160', '66765228', 'CAG', '+22']
	locus_name='fmr1'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chrX', '146993469', '147032647', '146993569', '146993628', 'CGG', '+10']
	locus_name='aff2'; locus_name = locus_name.lower();
	gLoc[locus_name] = ['chrX', '147582039', '148082293', '147582126', '147582212', 'CCG', '+15']
	#gLoc[locus_name] = ['chrX', '147582139', '148082193', '147582126', '147582212', 'CCG', '+15']
	#   gLoc['aff2'] = ['chrX', '148500619', '149000663', '148500606', '148500692', 'CCG', '+15', '6-35:200+']


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
	gLoc[locus_name] = ['chr2', '68139016','68339157', '68239079','68239126', 'TCTA', plusminus+'12', '7-17:7-17']

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
	gLoc[locus_name] = ['chr2','68239015','68239157','68239078','68239127','TCTA','+12','7-17:7-17',]
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
	gLoc[locus_name] = ['chr8','42536554','42536839','42536588','42536615','AAT','+9','9-22:9-22',]
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
	gLoc[locus_name] = ['chr8','16770880','16771227','16771051','16771172','ATAG','+30','20-34:20-34',]
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
			print glk, ('Wrong gene and repeat information: gene = [%d, %d], and repeat = [%d, %d]' % (genstr,genend,repstr,repend))
			#print glk, ('Wrong gene and repeat information: gene = [%d, %d], and repeat = [%d, %d]' % genstr,genend,repstr,repend)

	return gLoc;

def getDiseseGeneInRefGenomeLocation(hg):
	if hg=='' or hg=='hg38': return getDiseseGeneInRefGenomeLocation_hg38();
	elif hg=='hg19': return getDiseseGeneInRefGenomeLocation_hg19();
	else:
		logging.error('Error: not supported yet')

#from wiki: https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
#RefSeq chr:start-end would be used.
def getDiseseGeneInRefGenomeLocation_hg38():
        gLoc = {} # htt, atn1, atxn1, atxn2, atxn3, cacna1a, atxn7, tbp, fxn, dmpk, atxn8os, ppp2r2b, ar, fmr1, aff2
        #HTT at chr4:3074510-3243960               3074400   3076157
        gLoc['htt'] = ['chr4','3074510','3243960','3074877','3074939', 'CAG', '+19', '6-35:36-250']                # (CAG)* +1CAACAG 
        #ATN1 at         chr12:   6924463 -  6942321   6936729   6936773
        gLoc['atn1'] = ['chr12', '6924463', '6942321','6936717', '6936773', 'CAG', '+10', '6-35:49-88']   # 2(CAGCAA) (CAG)*
        #ar: (Q)* <+ 5non-Q + 7 Q>           6936346  /6937032
        gLoc['ar'] = ['chrX', '67544032', '67730619', '67545318', '67545386', 'CAG', '+22', '9-36:38-62']         # (CAG)* +1CAA
        gLoc['atxn1'] = ['chr6', '16299112', '16761490', '16327636', '16327722', 'CAG', '-14', '6-35:49-88']      # (CAG)* +2(ATGCTG) 11(CTG)
        gLoc['atxn2'] = ['chr12', '111452214', '111599676', '111598951', '111599019', 'CAG', '-13', '14-32:33-77'] # 10Q  (CAG)*
        #atxn3: chr:start-end for RefSeq is chr14:92058552-92106621. chr:start-end for UCSC is chr14:92038652-92106610
        #ATXN3 at chr14:92058552-92106621 
        gLoc['atxn3'] = ['chr14', '92038652', '92106610', '92071011', '92071052', 'CAG', '-14', '12-40:55-86']       # (CAG)* + CTGTTGCTGCTTTTGCTGCTG
        gLoc['atxn3'] = ['chr14', '92038652', '92106610', '92071011', '92071052', 'CAG', '-14', '12-40:55-150']
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
        gLoc['ppp2r2b'] = ['chr5', '146589505', '147081520', '146878729', '146878759', 'CAG', '-10', '7-28:66-78'] #may contain errors;  aff2, fxn, dmpk, atxn8os, ppp2r2b

        locus_name='ATXN10'; locus_name = locus_name.lower();    #$ atxn10
        #gLoc[locus_name] = ['chr22', '45671798', '45845307', '45795355', '45795424', 'ATTCT','+14','10-32',]
        gLoc[locus_name] = ['chr22', '45792355', '45798424', '45795355', '45795424', 'ATTCT','+14','10-32']

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
        gLoc[locus_name] = ['chr2','68011883','68012025','68011946','68011995','TCTA','+12','7-17:7-17',]
        locus_name='D10S1248'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','129294109','129294532','129294243','129294295','GGAA','+13','',]
        locus_name='D22S1045'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr22','37140244','37140530','37140286','37140342','ATT','+18','11-25:11-25',]
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
        gLoc[locus_name] = ['chr8','42681411','42681696','42681445','42681472','AAT','+9','9-22:9-22',]
        locus_name='D10S2325'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr10','12750946','12751246','12751051','12751127','ATAAG','+15','8-22:8-22',]
        locus_name='D11S554'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr11','44908682','44908921','44908734','44908830','AAAG','+23','',]
        locus_name='D14S306'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr14','37858877','37859246','37859085','37859187','TCTA','+25','',]
        locus_name='D18S535'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','40568810','40569365','40568861','40568930','TAGA','+17','',]
        locus_name='D18S849'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr18','57897345','57897656','57897387','57897526','TCTA','+35','28-42:28-42',]
        locus_name='D19S253'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr19','15617324','15617553','15617483','15617534','ATCT','+12','',]
        locus_name='D1S103'; locus_name = locus_name.lower();    #$
        gLoc[locus_name] = ['chr1','230700989','230701240','230701148','230701187','TG','+19','12-26:12-26',]
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
        gLoc[locus_name] = ['chr8','16913371','16913718','16913542','16913663','ATAG','+30','20-34:20-34',]
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
        gLoc[locus_name] = ['chr2','68011883','68012025','68011946','68011995','TCTA','+12','7-17:7-17',]
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
        gLoc[locus_name] = ['chr22','37140244','37140530','37140286','37140342','ATT','+18','11-25:11-25',]
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

        locus_name='MUC1'; locus_name = locus_name.lower();
        gLoc[locus_name] = ['chr1', '155186123', '155192868', '155188508', '155192051', 'GGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGCTGGGGGGGC', '+40', '']

        return gLoc;

def get_gLoc(dis_type, gLoc):
        dis_type = dis_type.lower();
        return gLoc[dis_type];



def get2Peaks(lengd, MinSup, commonoptions=None):
	return myPeakDetection.get2Peaks(lengd, MinSup, True, commonoptions) 
		
	
def getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, queryrep, commonOptions):
	unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv, match, mismatch, gap_in_perf, gap_in_read, gap_before_after);
	if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]

	unewstr = unewstr.replace('-', '')
	unewstr = unewstr.replace('N', ''); unewstr = unewstr.replace('n', '');

	if len(unewstr)>1:
		newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, hmmoptions, commonOptions)
	else:
		newstr, pre0, predstats = '', '', ''
	
	return [newstr, pre0, predstats]
	

def getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id, hgfn): #='hg38.fa'
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        fastafile = alignfolder + repeatgene + unique_file_id +'.fasta'

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

def getHMMOptions(repeatbeforeafter, repPat, forw_rerv, commonOptions):
        hmmoptions = []
        hmmoptions = myHMM.getTransition_start_emission_prob(repPat, commonOptions)
        return hmmoptions

def getRepeatForGivenGene2(commonOptions, specifiedOptions, moreOptions):
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
        MinSup = commonOptions['MinSup']

        len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
        logging.info("len_repPat="+str(len_repPat))

        #print repeatgene,
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/float(len_repPat) #3.0
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
        repregion_len_threhold = len_repPat; #3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]

        simplebeforeafter = len_repPat; #3;

        wrongalign = 0;
        hmmoptions = getHMMOptions(repeatbeforeafter, repPat, forw_rerv, commonOptions)

        repeats = [];
        for line in aligndata:
                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]

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

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                qrepadd = False;
                                if mdi=='M':
                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;
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
                                        if pos-1 >= repeat_start_end[0]-simplebeforeafter and pos-1 <= repeat_start_end[1]+simplebeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                        if pos-1 > repeat_start_end[1]+simplebeforeafter: longer = True;

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        print len(repeats), repeat_start_end, repeat_start_end[1]-repeat_start_end[0],gene_start_end, gene_start_end[1]-gene_start_end[0] 

        rptrue = []; rpfalse = []; orignial = [];
        for currep in repeats:
                newstr = currep[1]

                pre0 = 0; predstats=''
                if len(newstr)<commonOptions['MaxRep']*len_repPat:
                        if isRemInDel>0: 
                                pass 
                        else:
                                pass 
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2]+" "+str(commonOptions['MaxRep'])+" "+str(commonOptions['MaxRep']*len_repPat))
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr
                if currep[0]: 
                        rptrue.append(len(currep[1])/float(len_repPat)-simplebeforeafter*2/len_repPat);
                else:
                        rpfalse.append(len(currep[1])/float(len_repPat)-simplebeforeafter*2/len_repPat);

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
        logging.debug(trstr)

        p2, allocr = get2Peaks(rptrue, MinSup, commonoptions=commonOptions)

        if len(rpfalse)>0:
                flstr = 'fals ' + str(len(rpfalse)) + ' ['
                for rpf in rpfalse:
                        flstr = flstr + ('%.0f,' % rpf)
                flstr = flstr[:-1] + ']'
                logging.debug(flstr);

        logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions):
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
        MinSup = commonOptions['MinSup']

        len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
        logging.info("len_repPat="+str(len_repPat))
        #print repeatgene,
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/float(len_repPat) #3.0

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
        repregion_len_threhold = len_repPat #3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]
        repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
        if isExtend>0 and repeat_start_end[0]<1: repeat_start_end[0]=1

        wrongalign = 0;

        hmmoptions = getHMMOptions(repeatbeforeafter, repPat, forw_rerv, commonOptions)

        repeats = [];
        for line in aligndata:
                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]

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

                queryind = 0; hpadd = 0;
                queryrep = '';
                longer = False;

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                qrepadd = False;
                                if mdi=='M':
                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi =='I':
                                        qrepadd = True;
                                        queryind = queryind + 1;
                                elif mdi == 'D':
                                        pos = pos + 1;
                                elif mdi == 'S':
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi == 'H':
                                        if qrepadd: hpadd += 1; #pass;
                                elif mdi == 'P':
                                        if qrepadd: hpadd += 1; #pass;
                                else:
                                        logging.warning('Warning unknow CIGAR element ' + str(n1) + ' ' + mdi)
                                if qrepadd:
                                        if pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                        if pos-1 > repeat_start_end[1] + covermorebeforeafter: longer = True;

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        rptrue = []; rpfalse = []; orignial = [];
        for currep in repeats:
                newstr = currep[1]

                pre0 = 0; predstats=''
                if len(newstr)<commonOptions['MaxRep']*len_repPat:
                      newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, currep[1], commonOptions)
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2]+" "+str(commonOptions['MaxRep'])+" "+str(commonOptions['MaxRep']*len_repPat))
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr
                if currep[0]: 
                        rptrue.append(len(currep[1])/float(len_repPat)) #3.0);
                else:
                        rpfalse.append(len(currep[1])/float(len_repPat)) #3.0);

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
        logging.debug(trstr)

        p2, allocr = get2Peaks(rptrue, MinSup, commonoptions=commonOptions)

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
                logging.debug('\t'+str(currep[0]) + ' p:' + str(len(currep[1])) +'\t' + prestr+ (currep[1][:aaprinindex]))

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def fixsize2(p2, ind):
	if len(p2[ind])>1: pass
	elif len(p2[ind])==1:
		p2[ind].append(p2[ind][0])
	else: 
		p2[ind].append(0); 
		p2[ind].append(0);

def addSumForAGene(p2, myret, myretdetail, mstr, ind=2):
	fixsize2(p2, ind)
	myretdetail[mstr] = (('%10s' % mstr)+' '+str(p2));
	print mstr, p2
	logging.info(mstr+' '+str(p2)+'\n')
	sys.stdout.flush()
	myret[mstr] = p2[ind]
 
def getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions={}):
        if moreOptions.has_key('repeatgene'): 
           repeatgene = moreOptions['repeatgene'].lower(); #repeatgene.lower()
        else: 
           repeatgene = commonOptions['repeatgene'].lower()
           moreOptions['repeatgene'] = repeatgene
        
        gLoc = commonOptions['gLoc'];
        repeatgene = commonOptions['repeatgene'].lower()
        repeatgene = moreOptions['repeatgene'].lower()
        newinfo = commonOptions['specifiedGeneInfo']

        infospt = newinfo.split('/')
        repeatgene = repeatgene.lower()
        if gLoc.has_key(repeatgene):
           logging.info(repeatgene)
           mgloc = gLoc[repeatgene]
        else: mgloc = ['','','', '','', '','','']

        if len(infospt)<len(mgloc):
           logging.error("Error: wrong input for the gene of interes: %s whose length is %d less than the expected length %d" % (newinfo, len(infospt), len(mgloc)));
           sys.exit(101);

        for i in range(len(mgloc)):
           curnew = string.strip(infospt[i]);
           if not curnew=='':
              mgloc[i] = curnew
        errorstr = ''
        for i in range(len(mgloc)):
           if string.strip(mgloc[i])=='' and i<len(mgloc)-1:
              errorstr += ('Error no information for %d\n' % i)
        if not errorstr=='':
           logging.error(errorstr);
           print errorstr, mgloc, repeatgene
           sys.exit(102);

        print 'mgloc', mgloc

        gene_start_end = [int(mgloc[1]), int(mgloc[2])]
        repeat_start_end = [int(mgloc[3]), int(mgloc[4])]
        
        moreOptions['chr'] = mgloc[0]
        moreOptions['gene_start_end'] = gene_start_end
        moreOptions['repeat_orig_start_end'] = repeat_start_end
        moreOptions['repPat'] = mgloc[5]
        moreOptions['forw_rerv']= mgloc[6]
        moreOptions['mgloc'] = mgloc
        
        myHMM.produce_for_repPat(commonOptions, moreOptions)

        if not specifiedOptions["SepbamfileTemp"]==None: 
           specifiedOptions["bamfile"] = (specifiedOptions["SepbamfileTemp"] % moreOptions['chr'][3:])

        myret = {}; myretdetail = {}
        if testall:
           start_time = time.time();
           print 'p2bamself start'
           p2bamself = getRepeatForGivenGene2(commonOptions, specifiedOptions, moreOptions)
           memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
           addSumForAGene(p2bamself, myret, myretdetail, 'p2bamself', 2)
           end_time = time.time();
           print ('p2bamself end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()
        if (commonOptions['SplitAndReAlign'] in [0,2]) or testall:
           start_time = time.time();
           print 'p2bamhmm start'
           p2bamhmm = getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions) 
           memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
           addSumForAGene(p2bamhmm, myret, myretdetail, 'p2bamhmm', 2)
           end_time = time.time();
           print ('p2bamhmm end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()
        if (commonOptions['SplitAndReAlign'] in [1,2]) or testall:
           start_time = time.time();
           print 'p2sp start'
           moreOptions['fafqfile'] = specifiedOptions["bamfile"]
           moreOptions['fafqtype'] = 'bam'
           p2sp = myRepeatReAlignment.getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
           memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
           addSumForAGene(p2sp, myret, myretdetail, 'p2sp', 2)
           end_time = time.time();
           print ('p2sp end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

        return [myret, myretdetail];

def getRepeat(commonOptions, specifiedOptions):
        summary = {};

        gLoc = commonOptions['gLoc']
        moreOptions = {}

        glkeys = gLoc.keys(); glkeys.sort()
        for glk in glkeys:
                moreOptions['repeatgene'] = glk
                #print 'commonOptions', commonOptions
                #print 'specifiedOptions', specifiedOptions
                #print 'moreOptions', moreOptions
                summary[glk] = getRepeatForKnownGene(commonOptions, specifiedOptions, moreOptions)

        return summary;

