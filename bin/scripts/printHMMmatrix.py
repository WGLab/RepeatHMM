
def getCompRep(comprep):
   if comprep=='0': return comprep

   repsp = comprep.split('l')
   CompRepPat = []
   for curp in repsp:
      cursall = curp.split('/')
      CompRepPat.append({})
      CompRepPat[-1]['all'] = 0
      for curR in cursall:
         CompRepPat[-1][curR[0]] = 1;
         if len(curR)>1:
            CompRepPat[-1][curR[0]] = int(curR[1:])
         CompRepPat[-1]['all'] += CompRepPat[-1][curR[0]]
      for curR in cursall:
         CompRepPat[-1][curR[0]] = CompRepPat[-1][curR[0]]/float(CompRepPat[-1]['all'])
      del CompRepPat[-1]['all']
   print 'CompRepPatc', comprep, CompRepPat
   return CompRepPat


def getCompRepFromSimple(repPat):
   CompRepPat = []
   for i in range(len(repPat)):
      CompRepPat.append({})
      CompRepPat[-1][repPat[i]] = 1

   print 'CompRepPats', repPat, CompRepPat
   return CompRepPat



def get_len_repPat(repPat, commonOptions):
	if commonOptions['CompRep']=='0':
		len_repPat = len(repPat)
	else:
		len_repPat = len(commonOptions['CompRep'])
	return len_repPat

def printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob):
                print 'states:', states
                print 'obs_symbols:', obs_symbols
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(states)): print ('%4s%2d' % (states[sj], (sj+1))),
                                print ''
                        print ('%4s%2d' % (states[si], (si+1))),
                        sum = 0;
                        for sj in range(len(states)):
                                print ('%.4f' % trainsmat[si][sj]),; sum += trainsmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(obs_symbols)): print ('%4s%2d' % (obs_symbols[sj], (sj+1))),
                                print ''
                        print ('%4s%2d' % (states[si], (si+1))),
                        sum = 0;
                        for sj in range(len(obs_symbols)):
                                print ('%.4f' % emisionmat[si][sj]),; sum += emisionmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)): print ('%6s' % states[si]),
                print ''
                sum = 0;
                for sj in range(len(states)): print ('%.4f' % startprob[sj]),; sum += startprob[sj]
                print ('\t\tsum=%.4f' % sum)



