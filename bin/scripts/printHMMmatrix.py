

def printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob):
                print 'states:', states
                print 'obs_symbols:', obs_symbols
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(states)): print ('%6s' % states[sj]),
                                print ''
                        print ('%6s' % states[si]),
                        sum = 0;
                        for sj in range(len(states)):
                                print ('%.4f' % trainsmat[si][sj]),; sum += trainsmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(obs_symbols)): print ('%6s' % obs_symbols[sj]),
                                print ''
                        print ('%6s' % states[si]),
                        sum = 0;
                        for sj in range(len(obs_symbols)):
                                print ('%.4f' % emisionmat[si][sj]),; sum += emisionmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)): print ('%6s' % states[si]),
                print ''
                sum = 0;
                for sj in range(len(states)): print ('%.4f' % startprob[sj]),; sum += startprob[sj]
                print ('\t\tsum=%.4f' % sum)



