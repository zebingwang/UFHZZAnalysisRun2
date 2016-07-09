from ROOT import *

infile = TFile('/scratch/osg/dsperka/Sync_80X.root','READ')
#ext=''
#ext='_norefit'
ext='_short'

#infile = TFile('Sync_80X_HighMass.root','READ')
#ext='_highmass'

#infile = TFile('test.root','READ')
#ext='_test'

t = infile.Get('Ana/passedEvents')
nevents = t.GetEntriesFast()

n4e=0
n4mu=0
n2e2mu=0

with open('eventsBUF_80X'+ext+'.txt','w') as f:

  for i in xrange(nevents):
    t.GetEntry(i)

    if (not (t.passedZ4lSelection==1 and t.massZ2>12.0)): continue

    if (t.njets_pt30_eta4p7>0): 
      ptjet1=t.jet_pt[t.jet_iscleanH4l[0]]
      qg1=t.jet_QGTagger[t.jet_iscleanH4l[0]]
    else: 
      ptjet1=-1.0
      qg1=-1.0

    if (t.njets_pt30_eta4p7>1): 
      ptjet2=t.jet_pt[t.jet_iscleanH4l[1]]
      qg2=t.jet_QGTagger[t.jet_iscleanH4l[1]]
    else: 
      ptjet2=-1.0
      qg2=-1.0

    if (ext=='_short'):
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.2f}:{:.2f}:{:.2f}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.mass4lErr,t.massZ1,t.massZ2,t.mass4lREFIT,t.mass4lErrREFIT)+'\n')
    elif (ext=='_norefit'):
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:d}:{:.2f}:{:.2f}:{}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.massZ1,t.massZ2,t.D_bkg_kin,t.D_bkg,t.Dgg10_VAMCFM,t.Djet_VAJHU,t.D_g4,t.njets_pt30_eta4p7,ptjet1,ptjet2,t.EventCat)+'\n')
    else:
      #f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:d}:{:.2f}:{:.2f}:{}:{:.2f}:{:.2f}:{:.3f}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.massZ1,t.massZ2,t.D_bkg_kin,t.D_bkg,t.Dgg10_VAMCFM,t.Djet_VAJHU,t.D_g4,t.D_VBF1j,t.D_WHh,t.D_ZHh,t.njets_pt30_eta4p7,ptjet1,ptjet2,t.EventCat,t.mass4lREFIT,t.mass4lErrREFIT,t.genWeight)+'\n')
      f.write('{}:{}:{}:{:.2f}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:d}:{:.2f}:{:.2f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.2f}:{:.2f}:{:.3f}'.format(t.Run,t.LumiSect,t.Event,t.mass4l,t.massZ1,t.massZ2,t.D_bkg_kin,t.D_bkg,t.Dgg10_VAMCFM,t.Djet_VAJHU,t.D_g4,t.D_VBF1j_VAJHU,t.D_WHh_VAJHU,t.D_ZHh_VAJHU,t.njets_pt30_eta4p7,ptjet1,ptjet2,qg1,qg2,t.D_VBF2j,t.D_VBF1j,t.D_WHh,t.D_ZHh,t.EventCat,t.mass4lREFIT,t.mass4lErrREFIT,t.genWeight)+'\n')

    if (t.finalState==1): n4mu+=1  
    if (t.finalState==2): n4e+=1  
    if (t.finalState>=3): n2e2mu+=1
  
  f.close()

print nevents,'events in the file'
print '4e:',n4e,'4mu:',n4mu,'2e2mu:',n2e2mu
