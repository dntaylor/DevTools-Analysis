#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import itertools
import operator

import ROOT

logger = logging.getLogger("ThreePhotonAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class ThreePhotonAnalysis(AnalysisBase):
    '''
    WZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','threePhotonTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ThreePhotonTree')
        super(ThreePhotonAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')

        # trigger
        self.addPhotonTriggers()
        #self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        #self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        #self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # 3 photon
        self.addComposite('ggg')

        # z leptons
        self.addDiCandidate('gg12')
        self.tree.add(lambda cands: ROOT.TMath.ACos(2*cands['gg12'].M()/cands['ggg'].M()-1)/ROOT.TMath.Pi(), 'gg12_mPrime', 'F')
        self.tree.add(lambda cands: cands['gg12Boost'].Theta()/ROOT.TMath.Pi(), 'gg12_thetaPrime', 'F')
        self.addDiCandidate('gg13')
        self.tree.add(lambda cands: ROOT.TMath.ACos(2*cands['gg13'].M()/cands['ggg'].M()-1)/ROOT.TMath.Pi(), 'gg13_mPrime', 'F')
        self.tree.add(lambda cands: cands['gg13Boost'].Theta()/ROOT.TMath.Pi(), 'gg13_thetaPrime', 'F')
        self.addDiCandidate('gg23')
        self.tree.add(lambda cands: ROOT.TMath.ACos(2*cands['gg23'].M()/cands['ggg'].M()-1)/ROOT.TMath.Pi(), 'gg23_mPrime', 'F')
        self.tree.add(lambda cands: cands['gg23Boost'].Theta()/ROOT.TMath.Pi(), 'gg23_thetaPrime', 'F')
        self.addPhoton('g1',doId=True,doScales=True)
        self.addPhoton('g2',doId=True,doScales=True)
        self.addPhoton('g3',doId=True,doScales=True)

        # met
        self.addMet('met')

    ############################
    ### select WZ candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'g1' : None,
            'g2' : None,
            'g3' : None,
            'gg12': None,
            'gg13': None,
            'gg23': None,
            'ggg': None,
            'met': self.pfmet,
            'cleanJets': [],
            #'leadJet': Candidate(None),
            #'subleadJet': Candidate(None),
            #'dijet': DiCandidate(Candidate(None),Candidate(None)),
        }

        gs = self.getPassingCands('ElectronVeto',self.photons)

        best = ()

        for gCand in itertools.permutations(gs,3):
            g1, g2, g3 = gCand
            if not (g1.pt()>g2.pt() and g2.pt()>g3.pt()): continue
            if g1.pt()<32: continue
            if g2.pt()<20: continue
            if not self.passPhotonPreselection(g1): continue
            if not self.passPhotonPreselection(g2): continue
            if not best: best = gCand
            if g1.pt()>best[0].pt():
                best = gCand
            elif g1.pt()==best[0].pt() and g2.pt()>best[1].pt():
                best = gCand
            elif g1.pt()==best[0].pt() and g2.pt()==best[1].pt() and g3.pt()>best[2].pt():
                best = gCand
            

        if not best: return candidate

        g1, g2, g3 = best
        candidate['g1'] = g1
        candidate['g2'] = g2
        candidate['g3'] = g3
        candidate['gg12'] = DiCandidate(g1,g2)
        candidate['gg13'] = DiCandidate(g1,g3)
        candidate['gg23'] = DiCandidate(g2,g3)
        candidate['ggg'] = CompositeCandidate(g1,g2,g3)

        goodGs = self.getPassingCands('Photon',self.photons)
        candidate['cleanJets'] = self.cleanCands(self.jets,goodGs+[g1,g2,g3],0.4)

        # boosted candidates
        candidate['gg12Boost'] = candidate['gg12'].p4()
        candidate['gg13Boost'] = candidate['gg13'].p4()
        candidate['gg23Boost'] = candidate['gg23'].p4()
        candidate['gg12Boost'].Boost(-candidate['gg12'].BoostVector())
        candidate['gg13Boost'].Boost(-candidate['gg13'].BoostVector())
        candidate['gg23Boost'].Boost(-candidate['gg23'].BoostVector())

        #g1 = candidate['g1']
        #g2 = candidate['g2']
        #gg12 = candidate['gg12']

        #g1p4 = g1.p4()
        #g2p4 = g2.p4()
        #gg12p4 = gg12.p4()

        #print 'original'

        #print 'g1', g1p4.Px(), g1p4.Py(), g1p4.Pz(), g1p4.E(), g1p4.Theta()
        #print 'g2', g2p4.Px(), g2p4.Py(), g2p4.Pz(), g2p4.E(), g2p4.Theta()
        #print 'gg12', gg12p4.Px(), gg12p4.Py(), gg12p4.Pz(), gg12p4.E(), gg12p4.Theta()

        #print 'boosting'

        #g1p4.Boost(-gg12.BoostVector())
        #g2p4.Boost(-gg12.BoostVector())
        #gg12p4.Boost(-gg12.BoostVector())

        #print 'g1', g1p4.Px(), g1p4.Py(), g1p4.Pz(), g1p4.E(), g1p4.Theta()
        #print 'g2', g2p4.Px(), g2p4.Py(), g2p4.Pz(), g2p4.E(), g2p4.Theta()
        #print 'gg12', gg12p4.Px(), gg12p4.Py(), gg12p4.Pz(), gg12p4.E(), gg12p4.Theta()

        return candidate


    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        # accept MC, check trigger for data
        #if self.event.isData()<0.5: return True
        # use trigger for MC and data
        isData = self.event.isData()>0.5
        triggerNames = {
            'DoubleEG'       : [
                'DoublePhoton60',
                'Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90',
            ],
            'SinglePhoton' : [
                'Photon175',
            ],
        }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleEG', 
            'SinglePhoton',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    #def triggerEfficiencyMC(self,cands):
    #    return self.triggerEfficiency(cands,mode='mc')

    #def triggerEfficiencyData(self,cands):
    #    return self.triggerEfficiency(cands,mode='data')

    #def triggerEfficiency(self,cands,mode='ratio'):
    #    candList = [cands[c] for c in ['g1','g2','g3']]
    #    triggerList = []
    #    if mode=='data':
    #        return self.triggerScales.getDataEfficiency(triggerList,candList)
    #    elif mode=='mc':
    #        return self.triggerScales.getMCEfficiency(triggerList,candList)
    #    elif mode=='ratio':
    #        return self.triggerScales.getRatio(triggerList,candList)








def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('DoubleEG',version='80XPhoton'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='threePhotonTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    threePhotonAnalysis = ThreePhotonAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ThreePhotonTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       threePhotonAnalysis.analyze()
       threePhotonAnalysis.finish()
    except KeyboardInterrupt:
       threePhotonAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
