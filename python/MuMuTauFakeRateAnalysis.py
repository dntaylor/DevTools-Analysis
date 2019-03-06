#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MuMuTauFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class MuMuTauFakeRateAnalysis(AnalysisBase):
    '''
    MuMuTauFakeRate analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','muMuTauFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MuMuTauFakeRateTree')
        super(MuMuTauFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[0], 'triggerEfficiencyMC', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[1], 'triggerEfficiencyMCUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[2], 'triggerEfficiencyMCDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[0], 'triggerEfficiencyData', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[1], 'triggerEfficiencyDataUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[2], 'triggerEfficiencyDataDown', 'F')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('z1')
        self.addDetailedMuon('z1')
        self.addLepton('z2')
        self.addDetailedMuon('z2')

        # tau
        self.addDiLepton('mt')
        self.addLepton('t')
        self.addDetailedTau('t')
        self.tree.add(lambda cands: self.tauGenMatch(cands,'t'), 't_tauGenMatch2', 'F')
        self.addJet('tjet')
        self.addLepton('m')
        self.addDetailedMuon('m')

        # met
        self.addMet('met')


    ############################
    ### Additional functions ###
    ############################
    def passMuon(self,cand):
        if cand.pt()<3: return False
        if abs(cand.dxy())>=0.2: return False
        if abs(cand.dz())>=0.5: return False
        if not cand.isPFMuon(): return False
        if not (cand.isGlobalMuon() or cand.isTrackerMuon()): return False
        return True

    def passTau(self,cand):
        if cand.pt()<10: return False
        if abs(cand.dxy())>=0.2: return False
        if abs(cand.dz())>=0.5: return False
        if not cand.decayModeFinding(): return False
        return True

    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            't' : None,
            'tjet' : Candidate(None),
            'z' : None,
            'm': None,
            'mt': None,
            'met': self.pfmet,
        }

        # get leptons
        muons = [m for m in self.muons if self.passMuon(m)]
        taus = [t for t in self.taus if self.passTau(t)]
        leps = muons+taus
        if len(muons)<2: return candidate
        if len(taus)<1: return candidate


        # get the candidates
        zCand = []
        for mm in itertools.permutations(muons,2):
            m1 = mm[0]
            m2 = mm[1]
            # charge OS
            if m1.charge()==m2.charge(): continue
            # require lead m pt>25
            if m1.pt()<26: continue
            if m2.pt()<10: continue
            # iso
            if m1.relPFIsoDeltaBetaR04()>=0.25: continue
            if m2.relPFIsoDeltaBetaR04()>=0.25: continue
            # make composites
            z = DiCandidate(m1,m2)
            if not zCand: zCand = mm
            better = True
            if zCand[0].pt()>m1.pt():
                better = False
            elif zCand[1].pt()>m2.pt():
                better = False
            if better: zCand = mm
                
        if not zCand: return candidate

        m1 = zCand[0] if zCand[0].pt()>zCand[1].pt() else zCand[1]
        m2 = zCand[1] if zCand[0].pt()>zCand[1].pt() else zCand[0]

        z = DiCandidate(m1,m2)
        if z.M()>120: return candidate
        if z.M()<60: return candidate

        # highest pt tau with DR>0.8 from selected muons
        goodTaus = [t for t in taus if deltaR(t.eta(),t.phi(),m1.eta(),m1.phi())>0.8 and deltaR(t.eta(),t.phi(),m2.eta(),m2.phi())>0.8]
        if len(goodTaus)==0: return candidate

        t = goodTaus[0]

        # find if there is a muon nearby
        otherMuons = [m for m in muons if m!=m1 and m!=m2 and deltaR(m.eta(),m.phi(),t.eta(),t.phi())<0.8]
        otherMuonsByDR = sorted(otherMuons, key=lambda m: deltaR(m.eta(),m.phi(),t.eta(),t.phi()))
        nearM = otherMuonsByDR[0] if otherMuonsByDR else Candidate(None)

        candidate['z1'] = m1
        candidate['z2'] = m2
        candidate['t'] = t
        candidate['z'] = DiCandidate(m1,m2)
        candidate['m'] = nearM
        candidate['mt'] = DiCandidate(nearM,t) if len(otherMuonsByDR) else DiCandidate(Candidate(None),Candidate(None))

        result = self.tauGenMatch(candidate,'t')

        # match jet to tau
        dr = 999
        j = None
        for jet in self.jets:
            jt = DiCandidate(jet,t)
            if jt.deltaR()<dr:
                j = jet
                dr = jt.deltaR()
        if j:
            candidate['tjet'] = j

        return candidate

    def tauGenMatch(self,candidate,label):
        if self.event.isData(): return -1
        t = candidate[label]
        gtaus = [g for g in self.gen if abs(g.pdgId())==15]
        gelectrons = [g for g in self.gen if abs(g.pdgId())==11]
        gmuons = [g for g in self.gen if abs(g.pdgId())==13]
        bestG = None
        bestDR = 999.
        for g in self.gen:
            if not g.isPrompt(): continue
            if g.pt()<8: continue
            dr = deltaR(t.eta(),t.phi(),g.eta(),g.phi())
            if dr<bestDR:
                bestG = g
                bestDR = dr
        bestTauG = None
        bestTauDR = 999
        for g in gtaus:
            dr = deltaR(t.eta(),t.phi(),g.eta(),g.phi())
            if dr<bestTauDR:
                bestTauG = g
                bestTauDR = dr
        #if bestDR<0.2:
        #    dr = deltaR(t.eta(),t.phi(),bestG.eta(),bestG.phi())
        #    print bestG.pdgId(), bestG.numberOfDaughters(), bestG.daughter_1(), bestG.daughter_2(), bestG.mother_1(), dr
        #else:
        #    print 'No prompt'
        #if bestTauDR<0.2:
        #    dr = deltaR(t.eta(),t.phi(),bestTauG.eta(),bestTauG.phi())
        #    print bestTauG.pdgId(), bestTauG.numberOfDaughters(), bestTauG.daughter_1(), bestTauG.daughter_2(), bestTauG.mother_1(), dr
        #else:
        #    print 'No gen tau'
        if bestDR<0.2 and bestDR<bestTauDR:
            if abs(bestG.pdgId()) == 11: return 1
            elif abs(bestG.pdgId()) == 13: return 2
        elif bestTauDR<0.2 and bestTauDR<bestDR:
            if abs(bestTauG.pdgId()) == 15  and 11 in [abs(bestTauG.daughter_1()), abs(bestTauG.daughter_2())]: return 3
            elif abs(bestTauG.pdgId()) == 15  and 13 in [abs(bestTauG.daughter_1()), abs(bestTauG.daughter_2())]: return 4
            elif abs(bestTauG.pdgId()) == 15: return 5

        return 6



    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        #isData = self.event.isData()>0.5
        #if self.version=='76X':
        #    triggerNames = {
        #        'DoubleMuon'     : [
        #            'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
        #            'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
        #        ],
        #    }
        #else:
        #    triggerNames = {
        #        'DoubleMuon'     : [
        #            'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
        #            'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
        #        ],
        #    }
        ## the order here defines the heirarchy
        ## first dataset, any trigger passes
        ## second dataset, if a trigger in the first dataset is found, reject event
        ## so forth
        #datasets = [
        #    'DoubleMuon',
        #]
        #return self.checkTrigger(*datasets,**triggerNames)

        isData = self.event.isData()>0.5
        if self.version=='76X':
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
            }
        else:
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu24',
                    'IsoTkMu24',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'SingleMuon',
        ]
        result = self.checkTrigger(*datasets,**triggerNames)
        return result

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['z1','z2']]
        triggerList = ['Mu17_Mu8'] if self.version=='76X' else ['Mu17Mu8']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)







def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy',version='80XMuMuTauTauZSkim'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='muMuTauFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    muMuTauFakeRateAnalysis = MuMuTauFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MuMuTauFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       muMuTauFakeRateAnalysis.analyze()
       muMuTauFakeRateAnalysis.finish()
    except KeyboardInterrupt:
       muMuTauFakeRateAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
