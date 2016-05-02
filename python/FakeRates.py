import os
import sys
import logging
import math

import ROOT

import operator

class FakeRates(object):
    '''Class to access the fakerates for a given lepton ID.'''

    def __init__(self):
        self.fakehists = {'electrons':{},'muons':{},'taus':{}}
        self.fakehists_mc = {'electrons':{},'muons':{},'taus':{}}
        self.fakekey = '{num}_{denom}'
        # WZ fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        self.fake_rootfile = ROOT.TFile(fake_path)
        self.fakehists['electrons'][self.fakekey.format(num='WZMedium',denom='WZLoose')] = self.fake_rootfile.Get('e/medium/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='WZTight',denom='WZLoose')] = self.fake_rootfile.Get('e/tight/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='WZMedium',denom='WZLoose')] = self.fake_rootfile.Get('m/medium/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='WZTight',denom='WZLoose')] = self.fake_rootfile.Get('m/tight/fakeratePtEta')
        # H++ fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_hpp_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        self.fake_hpp_rootfile = ROOT.TFile(fake_path)
        self.fakehists['electrons'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_hpp_rootfile.Get('e/medium/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_hpp_rootfile.Get('e/tight/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_hpp_rootfile.Get('m/medium/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_hpp_rootfile.Get('m/tight/fakeratePtEta')
        # tau fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_w_tau_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        self.fake_tau_rootfile = ROOT.TFile(fake_path)
        self.fakehists['taus'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta_fromMC')

    def __exit__(self, type, value, traceback):
        self.__finish()

    def __del__(self):
        self.__finish()

    def __finish(self):
        self.fake_rootfile.Close()
        self.fake_hpp_rootfile.Close()
        self.fake_tau_rootfile.Close()

    def __get_fakerate(self,cand,pt,eta,num,denom):
        if cand[0] not in self.fakehists: return 0.
        key = self.fakekey.format(num=num,denom=denom)
        if key not in self.fakehists[cand[0]]: return 0.
        hist = self.fakehists[cand[0]][key]
        if pt > 100.: pt = 99.
        return hist.GetBinContent(hist.FindBin(pt,abs(eta)))

    def __get_fakerate_mc(self,cand,pt,eta,num,denom):
        if cand[0] not in self.fakehists_mc: return 0.
        key = self.fakekey.format(num=num,denom=denom)
        if key not in self.fakehists_mc[cand[0]]: return 0.
        hist = self.fakehists_mc[cand[0]][key]
        if pt > 100.: pt = 99.
        return hist.GetBinContent(hist.FindBin(pt,abs(eta)))

    def getFakeRate(self,rtrow,cand,num,denom):
        if cand[1]<0: return 0 # not defined
        pt  = getattr(rtrow,'{0}_rochesterPt'.format(cand[0]))[cand[1]] if cand[0]=='muons' else getattr(rtrow,'{0}_pt'.format(cand[0]))[cand[1]]
        eta = getattr(rtrow,'{0}_rochesterEta'.format(cand[0]))[cand[1]] if cand[0]=='muons' else getattr(rtrow,'{0}_eta'.format(cand[0]))[cand[1]]
        return self.__get_fakerate(cand,pt,eta,num,denom)

    def getFakeRateMC(self,rtrow,cand,num,denom):
        if cand[1]<0: return 0 # not defined
        pt  = getattr(rtrow,'{0}_rochesterPt'.format(cand[0]))[cand[1]] if cand[0]=='muons' else getattr(rtrow,'{0}_pt'.format(cand[0]))[cand[1]]
        eta = getattr(rtrow,'{0}_rochesterEta'.format(cand[0]))[cand[1]] if cand[0]=='muons' else getattr(rtrow,'{0}_eta'.format(cand[0]))[cand[1]]
        return self.__get_fakerate_mc(cand,pt,eta,num,denom)

