import os

from DevTools.Utilities.hdfsUtils import *
from DevTools.Utilities.utilities import *

import ROOT

job = {
    'vbf': 'crab_2018-11-20_Skim_MuMuTauTau_VBF_80X_v1',
    'gg' : 'crab_2018-11-20_Skim_MuMuTauTau_GGF_80X_v1',
}

amasses = {
    'vbf': [5,9,15,21],
    'gg' : [5,9,15],
}
hmasses = {
    'vbf': [125, 300, 750],
    'gg' : [200,250,400,500,1000],
}

for m in ['vbf','gg']:
    tag = 'GluGlu' if m=='gg' else 'VBF'
    samples = []
    for h in hmasses[m]:
        for a in amasses[m]:
            if h==125:
               samples += ['SUSY{tag}ToHToAA_AToMuMu_AToTauTau_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'.format(tag=tag,a=a)]
            else:
               samples += ['SUSY{tag}ToHToAA_AToMuMu_AToTauTau_M-{h}_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'.format(tag=tag,h=h,a=a)]
    
    for sample in samples:
        baseDir = '/hdfs/store/user/dntaylor/{}/{}'.format(job[m],sample)
        files = get_hdfs_root_files(baseDir)
        files = [f for f in files if 'lumi' in f]
    
        tree = ROOT.TChain('lumiTree/LumiTree')
        for fName in files:
            if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs',fName)
            tree.Add(fName)
    
    
        summedWeights = {-1: 0}
        for row in tree:
            nw = len(row.summedGenWeights)
            summedWeights[-1] += row.summedWeights
            for i in range(nw):
                if i not in summedWeights: summedWeights[i] = 0
                summedWeights[i] += row.summedGenWeights[i]
    
        dumpResults(summedWeights,'MuMuTauTau','{}/summedWeights'.format(sample))
