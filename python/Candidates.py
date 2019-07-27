import math

import ROOT

from utilities import deltaR, deltaPhi

##########################
### Basic tree wrapper ###
##########################
class Event(object):
    '''
    Wrap a tree.
    '''
    def __init__(self,tree):
        self.tree = tree

    def __getattr__(self,name):
        return lambda: self.get(name) # returns attribute as a function

    def get(self,var):
        return getattr(self.tree,var)

    def has(self,var):
        return hasattr(self.tree,var)

##############################
### Basic candidate access ###
##############################
class Candidate(object):
    '''
    Encapsulate access to an object in a TTree.
    '''
    def __init__(self,tree,entry=-1,collName=''):
        self.tree = tree
        self.collName = collName
        self.entry = entry

    def __getattr__(self,name):
        return lambda: self.get(name) # returns the attribute as a function

    def __repr__(self):
        return '<{}[{}] {} {:.2f}:{:.2f}:{:.2f}:{:.2f}>'.format(
            self.collName,
            self.entry,
            self.pdgId(),
            self.pt(),
            self.eta(),
            self.phi(),
            self.energy(),
        )

    def get(self,var):
        '''Default variable access from tree.'''
        if not self.tree: return 0
        varName = '{0}_{1}'.format(self.collName,var)
        if self.entry<0: # for flat trees
            return getattr(self.tree,varName)
        else: # for vector branches
            return getattr(self.tree,varName)[self.entry]

    def has(self,var):
        if not self.tree: return False
        varName = '{0}_{1}'.format(self.collName,var)
        return hasattr(self.tree,varName)

    def p4(self):
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiE(self.pt(), self.eta(), self.phi(), self.energy())
        return p4

############
### Muon ###
############
class Muon(Candidate):
    '''
    Muon object access.
    '''
    def __init__(self,tree,entry=-1,collName='muons',shift=None,shiftSigma=0.):
        super(Muon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

    def _p4(self):
        pt = self.get('rochesterPt')
        eta = self.get('eta')
        phi = self.get('phi')
        m = self.get('mass')
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiM(pt,eta,phi,m)
        if self.shiftSigma:
            if pt<100:
                if self.shift=='MuonEnUp' or self.shift=='MuonEnDown': p4 *= 1 + 0.002*self.shiftSigma
            else:
                if self.shift=='MuonEnUp' or self.shift=='MuonEnDown': p4 *= 1 + 0.05*self.shiftSigma
        else:
            if pt<100:
                if self.shift=='MuonEnUp': p4 *= 1+0.002
                if self.shift=='MuonEnDown': p4 *= 1-0.002
            else:
                if self.shift=='MuonEnUp': p4 *= 1+0.05
                if self.shift=='MuonEnDown': p4 *= 1-0.05
        return p4

    def pt(self):
        # instead, because the rochester was applied to the base and not to the energy shifts, calculate manually
        # TODO: need to correct MET
        #var = 'rochesterPt'
        #var = 'pt'
        #if self.shift=='MuonEnUp': var = 'pt_muonEnUp'
        #if self.shift=='MuonEnDown': var = 'pt_muonEnDown'
        #return self.get(var)
        p4 = self._p4()
        return p4.Pt()

    def energy(self):
        #var = 'rochesterEnergy'
        #var = 'energy'
        #if self.shift=='MuonEnUp': var = 'energy_muonEnUp'
        #if self.shift=='MuonEnDown': var = 'energy_muonEnDown'
        #return self.get(var)
        p4 = self._p4()
        return p4.Energy()

################
### Electron ###
################
class Electron(Candidate):
    '''
    Electron object access.
    '''
    def __init__(self,tree,entry=-1,collName='electrons',shift=None,shiftSigma=0.):
        super(Electron, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma =shiftSigma

    def pt(self):
        var = 'pt'
        if self.shift=='ElectronEnUp': var = 'pt_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'pt_electronEnDown'
        return self.get(var)

    def energy(self):
        var = 'energy'
        if self.shift=='ElectronEnUp': var = 'energy_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'energy_electronEnDown'
        return self.get(var)

###########
### Tau ###
###########
class Tau(Candidate):
    '''
    Tau object access.
    '''
    def __init__(self,tree,entry=-1,collName='taus',shift=None,shiftSigma=0.):
        super(Tau, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

    def _p4(self):
        pt = self.get('pt')
        eta = self.get('eta')
        phi = self.get('phi')
        m = self.get('mass')
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiM(pt,eta,phi,m)
        # energy correction (2016 ReReco)
        dm = self.get('decayMode')
        if dm==0:
            p4 *= 1-0.005
        elif dm==1:
            p4 *= 1+0.011
        elif dm==0:
            p4 *= 1+0.006
        # uncertainty
        if pt<400:
            if self.shift=='MuonEnUp': p4 *= 1+0.012
            if self.shift=='MuonEnDown': p4 *= 1-0.012
        else:
            if self.shift=='MuonEnUp': p4 *= 1+0.03
            if self.shift=='MuonEnDown': p4 *= 1-0.03
        return p4

    def pt(self):
        p4 = self._p4()
        return p4.Pt()

    def energy(self):
        p4 = self._p4()
        return p4.Energy()


###########
### Jet ###
###########
class Jet(Candidate):
    '''
    Jet object access.
    '''
    def __init__(self,tree,entry=-1,collName='jets',shift=None,shiftSigma=0.):
        super(Jet, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

    def pt(self):
        var = 'pt'
        if self.shift=='JetEnUp': var = 'pt_jetEnUp'
        if self.shift=='JetEnDown': var = 'pt_jetEnDown'
        return self.get(var)

    def energy(self):
        var = 'energy'
        if self.shift=='JetEnUp': var = 'energy_jetEnUp'
        if self.shift=='JetEnDown': var = 'energy_jetEnDown'
        return self.get(var)

##############
### Photon ###
##############
class Photon(Candidate):
    '''
    Photon object access.
    '''
    def __init__(self,tree,entry=-1,collName='photons',shift=None,shiftSigma=0.):
        super(Photon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

###################
### GenParticle ###
###################
class GenParticle(Candidate):
    '''
    Gen particle object access.
    '''
    def __init__(self,tree,entry=-1,collName='genParticles',shift=None,shiftSigma=0.):
        super(GenParticle, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

###########
### MET ###
###########
class Met(Candidate):
    '''
    Met object access.
    '''
    def __init__(self,tree,entry=0,collName='pfmet',shift=None,shiftSigma=0.):
        super(Met, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift
        self.shiftSigma = shiftSigma

    def __repr__(self):
        return '<{0} {1} {2:.2f}:{3:.2f}>'.format(
            self.collName,
            self.entry,
            self.et(),
            self.phi(),
        )

    def et(self):
        var = 'et'
        if self.shift=='ElectronEnUp':      var = 'et_electronEnUp'
        if self.shift=='ElectronEnDown':    var = 'et_electronEnDown'
        if self.shift=='MuonEnUp':          var = 'et_muonEnUp' # TODO, need to fix to work with Roch, but probably doesn't matter
        if self.shift=='MuonEnDown':        var = 'et_muonEnDown' # TODO
        if self.shift=='TauEnUp':           var = 'et_tauEnUp'
        if self.shift=='TauEnDown':         var = 'et_tauEnDown'
        if self.shift=='PhotonEnUp':        var = 'et_photonEnUp'
        if self.shift=='PhotonEnDown':      var = 'et_photonEnDown'
        if self.shift=='JetEnUp':           var = 'et_jetEnUp'
        if self.shift=='JetEnDown':         var = 'et_jetEnDown'
        if self.shift=='JetResUp':          var = 'et_jetResUp'
        if self.shift=='JetResDown':        var = 'et_jetResDown'
        if self.shift=='UnclusteredEnUp':   var = 'et_unclusteredEnUp'
        if self.shift=='UnclusteredEnDown': var = 'et_unclusteredEnDown'
        return self.get(var)

    def phi(self):
        var = 'phi'
        if self.shift=='ElectronEnUp':      var = 'phi_electronEnUp'
        if self.shift=='ElectronEnDown':    var = 'phi_electronEnDown'
        if self.shift=='MuonEnUp':          var = 'phi_muonEnUp'
        if self.shift=='MuonEnDown':        var = 'phi_muonEnDown'
        if self.shift=='TauEnUp':           var = 'phi_tauEnUp'
        if self.shift=='TauEnDown':         var = 'phi_tauEnDown'
        if self.shift=='PhotonEnUp':        var = 'phi_photonEnUp'
        if self.shift=='PhotonEnDown':      var = 'phi_photonEnDown'
        if self.shift=='JetEnUp':           var = 'phi_jetEnUp'
        if self.shift=='JetEnDown':         var = 'phi_jetEnDown'
        if self.shift=='JetResUp':          var = 'phi_jetResUp'
        if self.shift=='JetResDown':        var = 'phi_jetResDown'
        if self.shift=='UnclusteredEnUp':   var = 'phi_unclusteredEnUp'
        if self.shift=='UnclusteredEnDown': var = 'phi_unclusteredEnDown'
        return self.get(var)

    def p4(self):
        metP4 = ROOT.TLorentzVector()
        metP4.SetPtEtaPhiM(self.et(),0.,self.phi(),0)
        return metP4

############################
### Composite candidates ###
############################
class CompositeCandidate(object):
    '''
    Primary object for access to composite variables.
    '''
    def __init__(self,*objects):
        self.objects = objects

    def __getitem__(self,key):
        if isinstance(key,int): return self.objects[key]

    def __setitem__(self,key,value):
        if isinstance(key,int): self.objects[key] = value

    def __getattr__(self,name):
        try:
            return self.get(name)
        except:
            return self.get(name.capitalize()) # catch things like 'pt' instead of 'Pt'

    def __repr__(self):
        return '<CompositeCandidate {0}>'.format(
            ' '.join([obj.__repr__() for obj in self.objects])
        )

    def get(self,var):
        '''Default variable access from TLorentzVector'''
        vec = self.p4()
        return getattr(vec,var)

    def p4(self):
        p4 = ROOT.TLorentzVector()
        for obj in self.objects: p4 += obj.p4()
        return p4

    def st(self):
        return sum([obj.pt() for obj in self.objects])

###################
### Dicandidate ###
###################
class DiCandidate(CompositeCandidate):
    '''
    Dicandidate variable access.
    '''
    def __init__(self,obj0,obj1):
        super(DiCandidate, self).__init__(obj0,obj1)

    def deltaR(self):
        return deltaR(self.objects[0].eta(),
                      self.objects[0].phi(),
                      self.objects[1].eta(),
                      self.objects[1].phi())

    def deltaPhi(self):
        return deltaPhi(self.objects[0].phi(),
                       self.objects[1].phi())

    def deltaEta(self):
        return abs(self.objects[0].eta()-self.objects[1].eta())

##########################
### Candidate plus met ###
##########################
class MetCompositeCandidate(CompositeCandidate):
    '''
    Met + candidate variable specialization.
    '''
    def __init__(self,met,*objects):
        if not isinstance(met,Met):
            raise TypeError('First argument must be Met object')
        super(MetCompositeCandidate, self).__init__(met,*objects)
        self.met = met
        self.cands = objects

    def metP4(self):
        return self.met.p4()

    def candP4(self):
        p4 = ROOT.TLorentzVector()
        for cand in self.cands: p4 += cand.p4()
        return p4

    def mt(self):
        metP4 = self.metP4()
        candP4 = self.candP4()
        return math.sqrt(abs((candP4.Et()+metP4.Et())**2 - ((candP4+metP4).Pt())**2))

    def mt2(self):
        return self.mt()**2

    def Mt(self):
        return self.mt()

    def Mt2(self):
        return self.mt2()

    def deltaPhi(self):
        candP4 = self.candP4()
        return deltaPhi(self.met.phi(),candP4.Phi())

    def _x(self,c,o):
        dphi = deltaPhi(c.phi(),self.met.phi())
        codphi = deltaPhi(c.phi(),o.phi())
        if abs(codphi)==math.pi or codphi==0:
            x = 0
        else:
            x = c.pt()/(c.pt() + self.met.et()*(math.cos(dphi) - math.sin(dphi)/math.tan(codphi)))
        # note: large met and small deltaphi between c/o results in large negative values for denom
        # this doesnt work for our boosted topology in a->tt since met resolution is poor
        #if x<0:
        #    print x, c.pt(), self.met.et(), dphi, codphi
        return x

    def mcat(self,i1=0,i2=1):
        c1 = self.cands[i1]
        c2 = self.cands[i2]
        x1 = self._x(c1,c2)
        x2 = self._x(c2,c1)
        if x1*x2<=0: return 0.
        cp4 = ROOT.TLorentzVector()
        cp4 += self.cands[i1].p4()
        cp4 += self.cands[i2].p4()
        if len(self.cands)==2: # both candidates had neutrinos
            mcoll = cp4.M() / math.sqrt(x1*x2)
        else: # add in the other cands (assuming no contribution from met)
            op4 = ROOT.TLorentzVector()
            for i,c in enumerate(self.cands):
                if i in [i1,i2]: continue
                op4 += c.p4()
            ncp4 = ROOT.TLorentzVector()
            ncp4.SetPtEtaPhiM(cp4.Pt()/(x1*x2), cp4.Eta(), cp4.Phi(), cp4.M() / math.sqrt(x1*x2))
            tcp4 = ncp4+op4
            mcoll = tcp4.M()
        return mcoll

    def Mcat(self,*args):
        return self.mcat(*args)
