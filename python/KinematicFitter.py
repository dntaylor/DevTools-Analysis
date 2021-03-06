import os
import sys
import math
from array import array
from collections import OrderedDict
from scipy.optimize import minimize, minimize_scalar

import ROOT

ROOT.gROOT.SetBatch()

from DevTools.Analyzer.utilities import deltaR, deltaPhi


class COMMON(object):

    @staticmethod
    def printP4(label,p4,extra=''):
        print '{}: pt {:6.2f}, eta {:6.2f}, phi {:6.2f}, energy {:6.2f}, mass {:6.2f}{}'.format(
            label, p4.Pt(), p4.Eta(), p4.Phi(), p4.E(), p4.M(), ', {0}'.format(extra) if extra else '',
        )

    @staticmethod
    def printCov(label,cov,extra=''):
        print '{}: cov [{:.2f}, {:.2f}; {:.2f}, {:.2f}]{}'.format(
            label, cov[0][0], cov[0][1], cov[1][0], cov[1][1], ', {0}'.format(extra) if extra else '',
        )

    @staticmethod
    def printVec(label,vec,extra=''):
        print '{}: pt {:6.2f}, eta {:6.2f}, phi {:6.2f}{}'.format(
            label, vec.Pt(), vec.Eta(), vec.Phi(), ', {0}'.format(extra) if extra else '',
        )


class PDG(object):

    M_E = 0.000511
    M_MU = 0.105658
    M_TAU = 1.77686

    M_PI0 = 0.13957
    M_PIPLUS = 0.1344977
    M_RHO = 0.77526
    M_A1 = 1.260


class Particle(ROOT.TLorentzVector):

    def __init__(self,label,pt=None,eta=None,phi=None,energy=None,mass=None,px=None,py=None,pz=None):
        super(Particle, self).__init__()
        self.label = label
        self.status = -1
        if None not in [pt,eta,phi,energy]:
            self.SetPtEtaPhiE(pt,eta,phi,energy)
        elif None not in [pt,eta,phi,mass]:
            self.SetPtEtaPhiM(pt,eta,phi,mass)
        elif None not in [px,py,pz,energy]:
            self.SetPxPyPzE(px,py,pz,energy)
        else:
            print 'Cannot initialize particle'
            print pt, eta, phi, energy, mass, px, py, pz
        self.recoP4 = ROOT.TLorentzVector(self)
        self.dE = 0.
        self.dEta = 0.
        self.dPhi = 0.

    def setErrors(self,dE,dEta,dPhi):
        self.dE = dE
        self.dEta = dEta
        self.dPhi = dPhi

    def getStatus(self):
        return self.status

    def getPDGMass(self):
        return 0

    def getCov(self):
        dp = self.dE * self.P()/self.E()
        dpt = math.sin(self.Theta()) * dp
        cov = ROOT.TMatrixD(2,2)
        cov[0][0] = (math.cos(self.Phi())*dpt)**2
        cov[0][1] = math.sin(self.Phi())*math.cos(self.Phi())*dpt*dpt
        cov[1][0] = math.sin(self.Phi())*math.cos(self.Phi())*dpt*dpt
        cov[1][1] = (math.sin(self.Phi())*dpt)**2
        return cov

    def getChi2(self):
        dE = self.dE
        if not dE: return 0
        E_fit = self.E()
        E_reco = self.getRecoP4().E()
        return (E_fit-E_reco)**2 / dE**2

    def getRecoP4(self):
        return ROOT.TLorentzVector(self.recoP4)

    def scaleEnergy(self,X):
        valid = True
        if X<=0:
            #print 'Warning, Invalid X {0}, will not scale'.format(X)
            valid = False
        # scale energy assuming that X is the visible fraction of total energy
        beta = self.recoP4.Pt()/self.recoP4.E()
        newE = self.recoP4.E()/X if X else self.recoP4.E()
        newPt = beta*newE
        newEta = self.recoP4.Eta()
        newPhi = self.recoP4.Phi()
        self.SetPtEtaPhiE(newPt,newEta,newPhi,newE)
        self.status = int(valid)
        return valid

    def printP4(self):
        COMMON.printP4(self.label,self)

    def printRecoP4(self):
        COMMON.printP4(self.label + ' RECO', self.recoP4)

    def printCov(self):
        COMMON.printCov(self.label, self.getCov())

        
class Electron(Particle):

    def __init__(self,label,pt,eta,phi,energy):
        super(Electron, self).__init__(label,pt=pt,eta=eta,phi=phi,energy=energy)

    def getPDGMass(self):
        return PDG.M_E

class Muon(Particle):

    def __init__(self,label,pt,eta,phi,energy):
        super(Muon, self).__init__(label,pt=pt,eta=eta,phi=phi,energy=energy)

    def getPDGMass(self):
        return PDG.M_Mu

class Tau(Particle):

    def __init__(self,label,pt,eta,phi,energy,dm):
        super(Tau, self).__init__(label,pt=pt,eta=eta,phi=phi,energy=energy)
        self.decayMode = dm

    def getPDGMass(self):
        return PDG.M_TAU

    def getDecayMode(self):
        return self.decayMode

    def printP4(self):
        COMMON.printP4(self.label,self,extra='dm {0}'.format(self.decayMode))

    def printRecoP4(self):
        COMMON.printP4(self.label + ' RECO', self.recoP4,extra='dm {0}'.format(self.decayMode))


class HadronicTau(Tau):

    def __init__(self,label,pt,eta,phi,energy,dm):
        super(HadronicTau, self).__init__(label,pt,eta,phi,energy,dm)

class ElectronTau(Tau):

    def __init__(self,label,pt,eta,phi,energy):
        super(ElectronTau, self).__init__(label,pt,eta,phi,energy,dm=-11)

class MuonTau(Tau):

    def __init__(self,label,pt,eta,phi,energy):
        super(MuonTau, self).__init__(label,pt,eta,phi,energy,dm=-13)

class MET(Particle):

    def __init__(self,label,pt,phi,cov00,cov01,cov10,cov11):
        super(MET, self).__init__(label,pt=pt,eta=0,phi=phi,mass=0)
        self.cov = ROOT.TMatrixD(2,2)
        self.cov[0][0] = cov00
        self.cov[0][1] = cov01
        self.cov[1][0] = cov10
        self.cov[1][1] = cov11

    def getCov(self):
        return ROOT.TMatrixD(self.cov)

class Composite(tuple):

    def __init__(self, label, *particles):
        self.label = label
        self.uncertainty = 0 # uncertainty on mass
        self.shift = 0 # the sigma shift on the uncertainty
        self.status = -1
        self._updateP4()
        self.recoP4 = ROOT.TLorentzVector()
        for p in self:
            self.recoP4 += p.getRecoP4()

    def __new__(self, label, *particles):
        return tuple.__new__(self, tuple(particles))

    def __getattr__(self,name):
        self._updateP4() # TODO: Find a way to see if the underlying TLorentzVectors have been updated
        return lambda: getattr(self.p4,name)() # return attribute as a function

    def __contains__(self,item):
        return item in [p.label for p in self]

    def _updateP4(self):
        self.p4 = ROOT.TLorentzVector()
        for p in self:
            self.p4 += p

    def getStatus(self):
        return self.status

    def getRecoP4(self):
        return ROOT.TLorentzVector(self.recoP4)

    def getDeltaR(self):
        if len(self)!=2: return 0
        ap4 = self[0]
        bp4 = self[1]
        return deltaR(ap4.Eta(),ap4.Phi(),bp4.Eta(),bp4.Phi())

    def getRecoDeltaR(self):
        if len(self)!=2: return 0
        ap4 = self[0].getRecoP4()
        bp4 = self[1].getRecoP4()
        return deltaR(ap4.Eta(),ap4.Phi(),bp4.Eta(),bp4.Phi())

    def setUncertainty(self,uncertainty):
        self.uncertainty = uncertainty

    def setShift(self,shift):
        self.shift = shift

    def M(self):
        # hackish to account for uncertainty on composite object rather than individual components
        self._updateP4()
        if self.uncertainty:
            return self.p4.M()*(1+self.shift*self.uncertainty)
        else:
            return self.p4.M()

    def getChi2(self):
        if not self.uncertainty: return 0
        return self.shift**2 / self.uncertainty**2

    def getCosAlpha(self):
        cosa = math.cos(self[0].Angle(self[1].Vect()))
        return cosa

    def constrainEnergyFromMass(self,fixed,mass):
        valid = True
        if len(self)!=2:
            print 'Don\'t know how to handle constraints with more than 2 particles'
            valid = False
            self.status = int(valid)
            return valid

        a = self[0] if self[0].label==fixed else self[1]
        b = self[1] if self[0].label==fixed else self[0]

        # get new energy given constraints
        #   a : particle that was updated
        #   b : particle to constrain
        #   mass : mass constraint
        am = a.getPDGMass()
        bm = b.getPDGMass()
        cosa = math.cos(a.Angle(b.Vect()))
        D = a.E()/(a.P()*cosa)
        F = (mass**2 - 2*am**2)/(2*a.P()*cosa)
        G = am**2 * (1-D**2) + F**2
        if G<0:
            #print 'Warning, negative value in sqrt, will not constrain energy'
            #print a.E(), a.P(), am, cosa, mass
            #print D, F, G
            valid = False
            E = b.E()
        elif cosa>0:
            E = 1/(D**2-1) * (D*F + (G)**0.5)
        elif cosa<0:
            E = 1/(D**2-1) * (D*F - (G)**0.5)
        else:
            valid = False
            E = 0.

        # update energy
        beta = b.getRecoP4().Pt()/b.getRecoP4().E()
        #b.SetPtEtaPhiM(beta*E, b.getRecoP4().Eta(), b.getRecoP4().Phi(), b.getPDGMass())
        b.SetPtEtaPhiE(beta*E, b.getRecoP4().Eta(), b.getRecoP4().Phi(), E)

        self.status = int(valid)
        return valid


    def printP4(self):
        COMMON.printP4(self.label,self,extra='dR {:.2f}'.format(self.getDeltaR()) if len(self)==2 else '')

    def printRecoP4(self):
        COMMON.printP4(self.label + ' RECO', self.recoP4, extra='dR {:.2f}'.format(self.getRecoDeltaR()) if len(self)==2 else '')


class KinematicConstraints(object):

    def __init__(self):
        self.particles = OrderedDict()
        self.genparticles = OrderedDict()
        self.composites = OrderedDict()
        self.constraints = []

    def addParticle(self,label,particle):
        self.particles[label] = particle

    def addGenParticle(self,label,particle):
        self.genparticles[label] = particle

    def addComposite(self,label,*plabels, **kwargs):
        self.composites[label] = Composite(label, *[self.particles[p] for p in plabels])
        self.composites[label].setUncertainty(kwargs.pop('uncertainty',0))

    def getParticleStatus(self,label):
        return self.particles[label].getStatus()

    def getCompositeStatus(self,label):
        return self.composites[label].getStatus()

    def getParticle(self,label):
        return self.particles[label]

    def getComposite(self,label):
        return self.composites[label]

    def getRecoil(self):
        recoil = ROOT.TVector3()
        for p,part in self.particles.iteritems():
            recoil -= part.Vect()
        return recoil

    def getRecoRecoil(self):
        recoil = ROOT.TVector3()
        for p,part in self.particles.iteritems():
            recoil -= part.getRecoP4().Vect()
        return recoil

    def getRecoilCov(self):
        covrecoil = ROOT.TMatrixD(2,2)
        for p in self.particles:
            part = self.particles[p]
            if isinstance(part,MET):
                covrecoil += part.getCov()
            else:
                # simplification ignores this
                #covrecoil -= part.getCov()
                pass
        return covrecoil

    def getDeltaRecoil(self):
        drecoil = self.getRecoRecoil() - self.getRecoil()
        return drecoil

    def getChi2(self):
        covrecoil = self.getRecoilCov()
        covinv = covrecoil.Invert()
        recoildiff = self.getDeltaRecoil()
        rpx = recoildiff.Px()
        rpy = recoildiff.Py()
        recoildiffvec = ROOT.TVectorD(2)
        recoildiffvec[0] = rpx
        recoildiffvec[1] = rpy
        #recoilchi2 = recoildiffvec * (covinv * recoildiffvec) # not working
        recoilchi2 = rpx * (covinv[0][0] * rpx + covinv[0][1] * rpy) + rpy * (covinv[1][0] * rpx + covinv[1][1] * rpy)

        # add gaussian particle chi2 = (x_fit-x_reco)**2/sigma**2
        # will be used if we add in the muon pt resolution and shift within its uncertainty
        #hasUncE = [p for p in self.particles if self.particles[p].dE]
        #for p in hasUncE:
        #    recoilchi2 += self.particles[p].getChi2()

        hasUnc = [c for c in self.composites if self.composites[c].uncertainty]
        for c in hasUnc:
            recoilchi2 += self.composites[c].getChi2()

        return recoilchi2

    def getP4(self,*labels):
        p4 = ROOT.TLorentzVector()
        for p in labels:
            p4 += self.particles[p]
        return p4

    def getDeltaR(self,composite):
        return self.composites[composite].deltaR()

    def getDeltaR(self,a,b):
        ap4 = self.getP4(a)
        bp4 = self.getP4(b)
        return deltaR(ap4.Eta(),ap4.Phi(),bp4.Eta(),bp4.Phi())

    def addMassConstraint(self,particles1,particles2,mass=0):
        # require p4(particles1).M() - p4(particles2).M() - mass = 0
        self.constraints += [('mass',particles1,particles2,mass)]

    def getAllowedEnergyFractionRange(self,label):
        # TODO scale within uncertainties as well
        # scale of visible energy
        part = self.particles[label]
        if isinstance(part,ElectronTau) or isinstance(part,MuonTau):
            # leptonic decay of tau
            #print 'ElectronTau, MuonTau, unconstrained'
            return [0,1]
        elif isinstance(part,HadronicTau):
            mvis = part.getRecoP4().M()
            #if mvis==0:
            #    print 'Tau visibile mass is 0'
            #if mvis>PDG.M_TAU:
            #    print 'Visible mass greater than tau mass: {0}'.format(mvis)
            #if part.getDecayMode() in [1]:
            #    # h+pi0 has a wider window due to resolution effects
            #    # this is taken from the tau reco paper on the window of the tau mass
            #    lower = 0.3
            #    upper = 1.2*(part.getRecoP4().Pt()/100)
            #    if upper<1.2: upper = 1.2
            #    if upper>4.2: upper = 4.2
            #    return [lower**2/PDG.M_TAU**2, upper**2/PDG.M_TAU**2]
            #else:
            return [mvis**2/PDG.M_TAU**2,1]
        else:
            return [1,1]

    def getAlternativeAllowedEnergyFractionRange(self,label,constraint='mass'):
        part = self.particles[label]
        con = [c for c in self.constraints if c[0]==constraint][0]
        other = self.composites[con[1]] if label in self.composites[con[2]] else self.composites[con[2]]
        this  = self.composites[con[2]] if label in self.composites[con[2]] else self.composites[con[1]]
        mass = other.M()
        M = part.getPDGMass()
        E = part.getRecoP4().E()
        P = part.getRecoP4().P()
        B = P/E
        cosa = this.getCosAlpha()
        minX = M * (1/(B*cosa)**2 - 1)**0.5 / abs((mass**2 - 2*M**2)/(2*P*cosa))
        return [minX,1]

    def scaleParticleEnergy(self,label,X):
        part = self.particles[label]
        valid = part.scaleEnergy(X)
        # update other particles using constraints
        for c in self.constraints:
            if c[0]=='mass':
                c1 = self.composites[c[1]]
                c2 = self.composites[c[2]]
                mass = c[3]
                if label in c2:
                    tmp = c2
                    c2 = c1
                    c1 = tmp
                mass = mass + c2.shift # shift by the current uncertainty
                valid = c1.constrainEnergyFromMass(label,c2.M()-mass) and valid
        return valid

    def scaleParticleUncertainty(self,label,unc):
        part = self.particles[label]
        X = 1+(unc/part.getRecoP4().E())
        valid = part.scaleEnergy(X)
        return valid

    def printRecoil(self):
        COMMON.printVec('recoil', self.getRecoil())

    def printRecoRecoil(self):
        COMMON.printVec('recoil RECO', self.getRecoRecoil())

    def printRecoilCov(self):
        COMMON.printCov('recoil', self.getRecoilCov())

    def printKinematics(self):
        print 'Fitted Kinematics'
        for p,p4 in self.particles.iteritems():
            p4.printP4()
        for p,p4 in self.composites.iteritems():
            p4.printP4()
        self.printRecoil()
        print 'chi2 {:.2f}'.format(self.getChi2())

    def printRecoKinematics(self):
        print 'RECO kinematics'
        for p,p4 in self.particles.iteritems():
            p4.printRecoP4()
            #p4.printCov()
        for p,p4 in self.composites.iteritems():
            p4.printRecoP4()
        self.printRecoRecoil()
        #self.printRecoilCov()

    def printGenKinematics(self):
        if not self.genparticles: return
        print 'GEN kinematics'
        for p,p4 in self.genparticles.iteritems():
            p4.printRecoP4()

class KinematicFitter(KinematicConstraints):

    def __init__(self):
        super(KinematicFitter, self).__init__()
        self.minimizationParticles = []
        self.fitStatus = -1
        self.X = -1
        self.chi2 = -1
        self.atLowerBound = -1
        self.atUpperBound = -1

    def setMinimizationParticles(self,*labels):
        self.minimizationParticles = labels

    # attempt to use root
    #def build_func(self):

    #    other = self

    #    class Functor(ROOT.TPyMultiGenFunction):

    #        def __init__(self):
    #            super(Functor, self).__init__()

    #        def NDim(self):
    #            return 1

    #        def DoEval(self, args):
    #            for i in range(len(other.minimizationParticles)):
    #                other.scaleParticleEnergy(other.minimizationParticles[i],args[i])
    #            chi2 = other.getChi2()
    #            return chi2

    #    func = Functor()
    #    return func
    #

    def getFitStatus(self):
        return self.fitStatus

    def getFitX(self):
        return self.X

    def getFitChi2(self):
        return self.chi2

    def fit(self,useUncertainty=False):
        if len(self.minimizationParticles)!=1:
            print 'Can only fit 1 particle currently'
            return

        #useUncertainty = False

        # get allowed range of first fit variable
        X_range = self.getAllowedEnergyFractionRange(self.minimizationParticles[0])
        X_range_alt = self.getAlternativeAllowedEnergyFractionRange(self.minimizationParticles[0])
        #print 'ranges', X_range, X_range_alt
        X_range = X_range_alt
        if X_range[0]>X_range[1]: X_range = [X_range[1], X_range[0]]
        X = X_range[0] + (X_range[1]-X_range[0])/2 # split the difference as a test
        X = 1

        # this was for scaling an individual particles uncertainty
        #hasUncE   = [p for p in self.particles if self.particles[p].dE]
        ##hasUncEta = [p for p in self.particles if self.particles[p].dEta]
        ##hasUncPhi = [p for p in self.particles if self.particles[p].dPhi]
        #boundsE   = {p: [self.particles[p].dE*x   for x in [-5,5]] for p in hasUncE}
        ##boundsEta = {p: [self.particles[p].dEta*x for x in [-5,5]] for p in hasUncEta}
        ##boundsPhi = {p: [self.particles[p].dPhi*x for x in [-5,5]] for p in hasUncPhi}

        #X_unc_range = [X_range] + [boundsE[p] for p in hasUncE]

        # composites
        hasUnc = [c for c in self.composites if self.composites[c].uncertainty]
        bounds = {c: [self.composites[c].uncertainty*x for x in [-5,5]] for c in hasUnc}
        X_unc_range = [X_range] + [bounds[c] for c in hasUnc]

        # quick test
        #self.fit_func(X)

        # run fit
        # attempt to use root
        #minimizer = ROOT.Math.Factory.CreateMinimizer('GSLMultiMin')
        #minimizer.SetMaxFunctionCalls(1000)
        #minimizer.SetMaxIterations(1000)
        #minimizer.SetTolerance(1e-3)
        #func = self.build_func()
        #minimizer.SetFunction(func)
        #for p,label in enumerate(self.minimizationParticles):
        #    minimizer.SetVariable(p, label, X, 1e-3)
        #minimizer.Minimize()

        def func_unc(xi):
            X = xi[0]
            shift = xi[1:]
            #hasUncE   = [p for p in self.particles if self.particles[p].dE]
            #hasUncEta = [p for p in self.particles if self.particles[p].dEta]
            #hasUncPhi = [p for p in self.particles if self.particles[p].dPhi]
            #hasUnc = hasUncE + hasUncEta + hasUncPhi
            #hasUnc = hasUncE
            hasUnc = [c for c in self.composites if self.composites[c].uncertainty]
            valid = True
            #for s,p in zip(shift,hasUnc):
            #    valid = self.scaleParticleUncertainty(p,s) and valid
            for s,c in zip(shift,hasUnc):
                self.composites[c].setShift(s)
            valid = self.scaleParticleEnergy(self.minimizationParticles[0],X) and valid
            chi2 = self.getChi2()
            self.fitStatus = int(valid)
            if not valid: chi2 = 9999
            return chi2

        def func(X):
            valid = self.scaleParticleEnergy(self.minimizationParticles[0],X)
            chi2 = self.getChi2()
            self.fitStatus = int(valid)
            if not valid: chi2 = 9999
            return chi2

        # there is a problem with the minimization right now
        # there are invalid regions, and the function returns 9999, and cant find the derivative to improve
        if useUncertainty:
            # individual particles
            #x0 = [X] + [0]*len(hasUncE)
            # composite particles
            x0 = [X] + [0] * len(hasUnc)
            res = minimize(func_unc, x0, bounds=X_unc_range)
            self.X = res.x[0]
            self.atLowerBound = int(abs(self.X-X_range[0])<1e-3)
            self.atUpperBound = int(abs(self.X-X_range[1])<1e-3)
            self.chi2 = func_unc(res.x)
            if res.x[0]<X_range[0] or res.x[0]>X_range[1]:
                print 'result out of range', res.x, self.chi2, X_range
        else:
            #res = minimize_scalar(func, bounds=X_range, method='brent')
            res = minimize_scalar(func, bounds=X_range, method='bounded')
            self.X = res.x
            self.atLowerBound = int(abs(self.X-X_range[0])<1e-3)
            self.atUpperBound = int(abs(self.X-X_range[1])<1e-3)
            self.chi2 = func(res.x)
            if res.x<X_range[0] or res.x>X_range[1]:
                print 'result out of range', res.x, self.chi2, X_range
        #if not self.fitStatus:
        #    print 'Non valid fit'
        #    self.printRecoKinematics()
        #    self.printKinematics()
        #    self.printGenKinematics()
        #    print ''
        #part = self.particles['th']
        #mvis = part.getRecoP4().M()
        #if mvis>PDG.M_TAU:
        #    self.printRecoKinematics()
        #    self.printKinematics()
        #    self.printGenKinematics()
        #    print ''
        #print self.X, self.chi2
            

    def getGrid(self,**kwargs):
        '''Plot grid values'''
        npoints = kwargs.pop('npoints',100)
        useUncertainty = kwargs.pop('useUncertainty',False)

        X_range = self.getAllowedEnergyFractionRange(self.minimizationParticles[0])
        X_range_alt = self.getAlternativeAllowedEnergyFractionRange(self.minimizationParticles[0])
        X_range = X_range_alt

        #hasUncE   = [p for p in self.particles if self.particles[p].dE]
        ##hasUncEta = [p for p in self.particles if self.particles[p].dEta]
        ##hasUncPhi = [p for p in self.particles if self.particles[p].dPhi]
        #boundsE   = {p: [self.particles[p].dE*x   for x in [-5,5]] for p in hasUncE}
        ##boundsEta = {p: [self.particles[p].dEta*x for x in [-5,5]] for p in hasUncEta}
        ##boundsPhi = {p: [self.particles[p].dPhi*x for x in [-5,5]] for p in hasUncPhi}

        #X_unc_range = [boundsE[p] for p in hasUncE]

        # composites
        hasUnc = [c for c in self.composites if self.composites[c].uncertainty]
        bounds = {c: [self.composites[c].uncertainty*x for x in [-5,5]] for c in hasUnc}
        X_unc_range = [bounds[c] for c in hasUnc]

        xvals = []
        chi2vals = []
        for i in range(npoints):
            X = X_range[0] + (X_range[1]-X_range[0])*i/(npoints-1)
            if useUncertainty:
                def func_unc(xi):
                    shift = xi
                    #hasUncE   = [p for p in self.particles if self.particles[p].dE]
                    #hasUncEta = [p for p in self.particles if self.particles[p].dEta]
                    #hasUncPhi = [p for p in self.particles if self.particles[p].dPhi]
                    #hasUnc = hasUncE + hasUncEta + hasUncPhi
                    #hasUnc = hasUncE
                    hasUnc = [c for c in self.composites if self.composites[c].uncertainty]
                    valid = True
                    #for s,p in zip(shift,hasUnc):
                    #    valid = self.scaleParticleUncertainty(p,s) and valid
                    for s,c in zip(shift,hasUnc):
                        self.composites[c].setShift(s)
                    valid = self.scaleParticleEnergy(self.minimizationParticles[0],X) and valid
                    chi2 = self.getChi2()
                    self.fitStatus = int(valid)
                    if not valid: chi2 = 9999
                    return chi2

                x0 = [0]*len(hasUnc)
                res = minimize(func_unc, x0, bounds=X_unc_range)
                chi2 = func_unc(res.x)
            else:
                valid = self.scaleParticleEnergy(self.minimizationParticles[0],X)
                chi2 = self.getChi2()
                if not valid: chi2 = 9999
            xvals += [X]
            chi2vals += [chi2]

        return xvals, chi2vals

