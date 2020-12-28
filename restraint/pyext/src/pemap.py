#!/usr/bin/env python
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.isd
import IMP.container
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.restraints
from math import log
from collections import defaultdict
import itertools
import operator
import os
import math
#IMP.pmi.restraints.RestraintBase
class PEMAPRestraint(IMP.pmi.restraints.RestraintBase):
    
    def __init__(self,
                 root_hier,
                 mic_file,
                 sigma_init = 4.0,
                 weight=1.0,
                 slope = 0.0,
                 label = None):

        """
        Constructor:
        @param root_hier         Hierarchy of the system
        @param mic_file          File containing p1,r1,p2,r2,MIC
        @param weight            Weight of the restraint
        """

        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        
        rname = "PEMAPRestraint"
        super(PEMAPRestraint, self).__init__(
            self.m, name="PEMAPRestraint", label=label, weight=weight)

        self.weight = weight
        self.slope = slope
        self.sigma_init = sigma_init
        
        
        # Nuisance particles
        self.sigma_dictionary={}
        self.psi_dictionary={}
        
        self.sigma_is_sampled = True
        self.rs_sig = self._create_restraint_set("sigma")
        self._create_sigma('sigma_0', self.sigma_init)
        sigma_0 = self.sigma_dictionary['sigma_0'][0].get_particle_index()
        
        #self.psi_is_sampled = True
        #self.rs_psi = self._create_restraint_set("psi")
        # Create particle
        #psi_init = 0.01
        #self._create_psi('psi_0', psi_init)
        #psip_0 = self.psi_dictionary['psi_0'][0].get_particle_index()

        
        
        # Pairs of restrainted beads 
        self.pairs = []
       
        # Setup restraint
        self.pmr = IMP.isd.PEMAPRestraint(self.m,
                                          sigma_0)

        
        self.rs = IMP.RestraintSet(self.m, 'PEMAPRestraint')
        
        # Read file
        for line in open(mic_file):
            vals = line.split()
            if (vals[0]=="#"):
                continue
            prot1 = vals[0]
            resi1 = int(vals[2])
            prot2 = vals[1]
            resi2 = int(vals[3])
           
            mic = float(vals[4])
            if len(vals)==6:
                dist = float(vals[5])
            else:
                dist = 0.0
                
            # Set restraint paramters
            d_mic = self.get_d_mic(mic)
        
            ########################
            # Select particles
            ########################
            try:
                s0 = IMP.atom.Selection(root_hier,
                                        molecule = prot1,
                                        residue_index = resi1)
                s1 = IMP.atom.Selection(root_hier,
                                      molecule = prot2, 
                                      residue_index=resi2)
                
                p0 = s0.get_selected_particles()[0]
                p1 = s1.get_selected_particles()[0]

                ########################
                # Setup restraint
                ########################
                self.pmr.add_contribution(p0.get_index(), p1.get_index(), d_mic)
                self.pairs.append((p0, resi1, p1, resi2, dist))
            except:
                print("SimplifiedPEMAP: WARNING> residue %d of chain %s is not there (w/ %d %s)" % (resi1,prot1,resi2,prot2))
            
        rname = self.pmr.get_name()
        self.rs.add_restraint(self.pmr)

        self.restraint_sets = [self.rs] + self.restraint_sets[1:]
        
    def get_d_mic(self, mic):
        k = -0.014744
        n = -0.41
        if mic <0.6:
            return (math.log(mic)-n)/k
        else:
            return (math.log(0.6)-n)/k
        
    def get_line_params(self, pemap_file):
        import numpy as np
        import pylab as pl

        data = open(pemap_file)
        D = data.readlines()
        data.close()

        cc,ds = [],[]
        
        for d in D:
            d = d.strip().split()
            c=float(d[4])
            s=float(d[5])
            cc.append(c)
            ds.append(s)
        ccm,dsm=[],[]
        L = np.linspace(min(ds),max(ds),num=11)
        for i in xrange(1,len(L)):
            m = [x for x,j in enumerate(ds) if L[i-1]<=j<L[i]]
            if len(m)==0: continue
            c = max([cc[j] for j in m])
            ccm.append(cc[cc.index(c)])
            dsm.append(ds[cc.index(c)])
        ccm=np.array(ccm[1:])
        dsm=np.array(dsm[1:])
        A = np.vstack([dsm, np.ones(len(dsm))]).T
        k,n = np.linalg.lstsq(A, ccm)[0]
        print("Line parameters are: k=%f and n=%f" % (k,n))

        
        pl.plot(ds,cc,'ko')
        pl.plot(dsm,ccm,'ro')
        pl.plot(L,-0.0075*L+1.,'r-')
        pl.plot(L,k*L+n,'k-')
        pl.show()
        return k,n    

    def get_upper_bound(self,mic,k,n):   
        if mic <= 0.6:
            return (log(mic) - n) / k
        else: return (log(0.6) - n) / k
    
    def get_lower_bond(self,pearsoncc):
        return (pearsoncc-1.)/-0.0551

    def _create_psi(self, name,psiinit):
        """ Creates nuisances on the data uncertainty """

        m = self.root_hier.get_model()
        
        if name in self.psi_dictionary:
            return self.psi_dictionary[name][0]
        psiminnuis = 0.0000001
        psimaxnuis = 0.9999999
        psimin = 0.01
        psimax = 0.99
        psitrans = 0.1
        psi = IMP.pmi.tools.SetupNuisance(m,
                                          psiinit,
                                          psiminnuis,
                                          psimaxnuis,
                                          self.psi_is_sampled).get_particle()
        self.psi_dictionary[name] = (
            psi,
            psitrans,
            self.psi_is_sampled)

        self.rs_psi.add_restraint(IMP.isd.UniformPrior(self.m,
                                                       psi,
                                                       1000000000.0,
                                                       psimax,
                                                       psimin))

        #self.rs_psi.add_restraint(IMP.isd.JeffreysRestraint(self.m, psi))
        return psi
    
    def _create_sigma(self, name, sigmainit):
        """ This is called internally. Creates a nuisance
        on the structural uncertainty """
        if name in self.sigma_dictionary:
            return self.sigma_dictionary[name][0]

        m = self.root_hier.get_model()
        
        sigmaminnuis = 0.0000001
        sigmamaxnuis = 1000.0
        sigmamin = 0.0
        sigmamax = 4.0
        sigmatrans = 0.2
        sigma = IMP.pmi.tools.SetupNuisance(self.m,
                                            sigmainit,
                                            sigmaminnuis,
                                            sigmamaxnuis,
                                            self.sigma_is_sampled).get_particle()
        self.sigma_dictionary[name] = (sigma,
                                       sigmatrans,
                                       self.sigma_is_sampled)
        
        self.rs_sig.add_restraint(IMP.isd.UniformPrior(self.m,
                                                       sigma,
                                                       1000000000.0,
                                                       sigmamax,
                                                       sigmamin))
        self.rs_sig.add_restraint(IMP.isd.JeffreysRestraint(self.m, sigma))
        
        return sigma
    
    def get_particles_to_sample(self):
        ps = {}
        for sigma_name in self.sigma_dictionary:
            ps["PEMAP_Sigma_" +
               str(sigma_name) + self._label_suffix] =\
                   ([self.sigma_dictionary[sigma_name][0]], self.sigma_dictionary[sigma_name][1])

        print('PEMAP RESTRAINT> Number of nuisance particles:', len(ps))

        return ps
        
    def get_output(self):
        #self.m.update()
        
        output = super(PEMAPRestraint, self).get_output()

        for sigma_name in self.sigma_dictionary:
            output["PEMAPRestraint_" +
                   str(sigma_name) + self.label] = str(
                       self.sigma_dictionary[sigma_name][0].get_scale())
    
        '''
        satif = 0
        tot = 0
        for i in range(len(self.pairs)):

            p0 =self.pairs[i][0]
            p1 =self.pairs[i][2]
            resid0 =self.pairs[i][1]
            resid1 =self.pairs[i][3]
            d = self.pairs[i][4]
            #up = self.pairs[i][5]
            
            #ln = self.pairs[i][6]

            label = 'p0'+":"+str(resid0)+"_"+'p1'+":"+str(resid1)+'_'+str(d)
            output["SimplifiedPEMAP_Score_"+label]=str(self.weight*ln.unprotected_evaluate(None))

            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            dist_sim = IMP.core.get_distance(d0,d1)
            output["SimplifiedPEMAP_Distance_"+label]=str(dist_sim)

            if dist_sim < 1.5*d:
                satif += 1.
            tot += 1
        output["SimplifiedPEMAP_Satisfied"]=str(satif/tot)
        ''' 
        return output

class COMDistanceRestraint(IMP.pmi.restraints.RestraintBase):
    
    def __init__(self,
                 root_hier,
                 protein0,
                 protein1,
                 distance = 60.0,
                 strength = 1.0,
                 label = None,
                 weight = 1.0):

        """ Setup an upper-bound distance restraint
            between the center-of-mass of two proteins
        @ param root_hier 
        @ param protein1
        @ param protein2
        @ param distance
        @ param strength
        @ param label
        @ param weight
        """

        self.root_hier = root_hier
        
        model = self.root_hier.get_model()
      
        rname = "COMDistanceRestraint"
        super(COMDistanceRestraint, self).__init__(
            model, name="COMDistanceRestraint", label=label, weight=weight)

    
        self.strength = strength
        self.weight = weight
        self.distance = distance
        
        # Setup restraint
        self.rs = self._create_restraint_set()
        
        
        s0 = IMP.atom.Selection(root_hier,
                                molecule = protein0)
        
        p0 = s0.get_selected_particles()
        if len(p0)==0:
            print("COMDistanceRestraint: WARNING> cannot select protein %s)" % (protein0))
            exit()
        s1 = IMP.atom.Selection(root_hier,
                                molecule = protein1) 
        p1 = s1.get_selected_particles()
        if len(p1)==0:
            print("COMDistanceRestraint: WARNING> cannot select protein %s)" % (protein1))
            exit()

        # Get COMs
        self.com0 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), p0)
        self.com1 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), p1)

        coor0 = IMP.core.XYZ(self.com0).get_coordinates()
        coor1 = IMP.core.XYZ(self.com1).get_coordinates()

        d = IMP.algebra.get_distance(coor0, coor1)
    
        # Distance restraint
        hub = IMP.core.HarmonicUpperBound(self.distance,
                                          self.strength)
        
        df = IMP.core.DistancePairScore(hub)
        dr = IMP.core.PairRestraint(model, df, (self.com0, self.com1))
        self.rs.add_restraint(dr)

    def get_output(self):

        output = super(COMDistanceRestraint, self).get_output()
        
        return output



