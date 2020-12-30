"""@namespace IMP.pmi.restraints.pemap
Restraints for handling pemap data.
"""

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


class PEMAPRestraint(IMP.pmi.restraints.RestraintBase):
    
    def __init__(self,
                 root_hier,
                 mic_file,
                 sigma_init = 4.0,
                 weight=1.0,
                 label = None):

        """
        Constructor:
        @param root_hier         Hierarchy of the system
        @param mic_file          File containing p1,r1,p2,r2,MIC
        @param sigma_init        Initial value for noise parameter
        @param weight            Weight of the restraint
        @param label             Extra text to label the restraint so that it is
                                 searchable in the output  
        """

        self.root_hier = root_hier
        self.m = self.root_hier.get_model()
        
        rname = "PEMAPRestraint"
        super(PEMAPRestraint, self).__init__(
            self.m, name="PEMAPRestraint", label=label, weight=weight)

        
        self.mic_file = mic_file
        self.sigma_init = sigma_init
        self.weight = weight
        
        # Nuisance particles
        self.sigma_dictionary={}
        
        self.sigma_is_sampled = True
        self.rs_sig = self._create_restraint_set("sigma")
        self._create_sigma('sigma_0', self.sigma_init)
        sigma_0 = self.sigma_dictionary['sigma_0'][0].get_particle_index()
    
        # Pairs of restrainted beads 
        self.pairs = []
       
        # Setup restraint
        self.pmr = IMP.isd.PEMAPRestraint(self.m,
                                          sigma_0)

        
        self.rs = IMP.RestraintSet(self.m, 'PEMAPRestraint')
        
        # Read file
        self._read_mic_file()

        print('PEMAP restraints> Number of restrained pairs:', len(self.pairs))
        
        # Add pairs to restraint
        for (p0, p1, d_mic) in self.pairs:
            self.pmr.add_contribution(p0.get_index(), p1.get_index(), d_mic)
            
        rname = self.pmr.get_name()
        self.rs.add_restraint(self.pmr)

        self.restraint_sets = [self.rs] + self.restraint_sets[1:]

    def _read_mic_file(self):
        # Read file
        for line in open(self.mic_file):
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
            s0 = IMP.atom.Selection(self.root_hier,
                                    molecule = prot1,
                                    residue_index = resi1).get_selected_particles()
            s1 = IMP.atom.Selection(self.root_hier,
                                    molecule = prot2, 
                                    residue_index = resi2).get_selected_particles()
            if len(s0)>0 and len(s1)>0:
                
                p0 = s0[0]
                p1 = s1[0]

                self.pairs.append((p0, p1, dist))

            else:
                 print("PEMAP restraint: WARNING> residue %d of chain %s is not there (w/ %d %s)"
                       % (resi1,prot1,resi2,prot2))
        
        
    def get_d_mic(self, mic):
        k = -0.014744
        n = -0.41
        if mic <0.6:
            return (math.log(mic)-n)/k
        else:
            return (math.log(0.6)-n)/k
        
    def get_upper_bound(self,mic,k,n):   
        if mic <= 0.6:
            return (log(mic) - n) / k
        else: return (log(0.6) - n) / k
    
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

        return output

class COMDistanceRestraint(IMP.pmi.restraints.RestraintBase):
    
    def __init__(self,
                 root_hier,
                 protein0,
                 protein1,
                 distance = 60.0,
                 strength = 1.0,
                 label = None):

        """ Setup an upper-bound distance restraint
            between the center-of-mass of two proteins
        @ param root_hier    Hierarchy of the system 
        @ param protein1     Name of first protein being restrained
        @ param protein2     Name of second protein being restrained
        @ param distance     Target distance between centers-of-mass
        @ param strength     Strenght of the harmonic restraint
        @ param label        Extra text to label the restraint so that it is
                             searchable in the output  
        @ param weight
        """

        self.root_hier = root_hier
        
        model = self.root_hier.get_model()
      
        rname = "COMDistanceRestraint"
        super(COMDistanceRestraint, self).__init__(
            model, name="COMDistanceRestraint", label=label, weight=weight)

    
        self.strength = strength
        self.distance = distance
        
        # Setup restraint
        self.rs = self._create_restraint_set()
                
        s0 = IMP.atom.Selection(self.root_hier,
                                molecule = protein0)
        
        p0 = s0.get_selected_particles()
        if len(p0)==0:
            raise TypeError("COMDistanceRestraint: WARNING> cannot select protein %s)" % (protein0))
        s1 = IMP.atom.Selection(self.root_hier,
                                molecule = protein1) 
        p1 = s1.get_selected_particles()
        if len(p1)==0:
            raise TypeError("COMDistanceRestraint: WARNING> cannot select protein %s)" % (protein1))
            
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


