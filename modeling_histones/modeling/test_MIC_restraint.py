#import utils
import os
import unittest

import sys
sys.path.append(r'/Users/iecheverria/SOFTW/python-ihm/ihm/')

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from io import BytesIO as StringIO

import ihm.dumper
import ihm.format
import ihm.location
import ihm.representation
import ihm.startmodel
import ihm.dataset
import ihm.protocol
import ihm.analysis
import ihm.model
import ihm.restraint
import ihm.geometry



def _get_dumper_output(dumper, system):
    fh = StringIO()
    writer = ihm.format.CifWriter(fh)
    dumper.dump(system, writer)
    return fh.getvalue()


"""Test DerivedDistanceRestraintDumper"""
class MockObject(object):
    pass

system = ihm.System()

feat1 = MockObject()
feat1._id = 44
feat2 = MockObject()
feat2._id = 84
dataset = MockObject()
dataset._id = 97


dist = ihm.restraint.LowerBoundDistanceRestraint(25.0)
r1 = ihm.restraint.DerivedDistanceRestraint(dataset=dataset,
                                            feature1=feat1, feature2=feat2, distance=dist,
                                            probability=0.8)
r2 = ihm.restraint.DerivedDistanceRestraint(dataset=dataset,
                                            feature1=feat1, feature2=feat2, distance=dist,
                                            probability=0.4)
r3 = ihm.restraint.DerivedDistanceRestraint(dataset=dataset,
                                            feature1=feat1, feature2=feat2, distance=dist,
                                            probability=0.6)
rg = ihm.restraint.RestraintGroup((r2, r3))
system.restraints.extend((r1, r2)) # r2 is in restraints and groups
system.restraint_groups.append(rg)

dumper = ihm.dumper._DerivedDistanceRestraintDumper()
dumper.finalize(system) # assign IDs
