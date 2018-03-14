import sys
import numpy as np
import os
import glob
import shutil
import string
import generate_beam
import stat
import subprocess
import argparse

TEMPLATES_DIR = '/unix/pdpwa/jchappell/SIM_TEMPLATES'

spectrometergmad = '''! This file describes the basic spectrometer setup for the AWAKE experiment
! at CERN. 

option,physicsList="";
!option,ngenerate=100000;
include ../../data/beam/test/testBeam_0_to_10GeV.gmad;
include options.gmad;

! Careful here
precReg: cutsregion, prodCutProtons=10*um,
                 prodCutPhotons=10*um,
                 prodCutElectrons=10*um,
                 prodCutPositrons=10*um;


beamEnergy=1.3*GeV;
beam, particle="e-",
      energy=beamEnergy;
!      distrType="gauss",
!      sigmaX=250e-6,
!      sigmaY=250e-6,
!      sigmaE=0.3*GeV,
!      sigmaT=4e-12;

! Quad Geometry
quadStrength=4.678362;
quad_physical_length=285*mm;
quad_effective_length=310*mm;
quad_length_diff=quad_effective_length-quad_physical_length;
aper_quad=70*mm;

! beam line elements
def_aper=35*mm;

! Final iris
!iris1: ecol, xsize=20*mm, ysize=20*mm, material="G4_cu", l=0.6*mm, outerDiameter=35*mm;
iris_plasma_exit: ecol, xsize=20*mm, ysize=20*mm, l=0.6*mm, outerDiameter=50*mm;
m_iris_plasma: marker;
linePlasmaExit: line=(iris_plasma_exit, m_iris_plasma);

! Drift from exit of plasma cell
d1: drift,l=5.753*m, aper1=def_aper;

! Quadrupoles
mprequad:  marker;
d1q: drift, l=285*mm, aper1=def_aper; 	!for modelling without quads
qf0: quadrupole, l=285*mm, aper1=def_aper, 
k1=quadStrength;
d2: drift, l=210*mm, aper1=def_aper;
d2q: drift, l=285*mm, aper1=def_aper;	!for modelling without quads
qd0: quadrupole, l=285*mm, aper1=def_aper,
k1=-1*quadStrength;
mpostquad: marker;
d3: drift, l=1*mum, aper1=def_aper;

! Drift to spec
d4: drift, l=659.6*mm, aper1=def_aper;

! Simple Spectrometer model
! Dipole
!psz=1005*mm;
!sez=257.661*cm-psz;
!screenAngle=45*pi/180;
offset_val=0.125*m;

! Field
dipolefield: field, type="bmap3d", 
		    integrator="g4classicalrk4", 
		    magneticFile="bdsim3d:/unix/pdpwa2/ldeacon/data/fieldMaps/hb4/3D/bdsim_0.9/HB4_MAP3D_650A.tar.gz",
		    magneticInterpolator="nearest3d";
		    !x=-1*offset_val;

!dipole: rbend, l = 1*m, apertureType="rectangular", aper1=(0.4+offset_val)*m, aper2=35.8*mm, outerDiameter=1.2*m, fieldVacuum="dipolefield", offsetX=offset_val;
!dipole: rbend, l = 1*m, aper1=32*mm, outerDiameter=320*mm, fieldVacuum="dipolefield";! offsetX=offset_val;


!d5: drift, l = 1*mum, apertureType="rectangular", aper1=0.5*m, aper2=32*mm;

!d6: drift, l = 1*m, aper1=32*mm;

psz=1005*mm;
sez=257.661*cm-psz;
screenAngle=-45*pi/180;

mySpectrometer: awakespectrometer, l=4.4284*m-quad_length_diff/2,
                                   region="precReg",
                                   twindow=1*cm,
                                   windowmaterial="G4_Al",
                                   tscint=850*mum,
                                   windowScreenGap=65*mm,
                                   scintmaterial="lanex",
                                   angle=screenAngle,poleStartZ=psz,
                                   screenEndZ=sez,
                                   spec="vacuumChamberType=1&magnetGeometryType=1",
                                   screenPSize=25*um,
                                   B_string_def;

mSpectrIn : marker;
mSpectrOut: marker;

lSpectrometer: line = (mSpectrIn, mySpectrometer, mSpectrOut);

!specmeasureline: line=(d5);

awakespec: line=(d1, mprequad, qf0, d2, qd0, mpostquad, d3, d4, lSpectrometer);

!awakespec: line=(d1, mprequad, d1q, d2, d2q, mpostquad, d3, d4, lSpectrometer);	!modelling without quads

complete: line=(linePlasmaExit, awakespec);

!complete: line=(mySpectrometer);

use, period=complete;

!specmeasurelinePlace: placement, sequence="specmeasureline",
!				 referenceElement="d6",
!				 referenceElementNumber=0,
!				 x = (352.5+50)*mm,
!				 z = -1*m,
!				 !psi = 0.8;
!				 axisAngle = 1,
!				 axisY = 1,
!				 angle = screenAngle;
!
!sample, range = m_iris_plasma;
!sample, range = mprequad;
!sample, range = mpostquad;
!sample, all;


! Beam
!beam,   particle="e-",
!        energy=1.3*GeV,
!        sigmaE=0.4*GeV;

'''


def make_environment(current, resultsdir):

    print "Current is %s A." % (current)
    sg = os.path.join(resultsdir, 'spectrometer.gmad')
    fh = open(sg, "wb")
    b_string = 'B=0.5*(0.00154/(650/' + str(current) + '));'
    spectrometergmad1 = string.replace(spectrometergmad, 'B_string_def', b_string)
    fh.write(spectrometergmad1)
    fh.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        This script generates a number of BDSIM simulations each with a different 
        current value for the dipole within the spectrometer.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--I_min', dest='I_min', default=None,
                        help='''
        This is the minimum value of the current supplied to the dipole magnet
        within the spectrometer. Measured in Amperes.

        E.g. --I_min 40''')

    parser.add_argument('--I_max', dest='I_max', default=650,
                        help='''
            This is the maximum value of the current supplied to the dipole magnet
            within the spectrometer. Measured in Amperes. Default value is 650A.

            E.g. --I_min 650''')

    parser.add_argument('--I_step', dest='I_step', default=50,
                        help='''
            This is the step value you want to use to increment the current 
            for the spectrometer dipole, ranging from [I_min, I_max].
            Default value is 50A.''')

    arguments = parser.parse_args()

    currents = np.arange(float(arguments.I_min), float(arguments.I_max) + 0.0001,
                                float(arguments.I_step))

    cwd = os.getcwd()

    for current in currents:

        """Looping over the range of current values, creating a new 
        directory for each simulation."""

        os.chdir(cwd)
        print "Current: " + str(current) +  "A"
        res_dir = 'current_' + str(current) + 'A'
        print "Making directory: ", res_dir
        if os.path.isdir(res_dir) is False:
            os.mkdir(res_dir)
        else:
            print "Directory %s already exists. Stopping." % (res_dir)
            break
        print "Copy template files: "
        for filename in glob.glob(os.path.join(TEMPLATES_DIR, '*')):
            shutil.copy(filename, res_dir)
        make_environment(current, res_dir)
        os.chdir(res_dir)
        run_command = "bdsim --file=spectrometer.gmad --outfile=current__" \
                      + str(current) + "_" + " --batch"
        print run_command
        os.system(run_command)
        outfile = "current__" + str(current) + "__event.root"
        shutil.copy(outfile, cwd)


