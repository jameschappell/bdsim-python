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


'''Simple script used to generate simulations that use varying values of the quadrupole
strength for qd0.'''

'''Template for spectrometer.gmad'''

spectrometergmad='''! This file describes the basic spectrometer setup for the AWAKE experiment
! at CERN. 

option,physicsList="";
!option,ngenerate=100000;
include ../data/beam/test/testBeam_0_to_10GeV.gmad;
include options.gmad;

! Careful here
precReg: cutsregion, prodCutProtons=10*um,
                 prodCutPhotons=10*um,
                 prodCutElectrons=10*um,
                 prodCutPositrons=10*um;


beamEnergy=0.5*GeV;
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
		    !integrator="g4classicalrk4",
		    integrator="g4nystromrk4", 
		    magneticFile="bdsim3d:/unix/pdpwa/jchappell/data/fields/BDSIM_uniform_field__1.3598T__length_4.0m.dat",
		    magneticInterpolator="nearest3d";
		    !magneticFile="bdsim3d:/unix/pdpwa2/ldeacon/data/fieldMaps/hb4/3D/bdsim_0.9/HB4_MAP3D_540A.dat",
		    !magneticInterpolator="linear3d",
		    !x=0.13*m,
		    !z=-0.495*m;
		    

!dipole: rbend, l = 1*m, apertureType="rectangular", aper1=(0.4+offset_val)*m, aper2=35.8*mm, outerDiameter=1.2*m, fieldVacuum="dipolefield", offsetX=offset_val;
!dipole: rbend, l = 1*m, aper1=32*mm, outerDiameter=320*mm, fieldVacuum="dipolefield";! offsetX=offset_val;


!d5: drift, l = 1*mum, apertureType="rectangular", aper1=0.5*m, aper2=32*mm;

!d6: drift, l = 1*m, aper1=32*mm;

psz=1005*mm;
sez=1.676;
screenAngle=-45*pi/180;

mySpectrometer: awakespectrometer, l=4*m,
			                       region="precReg",
                                   twindow=1*cm,
                                   windowmaterial="G4_Al",
                                   tscint=850*mum,
                                   windowScreenGap=0*mm,
                                   scintmaterial="lanex",
                                   angle=screenAngle,poleStartZ=psz,
                                   screenEndZ=sez,
                                   spec="vacuumChamberType=1&magnetGeometryType=1",
                                   screenPSize=25*um,
                                   fieldAll="dipolefield";
				   !B=1.433;
				   
mSpectrIn : marker;
mSpectrOut: marker;

lSpectrometer: line = (mSpectrIn, mySpectrometer, mSpectrOut);

!specmeasureline: line=(d5);

awakespec: line=(d1, mprequad, qf0, d2, qd0, mpostquad, d3, d4, lSpectrometer);

!awakespec: line=(d1, mprequad, d1q, d2, d2q, mpostquad, d3, d4, lSpectrometer);	!modelling without quads

complete: line=(linePlasmaExit, awakespec);

!complete: line=(mySpectrometer);

use, period=complete;
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

def make_environment(s, resultsdir, beam_dir):

    '''Edits the spectrometergmad template to add the quadrupole correction
    factor and relevant beam file.'''

    print "Quadrupole correction factor is %f" %(s)
    sg = os.path.join(resultsdir, 'spectrometer.gmad')
    fh = open(sg, "wb")
    k1_string = 'k1=-' + str(s) + '*quadStrength;'
    spectrometergmad1 = string.replace(spectrometergmad, 'k1_corrected_def', k1_string)
    beamdir = beam_dir
    beam_input_string = 'include ../' + str(beamdir) + '/' + str(beamdir) + '.gmad'
    spectrometergmad2 = string.replace(spectrometergmad1, 'beam_include', beam_input_string)
    fh.write(spectrometergmad2)
    fh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    This script generates a number of BDSIM simulations each with a different 
    quadrupole correction factor for the second quadrupole.""",
    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--file', dest='file', default=None,
                        help = '''
    This is the path to the example beam .txt file that you are basing the new
    beam description on. It provides the position and momenta values for each
    particle.
    
    E.g. --file "<path to file>"''')

    parser.add_argument('--energy', dest='energy', default=None,
                        help='''
        This defines the mean energy of the gaussian beam that you would like to 
        generate. It is in units of GeV.

        E.g. --energy 5.0''')

    parser.add_argument('--spread_type', dest='type', default='percentage',
                        choices=['percentage', 'value'],
                        help='''
        This defines the type of energy spread that you are defining. The options
        are percentage or value. 

        percentage corresponds to:
        meanE +/- spread * meanE

        value corresponds to:
        meanE +/- spread

        value is defined in units of GeV.''')

    parser.add_argument('--spread', dest='spread', default=None,
                        help='''
        This defines the spread of the gaussian beam according to the type given
        in the --spread_type argument. percentage corresponds to a percentage of 
        the mean energy, whereas value corresponds to a particular energy spread
        value measured in GeV.''')

    parser.add_argument('--quad_corr_min', dest='quad_min', default=0.95,
                        help='''
        This is the minimum value of the quadrupole correction factor that you 
        want to scan over. Default value is 0.9''')

    parser.add_argument('--quad_corr_max', dest='quad_max', default=1.05,
                        help='''
            This is the maximum value of the quadrupole correction factor that you 
            want to scan over. Default value is 1.05''')

    parser.add_argument('--quad_corr_step', dest='quad_step', default=0.01,
                        help='''
            This is the step value you want to use to increment the quadrupole 
            correction factor between [quad_corr_min, quad_corr_max].
            Default value is 0.01.''')

    arguments = parser.parse_args()

    filename = arguments.file

    qd0_strengths = np.arange(float(arguments.quad_min), float(arguments.quad_max) + 0.0001,
                             float(arguments.quad_step))

    meanE = float(arguments.energy)

    spread = float(arguments.spread)

    cwd = os.getcwd()

    gen_beam_command = "python ${PYTHON_SCRIPTS_DIR}/generate_beam.py --file " + str(filename) + \
                       " --energy " + str(meanE) + \
                       " --spread_type " + str(arguments.type) + \
                       " --spread " + str(spread)

    os.system(gen_beam_command)

    if arguments.type == 'percentage':

        beam_dir = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread)

    if arguments.type == 'value':

        spread_fraction = spread / meanE

        beam_dir = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread_fraction)

    for strength in qd0_strengths:

        """Looping over the range of quadrupole correction factors, creating a new 
        directory for each simulation."""

        os.chdir(cwd)
        print "Quadrupole Correction Factor: %f" %(strength)
        res_dir = 'qd0_correction_' + str(strength)
        print "Making directory: ", res_dir
        if os.path.isdir(res_dir) is False:
            os.mkdir(res_dir)
        else:
            print "Directory %s already exists. Stopping." %(res_dir)
            break
        print "Copy template files: "
        for filename in glob.glob(os.path.join(TEMPLATES_DIR, '*')):
            shutil.copy(filename, res_dir)
        make_environment(strength, res_dir, beam_dir)
        os.chdir(res_dir)
        run_command = "bdsim --file=spectrometer.gmad --outfile=qd0_correction__" + str(strength) + "_" + " --batch"
        print run_command
        os.system(run_command)
        outfile = "qd0_correction__" + str(strength) + "__event.root"
        shutil.copy(outfile, cwd)
