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
import uniformfield

TEMPLATES_DIR = '/unix/pdpwa/jchappell/SIM_TEMPLATES'


'''Simple script used to generate simulations that use varying values of the quadrupole
strength for qd0.'''

'''Template for spectrometer.gmad'''

spectrometergmad='''! This file describes the basic spectrometer setup for the AWAKE experiment
! at CERN. 

option,physicsList="";
option,ngenerate=100000;
beam_include
include options.gmad;

! Careful here
precReg: cutsregion, prodCutProtons=10*um,
                 prodCutPhotons=10*um,
                 prodCutElectrons=10*um,
                 prodCutPositrons=10*um;


beamEnergydef;
beam, particle="e-",
      energy=beamEnergy,
      distrType="gauss",      
      sigmaXdef,
      sigmaYdef,
      sigmaEdef,
      sigmaT=4e-12,
      sigmaXpdef,
      sigmaYpdef;

! Quad Geometry
quadStrength=5.8;
quad_physical_length=285*mm;
quad_effective_length=310*mm;
quad_length_diff=quad_effective_length-quad_physical_length;
aper_quad=70*mm;

! beam line elements
def_aper=35*mm;

!Plasma cell
!dPlasma  : drift, aper=def_aper, l=10*m;
irisAper = 5*mm;
liris    = 0.6*mm;
iris1    : ecol, xsize=irisAper, ysize=irisAper, material="G4_Al",
           l=liris, outerDiameter=2*def_aper;
iris2    : iris1;
mIris1In : marker;
mIris1Out: marker;
mIris2In : marker;
mIris2Out: marker;

!Final plasma cell iris:
lineIris2: line=(mIris2In,iris2,mIris2Out);

lineIris2_l=liris;

d1: drift, l=850*mm, aper=def_aper;

!Laser dump
l_ldump1        = 246*mm;
w_ldump1        = 1*cm;
diagScreenAngle = 45*pi/180.0;
ldump1_s        : screen, l=l_ldump1, layerThicknesses={ 0.2*mm }, 
                  layerMaterials=["G4_Al"], screenXSize=w_ldump1, 
                  screenYSize=w_ldump1, angle=diagScreenAngle, 
                  aper=def_aper;
ldump1_d        : drift, l=l_ldump1, aper=def_aper;

d2: drift, aper=def_aper, l=469*mm;

bpm1: drift, aper=def_aper, l=248*mm;

d3: drift, aper=def_aper, l=84*mm;

!BTV1 (BI screen)
l_btv1 = 350*mm;
w_btv1 = w_ldump1;
btv1_s : screen, l=l_btv1, layerThicknesses={ 0.3*mm, 300*nm }, 
         layerMaterials=["G4_Si","G4_Ag"],screenXSize=w_btv1, 
         screenYSize=w_btv1, angle=diagScreenAngle, aper=def_aper;	     
btv1_d : drift, l=l_btv1, aper=def_aper;

dtab1_l = 280*mm;
dtab1   : drift, aper=def_aper, l=dtab1_l;

!OTR screen
l_otr  = 200*mm;
w_otr  = w_ldump1;
otr_s  : screen, l=l_otr, layerThicknesses={ 0.15*mm, 100*nm }, 
         layerMaterials=["G4_Si","G4_Ag"],screenXSize=w_otr, 
         screenYSize=w_otr, angle=diagScreenAngle, aper=def_aper;
otr_d  : drift, l=l_otr, aper=def_aper;

dtab2_l = 280*mm;
dtab2   : drift, aper=def_aper, l=dtab2_l;

!CTR screen
l_ctr = l_otr;
w_ctr = w_ldump1;
d4    : drift,l=137*mm, aper=def_aper;
ctr_s : screen, l=l_ctr, layerThicknesses={ 0.15*mm, 100*nm }, 
        layerMaterials=["G4_Si","G4_Ag"],screenXSize=w_ctr, 
        screenYSize=w_ctr, angle=diagScreenAngle, aper=def_aper;
ctr_d : drift, l=l_ctr, aper=def_aper;

dtab3_l = 280*mm;
dtab3   : drift, aper=def_aper,l=dtab3_l;

ctab3_l = l_otr;
ctab3   : drift, aper=def_aper,l=ctab3_l;

l_diagnosticTable=2*m;

dtab4: drift, aper=def_aper, l=l_diagnosticTable-dtab1_l - l_otr - 
       dtab2_l - l_ctr - dtab3_l - ctab3_l;

min        : marker;
mout       : marker;
motr_in    : marker;
motr_out   : marker;
mctr_in    : marker;
mctr_out   : marker;
mctab3_in  : marker;
mctab3_out : marker;


diagnosticTable: line=(dtab1, motr_in, otr_d, motr_out, dtab2, 
                 mctr_in, ctr_d, mctr_out, dtab3, mctab3_in, 
                 ctab3, mctab3_out, dtab4);

lined5l   = 345*mm-quad_length_diff/2;
d5shieldl = lined5l*0.99;
d5        : drift, aper=def_aper, l=345*mm;!lined5l-d5shieldl;
d5ColAper = 10*mm;
d5col     : ecol, l=d5shieldl, xsize=d5ColAper, ysize=d5ColAper, 
            outerDiameter=0.1*m, material="G4_W";

lined5: line=(d5,d5col);

! Quadrupoles
mprequad:  marker;
d1q: drift, l=285*mm, aper1=def_aper; 	!for modelling without quads
qf0: quadrupole, l=285*mm, aper1=def_aper, 
k1_x_corrected_def
d6:  drift,l=210*mm, aper=def_aper;
d2q: drift, l=285*mm, aper1=def_aper;	!for modelling without quads
qd0: quadrupole, l=285*mm, aper1=def_aper,
k1_y_corrected_def
mpostquad: marker;
d7:  drift,l=1005*mm, aper=def_aper;

offset_val=0.125*m;

! Field
dipolefield: field, type="bmap3d", 
		    !integrator="g4classicalrk4",
		    integrator="g4nystromrk4", 
		    fieldlocdef
		    fieldinterpolatordef
		    xoffsetdef
		    zoffsetdef
		    

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
                                   windowScreenGap=0.1*mm,
                                   scintmaterial="lanex",
                                   angle=screenAngle,poleStartZ=psz,
                                   screenEndZ=sez,
                                   spec="vacuumChamberType=1&magnetGeometryType=1",
                                   screenPSize=25*um,
                                   fieldAll="dipolefield";
				   
mSpectrIn : marker;
mSpectrOut: marker;

lSpectrometer: line = (mSpectrIn, mySpectrometer, mSpectrOut);

d8 : drift,l=2423*mm, aper=def_aper;
     !,tunnelOffsetX=tunOffX;

!BTV2 (BI screen)
diagScreenAngle = 45*pi/180.0;
l_btv2          = 448*mm;
w_btv2          = w_ldump1;
btv2_s          : screen, l=l_btv2, layerThicknesses={ 0.3*mm, 300*nm }, 
                  layerMaterials=["G4_Si","G4_Ag"],screenXSize=w_btv2, 
                  screenYSize=w_btv2, angle=diagScreenAngle, aper=def_aper;	     
btv2_d          : drift, l=l_btv2, aper=def_aper;

d9_l = 1711*mm;!123*mm;
d9   : drift, aper=def_aper, l=d9_l;

!Laser dump
l_ldump2   = l_ldump1;
w_ldump2   = w_ldump1;
ldump2_s_l = l_ldump2;
ldump2_s   : screen, l=ldump2_s_l, layerThicknesses={ 0.2*mm }, 
             layerMaterials=["G4_Al"], screenXSize=w_ldump2, 
             screenYSize=w_ldump2, angle=diagScreenAngle, aper=def_aper;
ldump2_d   : drift, l=ldump2_s_l, aper=def_aper;

d10: drift, aper=def_aper, l=8006*mm;!1670*mm-ldump2_s_l-d9_l;

!downstream of the screen (tunnel dimensions etc.)
tunRad_2     = 3.3*1000*m;
tunOffX_2    = 1.4*m;
beamCeilH_2  = 3.8*1000*m;
floorBeamH_2 = 1.6*1000*m;
thid9a       = tunRad_2*1000*m-tunRad_2+tunOffX_2*1000-
               tunOffX_2*1000;!+tunThi*1000;

!d9a                 : drift,l=ld9a, aper=def_aper;
!,tunnelOffsetX      = tunOffX_2, 
!tunnelRadius        = tunRad_2, tunnelThickness = thid9a, 
!                      beamlineCeilingHeight=beamCeilH_2, 
!floorBeamlineHeight = floorBeamH_2, tunnelType=0;

!d9b                 : drift,l=ld9b, aper=def_aper;
!,tunnelOffsetX      = tunOffX_2, 
!tunnelRadius        = tunRad_2, beamlineCeilingHeight=beamCeilH_2, 
!floorBeamlineHeight = floorBeamH_2, tunnelType=0;

!d9: line = (d9a, d9b);

!tunRad_3  = 1e3*(bpRad+bpThi+1*mm);
!doorWidth = 0*1e3*1*m;
!tunThi_3  = tunRad_2-tunRad_3+2*tunOffX_2-doorWidth;
!d10       : drift,l=0.5*m, aper=bpRad, tunnelRadius=tunRad_3, 
!            beamlineCeilingHeight=tunRad_3, 
!            floorBeamlineHeight=tunRad_3, tunnelType=1, 
!            tunnelThickness=tunThi_3, tunnelSoilThickness=1*mm;
!            !tunnelOffsetX = 0; 


doorWidth = 2*m;
endwall   : ecol, l=80*cm, xsize=def_aper, ysize=def_aper, 
            outerDiameter=2*1e-3*(tunRad_2+tunOffX_2-doorWidth), 
            material="concrete";
            !tunnelRadius=1000*10*m, beamlineCeilingHeight=1000*10*m, 
            !floorBeamlineHeight=1000*10*m, buildTunnel=0;
            !tunnelOffsetX=0, 

mbtv1_in    : marker;
mbtv1_out   : marker;
mldump_in   : marker;
mldump_out  : marker;
msptr_in    : marker;
msptr_out   : marker;
mbtv2_in    : marker;
mbtv2_out   : marker;
mldump2_in  : marker;
mldump2_out : marker;
mendwall    : marker;

mqf0: marker;

upstream : line=(lineIris2, d1, mldump_in, mldump_out,d2, 
           bpm1, d3, mbtv1_in, btv1_d, mbtv1_out, d4, 
           diagnosticTable, d5, mqf0, qf0, d6, qd0,d7,lSpectrometer);

downstream : line=(d8, mbtv2_in, btv2_d, mbtv2_out, d9, mldump2_in,mldump2_out);
             !ldump2_d, mldump2_out);

endwallLine: line=(d10, mendwall, endwall);

awakesbl: line=(min, upstream, downstream, mout);

use, period=awakesbl;
!use, period=lineIris2;

sample, range=mIris2In;
sample, range=mIris2Out;

sample, range = min;
sample, range = mout;
sample, range = motr_in;
sample, range = motr_out;
sample, range = mctr_in;
sample, range = mctr_out;
sample, range = mctab3_in;
sample, range = mctab3_out;

sample, range = mqf0;

sample, range = mSpectrIn ;
sample, range = mSpectrOut;

sample, range = mbtv1_in;
sample, range = mbtv1_out;
sample, range = mldump_in;
sample, range = mldump_out;

sample, range = mbtv2_in;
sample, range = mbtv2_out;
sample, range = mldump2_in;
sample, range = mldump2_out;
!sample, range = mendwall;
'''

subscript = '''
#!/bin/bash

#PBS -l walltime=02:00:00
#PBS -l mem=1G
#PBS -l nodes=1
#PBS -q medium
name
#PBS -m be
#PBS -M james.chappell.17@ucl.ac.uk
work_dir
#PBS -e ./logs/
#PBS -o ./logs/

scl enable devtoolset-3 bash <<EOF
source /unix/pbt/software/scripts/dev/bds_pbt-dev.bashrc
source /unix/pdpwa/jchappell/setup.sh
change_dir

sub_line
EOF
'''

def make_environment(meanE, spread, emitt, x_strength, y_strength, resultsdir,
                     fieldtype,
                     fieldstrength):

    '''Edits the spectrometergmad template to add the quadrupole correction
    factor and relevant beam file.'''

    emitt_factor = np.sqrt(emitt/100)
    gamma = meanE/(0.511e-3)
    orig_gamma = 1.3/(0.511e-3)
    gamma_factor = np.sqrt(orig_gamma/gamma)
    baselinex = 0.2
    baselinexp = 2

    quad_array = [format(x_strength, '.2f'), format(y_strength, '.2f')]
    print "Quadrupole correction factor is %s" %(quad_array)
    sg = os.path.join(resultsdir, 'spectrometer.gmad')
    fh = open(sg, "wb")

    beamenergy_string = 'beamEnergy=' + str(meanE) + '*GeV;'
    spectrometergmade = string.replace(spectrometergmad, 'beamEnergydef',
                                       beamenergy_string)

    beamespread_string = 'sigmaE=' + str(spread)
    spectrometergmade1 = string.replace(spectrometergmade, 'sigmaEdef',
                                        beamespread_string)

    beamxspread_string = 'sigmaX=' + str(emitt_factor * gamma_factor *
                                         baselinex)
    spectrometergmade2 = string.replace(spectrometergmade1, 'sigmaXdef',
                                        beamxspread_string)

    beamyspread_string = 'sigmaY=' + str(emitt_factor * gamma_factor *
                                         baselinex)
    spectrometergmade3 = string.replace(spectrometergmade2, 'sigmaYdef',
                                        beamyspread_string)

    beamxpspread_string = 'sigmaXp=' + str(emitt_factor * gamma_factor *
                                         baselinexp)
    spectrometergmade4 = string.replace(spectrometergmade3, 'sigmaXpdef',
                                        beamxpspread_string)

    beamypspread_string = 'sigmaYp=' + str(emitt_factor * gamma_factor *
                                         baselinexp)
    spectrometergmade5 = string.replace(spectrometergmade4, 'sigmaYpdef',
                                        beamypspread_string)

    k1_y_string = 'k1=-' + str(y_strength) + '*quadStrength;'
    spectrometergmad1 = string.replace(spectrometergmade5,
                                       'k1_y_corrected_def', k1_y_string)



    k1_x_string = 'k1=' + str(x_strength) + '*quadStrength;'
    spectrometergmad3 = string.replace(spectrometergmad1,
                                       'k1_x_corrected_def', k1_x_string)

    # Change field definition depending on field type

    if fieldtype == "uniform":

        fieldmap_string = 'magneticFile="bdsim3d:/unix/pdpwa/jchappell/data' \
                         '/fields/BDSIM_uniform_field__' \
                         + str(fieldstrength) + 'T__length_4.0m.dat",'
        spectrometergmad4 = string.replace(spectrometergmad3, 'fieldlocdef',
                                           fieldmap_string)

        fieldinterpolator_string = 'magneticInterpolator="nearest3d";'
        spectrometergmad5 = string.replace(spectrometergmad4,
                                           'fieldinterpolatordef',
                                           fieldinterpolator_string)

        fieldoffsetx_string = ''
        spectrometergmad6 = string.replace(spectrometergmad5, 'xoffsetdef',
                                           fieldoffsetx_string)

        fieldoffsetz_string = ''
        spectrometergmad7 = string.replace(spectrometergmad6, 'zoffsetdef',
                                           fieldoffsetz_string)

    else:

        fieldmap_string = 'magneticFile="bdsim3d:/unix/pdpwa2/ldeacon/data' \
                          '/fieldMaps/hb4/3D/bdsim_0.9/HB4_MAP3D_' \
                          + str(fieldstrength) + 'A.dat",'
        spectrometergmad4 = string.replace(spectrometergmad3, 'fieldlocdef',
                                           fieldmap_string)

        fieldinterpolator_string = 'magneticInterpolator="linear3d",'
        spectrometergmad5 = string.replace(spectrometergmad4,
                                           'fieldinterpolatordef',
                                           fieldinterpolator_string)

        fieldoffsetx_string = 'x=0.13*m,'
        spectrometergmad6 = string.replace(spectrometergmad5, 'xoffsetdef',
                                           fieldoffsetx_string)

        fieldoffsetz_string = 'z=-0.495*m;'
        spectrometergmad7 = string.replace(spectrometergmad6, 'zoffsetdef',
                                           fieldoffsetz_string)

    fh.write(spectrometergmad7)
    fh.close()

    ss = os.path.join(resultsdir, 'sub_script.bash')
    fs = open(ss, "wb")

    name_string = '#PBS -N BDSIMquadscan'
    subscript1 = string.replace(subscript, 'name', name_string)

    work_dir_string = '#PBS -wd ' + os.getcwd() +'/' + res_dir
    subscript2 = string.replace(subscript1, 'work_dir', work_dir_string)

    change_dir_string = 'cd ' + os.getcwd() + '/' + res_dir
    subscript3 = string.replace(subscript2, 'change_dir', change_dir_string)

    sub_line_str = "bdsim --file=spectrometer.gmad --outfile=x_" \
                   + format(x_strength, '.2f') + "_y_" + \
                   format(y_strength, '.2f') + " --batch"
    subscript4 = string.replace(subscript3, 'sub_line', sub_line_str)

    fs.write(subscript4)
    fs.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    This script generates a number of BDSIM simulations each with a different 
    quadrupole correction factor for the second quadrupole.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--energy', dest='energy', default=None,
                        help='''
        This defines the mean energy of the gaussian beam that you would like to 
        generate. It is in units of GeV.

        E.g. --energy 5.0''')

    parser.add_argument('--spread', dest='spread', default=None,
                        help='''
        This defines the spread of the gaussian beam as a percentage of the 
        mean energy. e.g. --spread 0.1 is a 10% energy spread around the mean.
        .''')

    parser.add_argument('--emitt', dest='emitt', default=None,
                        help='''
            This defines the normalised emittance of the beam. It is measured 
            in mm*mrad.
            .''')

    parser.add_argument('--scan_type', dest='scan_type', default='quick',
                        choices=['quick', 'fine'], help='''
        This option allows you to choose between a quick scan or a fine scan. 
        
        A quick scan means that the quadrupole coefficients will be the same for
        both quadrupoles, varying from 0.0 to 1.0. 
        e.g. (0.0, 0.0), (0.1, 0.1), (0.2, 0.2) etc.
        
        A fine scan means you can define the quadrupole coefficients separately 
        through the later options.
        e.g. for y varying between 0.1 and 0.2 in steps of 0.01
             and x varying between 0.5 and 0.7 in steps of 0.02
        You would give the following arguments:
        
        --y_quad_corr_min 0.1 --y_quad_corr_max 0.2 --y_quad_corr_step 0.01
        --x_quad_corr_min 0.5 --x_quad_corr_max 0.7 --x_quad_corr_step 0.02
        
        The quadrupole coefficients would be:
        (0.5, 0.1), (0.5, 0.11), (0.5, 0.12), ..., (0.5, 0.2)
        (0.52, 0.1), (0.52, 0,11), (0.52, 0.12), ..., (0.52, 0.2)
        and so on...
        
        ''')

    parser.add_argument('--y_quad_corr_min', dest='y_quad_min', default=None,
                        help='''
        This is the minimum value of the quadrupole correction factor for the 
        focusing in the y-direction that you want to scan over. Default value 
        is 0.9''')

    parser.add_argument('--y_quad_corr_max', dest='y_quad_max', default=None,
                        help='''
        This is the maximum value of the quadrupole correction factor for the
        focusing in the y-direction that you want to scan over. Default value 
        is 1.05''')

    parser.add_argument('--y_quad_corr_step', dest='y_quad_step', default=None,
                        help='''
        This is the step value you want to use to increment the quadrupole 
        correction factor for the focusing in the y-direction between 
        [quad_corr_min, quad_corr_max]. Default value is 0.01.''')

    parser.add_argument('--x_quad_corr_min', dest='x_quad_min', default=None,
                        help='''
        This is the minimum value of the quadrupole correction factor for the 
        focusing in the x-direction that you want to scan over. Default value 
        is 0.9''')

    parser.add_argument('--x_quad_corr_max', dest='x_quad_max', default=None,
                        help='''
        This is the maximum value of the quadrupole correction factor for the 
        focusing in the y-direction that you want to scan over. Default value 
        is 1.05''')

    parser.add_argument('--x_quad_corr_step', dest='x_quad_step', default=None,
                        help='''
        This is the step value you want to use to increment the quadrupole 
        correction factor for the focusing in the x-direction between 
        [quad_corr_min, quad_corr_max]. Default value is 0.01.''')

    parser.add_argument('--field_type', dest='field_type', default='uniform',
                        choices=['uniform', 'custom'], help='''
        This defines whether you want to use a uniform field or a custom 
        field map to model the magnetic field of the spectrometer dipole. 
        Default value is uniform.''')

    parser.add_argument('--field_strength', dest='field_strength', help='''
        If you are using a uniform field map, you can use any value of the 
        field strength. However, if you want to use a custom field map 
        please enter the current value you want. You can choose between:
        
        Field Strength /T           Current /A
            0.1265                      40
            0.3176                      100
            0.5403                      170
            0.7615                      240
            1.0043                      320
            1.1937                      400
            1.3598                      540
            1.4333                      650
            
        The field strength is measured in Tesla.''')

    arguments = parser.parse_args()

    meanE = float(arguments.energy)

    spread = float(arguments.spread)

    emitt = float(arguments.emitt)

    cwd = os.getcwd()

    # Generate a beam file for use in the simulation.

    #gen_beam_command = "python ${PYTHON_GIT_DIR}/generate_beam.py --file " +
    #  str(filename) + \
    #                   " --energy " + str(meanE) + \
    #                   " --spread_type " + str(arguments.type) + \
    #                   " --spread " + str(spread)
    #print gen_beam_command

    #os.system(gen_beam_command)

    # If uniform field required, generate a field map for use in the simulation.

    if arguments.field_type == 'uniform':

        filename = '/unix/pdpwa/jchappell/data/fields/BDSIM_uniform_field__' \
                   + str(arguments.field_strength) + 'T__length_4.0m.dat'

        if os.path.isfile(filename):

            print "Uniform field map already exists. No need to generate a " \
                  "new one."

        else:

            gen_field_command = "python ${" \
                                "PYTHON_GIT_DIR}/uniformfield.py" +  " " \
                                "--strength " + str(arguments.field_strength)\
                                + " --length 4.0"

            print "Uniform field map required. Generating field with strength " +\
                  str(arguments.field_strength) + "T."

            #os.system(gen_beam_command)

    # Generate beam directory name depending on spread type.

    #if arguments.type == 'percentage':

    #    beam_dir = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' +  \
    #               str(spread)

    #elif arguments.type == 'value':

    #    spread_fraction = spread / meanE

    #    beam_dir = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + \
    #               str(spread_fraction)

    #else:

    #   beam_dir = 'e_beam_max_' + str(meanE) + '_stepwise'

    # Loop over different quad strengths.

    if arguments.scan_type == 'quick':

        qd0_strengths = np.linspace(0, 1.0, 11)

        for strength in qd0_strengths:

            x_strength = strength

            y_strength = strength

            """Looping over the range of quadrupole correction factors, creating
            a new directory for each simulation."""

            os.chdir(cwd)
            quad_array = [format(x_strength, '.2f'), format(y_strength, '.2f')]
            print "Quadrupole Correction Factor: %s" % (quad_array)
            res_dir = 'x_' + format(x_strength, '.2f') + '_y_' + \
                      format(y_strength, '.2f')
            print "Making directory: ", res_dir
            if os.path.isdir(res_dir) is False:
                os.mkdir(res_dir)
            else:
                print "Directory %s already exists. Stopping." % (res_dir)
                break
            print "Copy template files: "
            for filename in glob.glob(os.path.join(TEMPLATES_DIR, '*')):
                shutil.copy(filename, res_dir)
            field_type = str(arguments.field_type)
            field_strength = arguments.field_strength
            make_environment(meanE, spread, emitt, x_strength, y_strength,
                             res_dir, field_type, field_strength)
            os.chdir(res_dir)
            os.mkdir('logs')
            run_command = "qsub sub_script.bash"
            print run_command
            os.system(run_command)
            # outfile = "x_" + str(x_strength) + "_y_" + str(y_strength)\
            #          + ".root"
            # shutil.copy(outfile, cwd)


    else:

        y_qd0_strengths = np.arange(float(arguments.y_quad_min),
                                    float(arguments.y_quad_max) + 0.0001,
                                    float(arguments.y_quad_step))

        x_qf0_strengths = np.arange(float(arguments.x_quad_min),
                                    float(arguments.x_quad_max) + 0.0001,
                                    float(arguments.x_quad_step))

        for y_strength in y_qd0_strengths:

            for x_strength in x_qf0_strengths:

                """Looping over the range of quadrupole correction factors, creating a new 
                directory for each simulation."""

                os.chdir(cwd)
                quad_array = [format(x_strength, '.2f'), format(y_strength, '.2f')]
                print "Quadrupole Correction Factor: %s" %(quad_array)
                res_dir = 'x_' + format(x_strength, '.2f') + '_y_' + \
                          format(y_strength, '.2f')
                print "Making directory: ", res_dir
                if os.path.isdir(res_dir) is False:
                    os.mkdir(res_dir)
                else:
                    print "Directory %s already exists. Stopping." %(res_dir)
                    break
                print "Copy template files: "
                for filename in glob.glob(os.path.join(TEMPLATES_DIR, '*')):
                    shutil.copy(filename, res_dir)
                field_type = str(arguments.field_type)
                field_strength = arguments.field_strength
                make_environment(meanE, spread, emitt, x_strength,
                                 y_strength,  res_dir, field_type,
                                 field_strength)
                os.chdir(res_dir)
                os.mkdir('logs')
                run_command = "qsub sub_script.bash"
                print run_command
                os.system(run_command)
                #outfile = "x_" + str(x_strength) + "_y_" + str(y_strength)\
                #          + ".root"
                #shutil.copy(outfile, cwd)
