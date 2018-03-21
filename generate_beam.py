import produce_beam
import os
import argparse
import shutil
import string

e_beamgmad = """
beam,  particle="e-",
       distrType="userfile",
       distrFileFormat="x[mm]:xp[mrad]:y[mm]:yp[mrad]:z[mm]:E[MeV]",
       distrFile_edit";

option_ngenerate_edit
option_nperfile_edit"""

def generate_beam_env(filename, meanE, dist, spread_percent=None, spread_number=None):

    """This command generates a text file describing a gaussian beam with mean = meanE and
        sigma = spread_percent * meanE or sigma = spread_number depending on the user's input. The beam
        position and momenta values used are those defined in filename. If both spread_percent and spread_number
        are defined in the input, it defaults to using the value specified by spread_percent."""

    open_file = filename

    if spread_percent is None and spread_number is not None:

        spread_fraction = spread_number/meanE
        energy_spread = spread_fraction

    elif spread_percent is not None and spread_number is None:

        energy_spread = spread_percent

    elif spread_percent is None and spread_number is None:

        print "WARNING: No energy spread defined. Using preset value of 10%."
        energy_spread = 0.1

    else:

        print "WARNING: Energy spread over-defined. Defaulting to using given percentage spread."
        energy_spread = spread_percent

    particle_number = produce_beam.produce_beam(open_file, meanE, energy_spread, dist)

    return particle_number


def make_gmad(meanE, spread, resultsdir, particle_number):

    file_name = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread) + '.gmad'
    sg = os.path.join(resultsdir, file_name)
    fh = open(sg, "wb")

    distfile_string = 'distrFile="../e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread) +\
                      '/e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread) + '.txt'
    e_beamgmad1 = string.replace(e_beamgmad, 'distrFile_edit', distfile_string)
    ngenerate_string = 'option,ngenerate=' + str(particle_number) + ';'
    e_beamgmad2 = string.replace(e_beamgmad1, 'option_ngenerate_edit', ngenerate_string)
    nperfile_string = '!option,nperfile=' + str(particle_number) + ';'
    e_beamgmad3 = string.replace(e_beamgmad2, 'option_nperfile_edit', nperfile_string)
    fh.write(e_beamgmad3)
    fh.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = """
    This script generates a gaussian beam based on the input file given as the 
    first argument when the script is run. The second argument should define the
    mean energy. The third argument should define whether the energy spread is
    to be defined using a percentage or energy value. The fourth argument should
    then be the desired energy spread value. 
    
    E.g. for generating a beam with mean 5GeV and percentage energy spread of 10%: 
    $ python generate_beam.py 5 percentage 0.1 
    
    for generating a beam with mean 5GeV and energy spread of 0.5GeV
    $ python generate_beam.py 5 value 0.5""",
    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--file', dest = 'file', default=None,
                        help = '''
    This is the path to the example beam .txt file that you are basing the new
    beam description on. It provides the position and momenta values for each
    particle.
    
    E.g. --file "<path to file>"''')

    parser.add_argument('--energy', dest = 'energy', default=None,
                        help = '''
    This defines the mean energy of the gaussian beam that you would like to 
    generate. It is in units of GeV.
    
    E.g. --energy 5.0''')

    parser.add_argument('--spread_type', dest = 'type', default='percentage',
                        choices=['percentage', 'value', 'stepwise'],
                        help = '''
    This defines the type of energy spread that you are defining. The options
    are percentage, value or stepwise. 
    
    percentage corresponds to:
    meanE +/- spread * meanE
    
    value corresponds to:
    meanE +/- spread
    
    value is defined in units of GeV.
    
    stepwise corresponds to producing an equally spaced energy array, ranging 
    from 0 to the value given by the argument --energy. If using this, the 
    value given to the argument '--spread' becomes irrelevant.''')

    parser.add_argument('--spread', dest = 'spread', default=None,
                        help = '''
    This defines the spread of the gaussian beam according to the type given
    in the --spread_type argument. percentage corresponds to a percentage of 
    the mean energy, whereas value corresponds to a particular energy spread
    value measured in GeV.''')

    arguments = parser.parse_args()

    filename = arguments.file

    meanE = float(arguments.energy)

    spread = float(arguments.spread)

    if arguments.type == 'percentage':

        particle_number = generate_beam_env(filename, meanE, dist, spread_percent=spread)

    elif arguments.type == 'value':

        particle_number = generate_beam_env(filename, meanE, dist, spread_number=spread)

    cwd = os.getcwd()

    print "Generating gaussian beam with mean = " + str(meanE) + "GeV and spread = " + str(spread) + "GeV."
    res_dir = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(spread)
    txt_file_name = res_dir + '.txt'
    print "Making directory: ", res_dir
    if os.path.isdir(res_dir) is False:
        os.mkdir(res_dir)
    else:
        print "Directory %s already exists. Stopping." % (res_dir)
        exit()
    make_gmad(meanE, spread, res_dir, particle_number)
    shutil.copy(txt_file_name, res_dir)
    os.remove(txt_file_name)







