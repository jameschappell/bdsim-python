import numpy as np


def read_file(filename):

    """This function reads the example beam file given by filename into the program and creates an array
       of the data."""

    open_file = filename
    #f = open(open_file, 'r')
    #g = f.readlines()

    array = []

    with open(open_file) as f:

        for line in f:

            array_numbers = map(float, line.split())
            array.append(array_numbers)

    return array


def generate_beam_energies(meanE, energy_spread, array, dist):

    """This function generates a sample of energies defined by the chosen mean energy and energy spread
        for a chosen number of particles defined by the length of the array that is given to the function.
        The default distribution is stepwise. This produces an equally-spaced array with energies ranging
        from 0 to the meanE. Other option is gaussian. """

    particle_number = len(array)

    if dist == "gaussian":

        energy_values = np.random.normal(meanE, energy_spread*meanE,
                                         particle_number)

    else:

        energy_values = np.linspace(0., meanE, particle_number)

    return energy_values


def replace_energies(array, energy_values):

    """This function replaces the original energy values from the file given to read_file with those generated
        by generate_beam_energies."""

    for i in range(0, len(array)):

        new_energy = energy_values[i]
        array[i][5] = new_energy

    return array


def produce_beam(filename, meanE, energy_spread, dist):

    """This function generates a text file of particle energies of a given mean energy and spread with
        position and momentum distributions according to the original file that it reads. It replaces the
        supplied energy values with a Gaussian distribution centred on meanE with energy spread defined via:
            meanE +/- spread_percent * meanE"""

    array = read_file(filename)

    energy_GeV_to_MeV = meanE*1000

    energy_values = generate_beam_energies(energy_GeV_to_MeV, energy_spread, array, dist)

    output = replace_energies(array, energy_values)

    output_file_name = 'e_beam_mean_' + str(meanE) + '_spread_percentage_' + str(energy_spread) + '.txt'

    np.savetxt(output_file_name, output, fmt='%.3f')

    return len(array)





