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


def produce_strange_beam(peak, gauss_wid_low, gauss_wid_high, n):

    """ This function generates a lopsided beam with gaussians of different
    widths on either side of the peak energy."""

    beam = []
    # frac = gauss_wid_high/gauss_wid_low
    frac = 0.5

    num_part_low = int((1 - frac) * n)
    lower_part = abs(np.random.normal(0, gauss_wid_low*peak, num_part_low))

    for i in range(0, len(lower_part)):
        beam.append(abs(peak - lower_part[i]))

    num_part_high = int(frac * n)
    upper_part = abs(np.random.normal(0, gauss_wid_high*peak, num_part_high))

    for j in range(0, len(upper_part)):
        beam.append(peak + upper_part[j])

    return beam


def generate_beam_energies(meanE, energy_spread, energy_spread_low, array,
                           dist):

    """This function generates a sample of energies defined by the chosen mean energy and energy spread
        for a chosen number of particles defined by the length of the array that is given to the function.
        The default distribution is stepwise. This produces an equally-spaced array with energies ranging
        from 0 to the meanE. Other option is gaussian. """

    particle_number = len(array)

    if dist == "gaussian":

        energy_values = abs(np.random.normal(meanE, energy_spread*meanE,
                                         particle_number))

    elif dist == "lopsided":

        energy_values = produce_strange_beam(meanE, energy_spread_low,
                                             energy_spread, particle_number)


    else:

        energy_value_start = np.linspace(meanE - energy_spread*meanE, meanE +
                                         energy_spread*meanE,
                                         particle_number)

        energy_values = []

        for j in energy_value_start:

            energy = abs(np.random.normal(j, 0.01, 1))
            energy_values.append(energy)

    return energy_values


def replace_energies(array, energy_values):

    """This function replaces the original energy values from the file given to read_file with those generated
        by generate_beam_energies."""

    for i in range(0, len(array)):

        new_energy = energy_values[i]
        array[i][5] = new_energy

    return array


def produce_beam(filename, meanE, energy_spread, energy_spread_low, dist):

    """This function generates a text file of particle energies of a given mean energy and spread with
        position and momentum distributions according to the original file that it reads. It replaces the
        supplied energy values with a Gaussian distribution centred on meanE with energy spread defined via:
            meanE +/- spread_percent * meanE"""

    array = read_file(filename)

    energy_GeV_to_MeV = meanE*1000

    energy_values = generate_beam_energies(energy_GeV_to_MeV, energy_spread, array, dist)

    output = replace_energies(array, energy_values)

    if dist == 'gaussian':

        output_file_name = 'e_beam_mean_' + str(meanE) + '_spread_percentage_'\
                           + str(energy_spread) + '.txt'

    elif dist == "lopsided":

        output_file_name = 'e_beam_mean_' + str(meanE) + \
                           '_upperspread_' + str(energy_spread) + \
                           '_lowerspread' + str(energy_spread_low) + '.txt'

    else:

        output_file_name = 'e_beam_max_' + str(meanE) + \
                           '_stepwise.txt'

    np.savetxt(output_file_name, output, fmt='%.3f')

    return len(array)


