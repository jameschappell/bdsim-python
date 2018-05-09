import os
import glob
import argparse


def find_files(dir, ext):

    os.chdir(dir)

    file_list = []

    ext_string = '*.' + ext

    for file_name in glob.glob(ext_string):

        file_list.append(file_name)

    return file_list


if __name__ == "__main__":

    cwd = os.getcwd()

    os.chdir(cwd)

    parser = argparse.ArgumentParser(description="""
           This script generates text files that access ROOT data files and 
           return the x- and y-positions of particles hitting the AWAKE 
           spectrometer screen. It uses a C++ script to access the event 
           data that was written by David Cooke (UCL).""",
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)

    parser.add_argument('--file', dest='file', default=None,
                        help='''
           This is the directory in which the ROOT files reside.''')

    arguments = parser.parse_args()

    file_loc = arguments.file

    string = cwd + '/' + arguments.file

    # find all .root files in the target directory

    file_list = find_files(string, 'root')

    length = len(file_list)

    i = 1

    # apply the analysis script to extract the relevant data

    for file_name in file_list:

        print 'Analysing file %i of %i total files' % (i, length)

        output_file_name = os.path.splitext(file_name)[0]

        output_file = cwd + '/' + output_file_name + '.txt'

        string = cwd + '/' + file_name

        execute_string = '$QUAD_SCAN_ANALYSIS/data_extract ' + string + ' > ' +\
                         output_file

        os.system(execute_string)

        i = i + 1

    # Compress output

    tar_name = 'quadscan_output_for_analysis.tar.bz2'

    tar_command = 'tar jcvf ' + tar_name + ' x*.txt'

    os.system(tar_command)





