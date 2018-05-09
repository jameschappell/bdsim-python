import os
import glob
import argparse
import shutil

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

    # find all directories within current directory

    string_1 = cwd + '/' + arguments.file

    dir_list = [name for name in os.listdir(string) if os.path.isdir(
        os.path.join(string, name))]

    length_dir = len(dir_list)

    k = 1

    for j in range(0, length_dir):

        print 'Entering directory %i of %i' % (k, length_dir)

        chdir_name = string_1 + '/' + dir_list[j]

        os.chdir(chdir_name)

        cwd_2 = os.getcwd()

        # find all .root files in the target directory

        file_list = find_files(cwd_2, 'root')

        length = len(file_list)

        i = 1

        # apply the analysis script to extract the relevant data

        for file_name in file_list:

            print 'Analysing file %i of %i total files' % (i, length)

            output_file_name = os.path.splitext(file_name)[0]

            output_file = cwd + '/' + output_file_name + '.txt'

            string_2 = cwd_2 + '/' + file_name

            execute_string = '$QUAD_SCAN_ANALYSIS/data_extract ' + string_2 + \
                             ' > ' +\
                             output_file

            os.system(execute_string)

            i = i + 1

        k = k + 1

    # Compress output

    tar_name = 'quadscan_output_for_analysis.tar.bz2'

    tar_command = 'tar jcvf ' + tar_name + ' x*.txt'

    os.system(tar_command)





