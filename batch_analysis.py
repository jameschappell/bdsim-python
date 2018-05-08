import os
import glob
import argparse



def find_files(dir):

    os.chdir(dir)

    file_list = []

    for file_name in glob.glob("*.root"):

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

    file_list = find_files(string)

    print file_list





