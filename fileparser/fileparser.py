#!/bin/python

import argparse
import os
import sys
import time

parser = argparse.ArgumentParser(description="This program extracts data from Gaussian log files.")
parser.add_argument('files', metavar='file', type=str, nargs='*', help='The name of single file.')
parser.add_argument('-i', metavar='list_file', action="store", nargs='*', dest="listfiles", type=str, help='A file with a listing of other files.')
parser.add_argument('-f', metavar='folder', action="store", nargs='*', dest="folders", type=str, help='A folder with a collection of files.')
parser.add_argument('-o', metavar='output', action="store", dest="outputfile", type=str, help='The output file.')
parser.add_argument('-e', action="store_true", dest="error", default=False, help='Toggles showing error messages.')
parser.add_argument('-p', action="store_true", dest="paths", default=False, help='Toggles showing paths to files.')
parser.add_argument('-R',  action="store_true", dest="rel", default=False, help='Toggles showing relative paths.')
parser.add_argument('-v', action="store_true", dest="verbose", default=False, help='Toggles showing all messages.')

class Output(object):
    def __init__(self, args):
        self.errors = []
        self.outputfilename = self.check_output_file(args.outputfile)
        self.error = args.error | args.verbose
        self.paths = args.paths | args.verbose | args.rel
        self.rel = args.rel
        self.files = self.check_input_files(args.files
                            + self.convert_files(args.listfiles)
                            + self.convert_folders(args.folders))

    def check_input_files(self, filelist):
        files = []
        for x in filelist:
            if not os.path.isfile(x):
                path = os.path.relpath(x) if self.rel else os.path.abspath(x)
                self.errors.append("Invalid filename:  '" + path + "'")
            else:
                files.append(x)
        return files

    def check_output_file(self, filename):
        if filename and os.path.isfile(filename):
            oldname = filename
            filename = "output_%i.txt"% time.time()
            self.errors.append("'" + oldname + "' already in use, using: 'output_%i.txt'"% time.time())
        return filename

    def convert_files(self, filenames):
        if filenames:
            files = list()
            for filename in filenames:
                if os.path.isfile(filename):
                    f = open(filename, 'r')
                    files += [x.strip() for x in f if x.strip()]
                    f.close()
            return files
        else:
            return list()

    def convert_folders(self, folders):
        if folders:
            files = list()
            for folder in folders:
                if os.path.isdir(folder):
                    path = os.path.relpath(folder) if self.rel else os.path.abspath(folder)
                    files += [os.path.join(path,x) for x in os.listdir(folder)]
            return files
        else:
            return list()

    def find_lines(self, filename):
        f = open(filename, 'r')
        flines = f.readlines()
        occline, virtualline, excitedline, timeline = ['']*4
        try:
            if "Entering Gaussian System" not in flines[0]:
                return occline, virtualline, excitedline, timeline
        except IndexError:
            return occline, virtualline, excitedline, timeline
    
        occbool = False
        skip = 300
        L = len(flines)
        for i, line in enumerate(flines[skip:]):
            if occline and virtualline and excitedline:
                break
            elif "1:" in line:
                excitedline = line
            elif "occ. eigenvalues" in line:
                occbool = True
            elif "virt. eigenvalues" in line and occbool:
                virtualline = line
                occline = flines[(i-1+skip)%L]
                occbool = False
            else:
                occbool = False
        for line in flines[::-1]:
            if 'Job cpu time' in line:
               timeline=line
               break
        f.close()
        return occline, virtualline, excitedline, timeline

    def clean_lines(self, (occline, virtualline, excitedline, timeline)):
        occ = occline.strip().split()[-1]
        try:
            virtual = virtualline.strip().split()[4]
        except:
            virtual = "---"
        try:
            excited = excitedline.strip().split()[4]
            float(excited)
            assert excited != "09"
        except:
            excited = "---"
        time = self.convert_time(timeline.strip().split()[3:-1][0::2])
        occ, virtual = self.convert_values((occ, virtual))
        return occ, virtual, excited, time

    def convert_time(self, time):
        con = (24., 1., 1/60., 1/3600)
        return str(sum(float(x)*con[i] for i, x in enumerate(time)))

    def convert_values(self, values):
        con = (27.2117, 27.2117)
        return tuple(str(float(x)*con[i]) if x != "---" else "---" for i, x in enumerate(values))

    def path(self, filename):
        return os.path.relpath(filename) if self.rel else os.path.abspath(filename)
        
    def filename(self, path, name):
        return path if self.paths else os.path.splitext(os.path.basename(name))[0]

    def write_file(self):
        if self.outputfilename:
            self.outputfile = open(self.outputfilename,'w')
        else:
            self.outputfile = None
        self.write("Filename, Occ, Virtual, Excited, Time")
        for filename in self.files:
            lines = self.find_lines(filename)
            if not any(lines):
                path = self.path(filename)
                self.errors.append("Invalid file type:  '" + path + "'")
                continue
            elif all(lines) or all(lines[:2] + tuple(lines[3])):    
                ovft = self.clean_lines(lines)
            elif any(lines):
                lines = [x if x else "---" for x in lines]
                ovft = self.clean_lines(lines)
            path = self.path(filename)
            filename = self.filename(path, filename)
            self.write(', '.join([x for x in (filename,) + ovft if x]))
        if self.error:
            self.write("\n---- Errors (%i) ----" % len(self.errors))
            for error in self.errors:
                self.write(error)
        if self.outputfile == None:
            raw_input("<Press Enter>")
        
    def write(self, line):
        if self.outputfile:
            self.outputfile.write(line+'\n')
        else:
            print line

if __name__ == '__main__':
    if len(sys.argv) > 1:
        a = Output(parser.parse_args(sys.argv[1:]))
        a.write_file()
    else:
        args = raw_input('Arguments: ')
        a = Output(parser.parse_args(args.strip().split()))
        a.write_file()
