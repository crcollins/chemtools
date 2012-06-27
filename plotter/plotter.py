import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plot
import math
np.seterr(all="ignore")

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="This program generates homo/lumo fits and graphs.")
parser.add_argument('files', metavar='file', type=str, nargs='*', help='The name of single file.')
parser.add_argument('-o', metavar='output', action="store", dest="output", type=str, help='The base output folder.')
parser.add_argument('-D', action="store_true", dest="data", default=False, help='Toggles writing a file with the data.')
parser.add_argument('-F', action="store_true", dest="folder", default=False, help='Toggles writing files to a seperate folder.')
parser.add_argument('-E', action="store_true", dest="error", default=False, help='Toggles showing error messages.')

class Output(object):
    def __init__(self, args):
        self.args = args
        self.errors = []

    def parse_file(self, filename):
        out = []
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    out.append([float(x.strip()) for x in line.replace(' ', '').split(',') if x])
        return out

    def write(self, line):
        if self.outputfile:
            self.outputfile.write(line)
        else:
            print line.rstrip('\n')

    def write_file(self):
        for filename in self.args.files:
            if self.args.output:
                dirname = self.args.output
            else:
                dirname = os.path.dirname(filename)
            base, ext = os.path.splitext(os.path.basename(filename))
            if self.args.folder:
                folderpath = os.path.join(dirname, base)
                if not os.path.exists(folderpath):
                    os.makedirs(folderpath)
                base = os.path.join(folderpath, '')

            if self.args.data:
                self.outputfile = open("%sdata.txt" % base, 'w')
            else:
                self.outputfile = None

            try:
                self.run(base, *self.parse_file(filename))
            except Exception as ex:
                self.errors.append(ex)
            if self.outputfile:
                self.outputfile.close()
                self.outputfile = None

        if self.errors:
            self.write("\n---- Errors (%i) ----" % len(self.errors))
            for error in self.errors:
                self.write(str(error))
        if self.outputfile == None or not self.args.error:
            raw_input("<Press Enter>")

    def run(self, filename, datax, datahomo, datalumo, datagap):
        def homofunc(x, a, b):
            return a * np.sqrt(1-b*np.cos(math.pi/(x+1)))

        x = np.array(datax)
        maxx = max(datax)
        if maxx > 1:
            x = 1. / x

        homoy = np.array(datahomo)
        (homoa, homob), var_matrix = curve_fit(homofunc, x, homoy, p0=[-8, -.8])
        self.write("Homo\n")
        self.write("A: %f, B: %f\n" % (homoa, homob))
        self.write("limit: %f\n" % homofunc(0, homoa, homob))
        self.write("\n")

        lumofunc = lambda x,a,b: homofunc(x,a,b) + homofunc(x, homoa, homob)
        lumoy = np.array(datalumo)
        (lumoa, lumob), var_matrix = curve_fit(lumofunc, x, lumoy, p0=[5, -.8])
        self.write("Lumo\n")
        self.write("A: %f, B: %f\n" % (lumoa, lumob))
        self.write("limit: %f\n" % lumofunc(0, lumoa, lumob))
        self.write("\n")

        gapfunc = lambda x,a,b: homofunc(x,a,b)+lumofunc(x, lumoa, lumob)
        gapy = np.array(datagap)
        (gapa, gapb), var_matrix = curve_fit(gapfunc, x, gapy, p0=[11, -.8])
        self.write("Gap\n")
        self.write("A: %f, B: %f\n" % (gapa, gapb))
        self.write("limit: %f" % gapfunc(0, gapa, gapb))

        plot.plot(x, homoy, 'ro')
        plot.plot(np.linspace(0,maxx,20),homofunc(np.linspace(0,maxx,20), homoa, homob),'r')
        plot.plot(x, lumoy, 'ro')
        plot.plot(np.linspace(0,maxx,20),lumofunc(np.linspace(0,maxx,20), lumoa, lumob),'g')

        plot.ylabel("Eg in eV")
        plot.xlabel("1/N")
        plot.savefig('%shomolumo.eps' % filename)

        plot.clf()
        plot.plot(x, gapy, 'ro')
        plot.plot(np.linspace(0,maxx,20),gapfunc(np.linspace(0,maxx,20), gapa, gapb),'r')
        plot.ylabel("Eg in eV")
        plot.xlabel("1/N")
        plot.savefig('%sgap.eps' % filename)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        a = Output(parser.parse_args(sys.argv[1:]))
        a.write_file()
    else:
        args = raw_input('Arguments: ')
        a = Output(parser.parse_args(args.strip().split()))
        a.write_file()
