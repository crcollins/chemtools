#!/usr/bin/python
import operator
import os
import argparse
import sys

messages = []
for name in ("bjob.txt", "tjob.txt", "gjob.txt"):
    messages.append(''.join(open("data/"+name, "r").readlines())) 

parser = argparse.ArgumentParser(description="This program writes Gaussian job files from molecule names.")
parser.add_argument('names', metavar='name', type=str, nargs='*', default=list(), help='The name of the molecule to create.')
parser.add_argument('-i', metavar='list_file', action="store", nargs='*', default=list(), type=str, dest="listfiles", help='A file with a listing of molecules to make job files for.')
parser.add_argument('-o', metavar="outfolder", action="store", default='.', dest="outfolder", type=str, help='The folder to output all the job files.')
parser.add_argument('-f', metavar='folder', action="store", nargs='*', default=list(), dest="folders", type=str, help='A folders to read molecules from.')

parser.add_argument('-n', type=int, action="store", dest="nodes", default=1, help="The amount of nodes to use for the calculation.(1 node by default)(This will be converted for Blacklight jobs.)")

parser.add_argument('-e', action="store", default="none@domain.com", dest="email", help="The email to send the callback to. (none@domain.com by default)")
parser.add_argument('-t', type=int, action="store", default=16, dest="time",  help="The amount of time, in hours, to use for the calculation . (16 hours by default)")

parser.add_argument('-E', action="store_true", dest="error", default=False, help='Toggles showing error messages.')
parser.add_argument('-V', action="store_true", dest="verbose", default=False, help='Toggles showing all messages.')
parser.add_argument('-B', action="store_true", dest="blacklight", default=False, help='Toggles writing Blacklight job files.')
parser.add_argument('-T', action="store_true", dest="trestles", default=False, help='Toggles writing Trestles job files.')
parser.add_argument('-G', action="store_true", dest="gordon", default=False, help='Toggles writing Gordon job files.')
parser.add_argument('-A', action="store_true", dest="all", default=False, help='Toggles writing all job files.')


class Output(object):
    def __init__(self, args):
        self.errors = []
        self.error = args.error | args.verbose
        self.email = args.email
        self.time = args.time
        self.names = [os.path.splitext(x)[0] for x in args.names]
        self.bjob = args.blacklight | args.all
        self.tjob = args.trestles | args.all
        self.gjob = args.gordon | args.all
        self.outfolder = args.outfolder
        self.nodes = args.nodes

        self.folders = args.folders
        for folder in self.folders:
            try:
                self.names += [os.path.splitext(x)[0] for x in os.listdir(folder) if os.path.isfile(x)]
            except:
                self.errors.append((folder, "Bad Folder Name"))
        self.listfiles = args.listfiles
        for filename in self.listfiles:
            try:
                f = open(filename, "r")
                self.names += [os.path.splitext(line.strip())[0] for line in f]
            except:
                self.errors.append((filename, "Bad File Name"))
        self.args = args

        j = ["b","t","g"]
        if not os.path.isdir(self.outfolder):
            self.errors.append((self.outfolder, "Bad Output Folder Name"))
        else:
            for name in self.names:
                try:
                    self.write_file(name)
                except Exception as (num, message, ):
                    self.errors.append((name, message))
                print name, "---- Done"
        
        if self.error:
            print "\n---- Errors ----"
            print "\n".join([" - ".join(x) for x in self.errors])

    def write_file(self, name):
        j = ["b","t","g"]
        factor = self.get_time_factor(name)
        time = "%d:00:00" %(self.time*factor)
        for i, x in enumerate([self.bjob, self.tjob, self.gjob]):
            if not x:
                continue
            f = open(os.path.join(self.outfolder, name+".%sjob"%j[i]), "w")
            use = { 
                    "name":name,
                    "email":self.email,
                    "nodes":self.nodes,
                    "ncpus":self.nodes*16,
                    "time":time
                    }
            f.write(messages[i].format(**use))
    
    def get_time_factor(self, name):
        s = name.split("_")
        s1 = [int(x[1:]) for x in s if (x[0] in "xyzn") and x[1].isdigit()]
        return reduce(operator.mul, s1) if s1 else 1


if __name__ == '__main__':
    if len(sys.argv) > 1:
        out = Output(parser.parse_args(sys.argv[1:]))
    else:
        args = raw_input('Arguments: ')
        out = Output(parser.parse_args(args.strip().split()))
        paw = raw_input("<Press Enter>")
