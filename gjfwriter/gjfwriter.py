#!/usr/bin/env python

import math
import os
import copy
import argparse
import sys
import Image
import ImageDraw

parser = argparse.ArgumentParser(description="This program writes Gaussian .gjf files from molecule names.")
parser.add_argument('names', metavar='name', type=str, nargs='*', default=list(), help='The name of the molecule to create.')
parser.add_argument('-i', metavar='list_file', action="store", nargs='*', default=list(), dest="listfiles", type=str, help='A file with a listing of molecules to make.')
parser.add_argument('-f', metavar='folder', action="store", default=".", dest="folder", type=str, help='A folder to output the files.')
parser.add_argument('-n', type=int, action="store", default=1, help="The length of the chain. (NOT IMPLEMENTED)")

parser.add_argument('-x', type=int, action="store", default=1, help="The amount of molecules to stack on the x axis.")
parser.add_argument('-y', type=int, action="store", default=1, help="The amount of molecules to stack on the y axis.")
parser.add_argument('-z', type=int, action="store", default=1, help="The amount of molecules to stack on the z axis.")

parser.add_argument('-b', action="store", dest="basis", default="b3lyp/6-31g(d)", help="The basis functional to use for the calculation. (b3lyp/6-31g(d) by default)")
parser.add_argument('-m', action="store", dest="mem", default ="59GB",  help="The amount of memory to use for the calculation. (59GB by default)")

parser.add_argument('-d', type=int, action="store", default=0, help="Used to scale an output image. (0 by default, meaning no picture)")

parser.add_argument('-T', action="store_true", dest="t", default=False, help="Toggles to use the TDDFT method.")
parser.add_argument('-E', action="store_true", dest="error", default=False, help='Toggles showing error messages.')
parser.add_argument('-V', action="store_true", dest="verbose", default=False, help='Toggles showing all messages.')
parser.add_argument('-L', action="store_true", dest="longname", default=False, help='Toggles showing the long name.')

DB = [x for x in os.listdir("data")]
CORES = [x for x in DB if len(x) == 3]
OTHERS = [x for x in DB if len(x) == 1]

##############################################################################

class Atom(object):
	def __init__(self, x, y, z, element, parent=None):
		self.parent = parent
		self.element = element
		self.x, self.y, self.z = float(x),float(y),float(z)

	@property
	def xyz(self):
		return self.x, self.y, self.z

	@property
	def id(self):
		return self.parent.index(self)+1

	@property
	def mol2(self):
		return "{0} {1}{0} {2} {3} {4} {1}".format(self.id, self.element, *self.xyz) 

	@property
	def gjf(self):
		return self.element + " %f %f %f" %(self.xyz)

	def __str__(self):
		return self.element + " %f %f %f" %(self.xyz)


class Bond(object):
	def __init__(self, atoms, type_, parent=None):
		self.parent = parent

		self.atoms = atoms
		self.type = type_

	def connection(self):
		if self.atoms[0].element[0] in "~*+":
			b = self.atoms[0].element[:2]
		else:
			b = self.atoms[1].element[:2]
		return ''.join([x for x in b if x in "~*+"])

	@property
	def length(self):
		return sum((x-y)**2 for (x,y) in zip(*tuple(a.xyz for a in self.atoms))) ** .5

	@property
	def id(self):
		return self.parent.index(self)+1

	@property
	def mol2(self):
		return "%d %d %d %s" %(self.id, self.atoms[0].id, self.atoms[1].id, self.type)


class Molecule(object):
	def __init__(self, fragments):
		self.atoms = []
		self.bonds = []
		self.clean_input(fragments)
	
	def clean_input(self, fragments):
		try:
			for frag in fragments:
				for atom in frag.atoms:
					atom.parent = self.atoms
					self.atoms.append(atom)
				for bond in frag.bonds:
					bond.parent = self.bonds
					self.bonds.append(bond)
		except AttributeError:
			for atom in fragments[0]:
				atom.parent = self.atoms
				self.atoms.append(atom)
			for bond in fragments[1]:
				bond.parent = self.bonds
				self.bonds.append(bond)

	def rotate(self, theta, point, offset):
		#point = (x,y); offset = (x,y)
		ct = math.cos(theta)
		st = math.sin(theta)
		for atom in self.atoms:
			x = atom.x - point[0]
			y = atom.y - point[1]
			atom.x = ct*x - st*y + offset[0]
			atom.y = st*x + ct*y + offset[1]

	def displace(self, x, y, z):
		for atom in self.atoms:
			atom.x += x
			atom.y += y
			atom.z += z

	def bounding_box(self):
		minx, miny, minz = self.atoms[0].xyz
		maxx, maxy, maxz = self.atoms[0].xyz
		for atom in self.atoms[1:]:
			minx = min(atom.x, minx)
			miny = min(atom.y, miny)
			minz = min(atom.z, minz)

			maxx = max(atom.x, minx)
			maxy = max(atom.y, miny)
			maxz = max(atom.z, minz)
		return (minx, miny, minz), (maxx, maxy, maxz)

	def open_ends(self):
		openbonds = []
		for x in self.bonds:
			if any(True for atom in x.atoms if atom.element[0] in "+*~"):
				openbonds.append(x)
		return openbonds
	
	def next_open(self, conn="~*+"):
		bonds = self.open_ends()
		for x in conn:
			for bond in bonds:
				if x in [atom.element[0] for atom in bond.atoms]:
					return bond
		try:
			for x in conn:
				for bond in bonds:
					if x in [atom.element[1] for atom in bond.atoms if len(atom.element)>1]:
						return bond
		except:
			pass
		
	def close_ends(self):
		for atom in self.atoms:
			if atom.element[0] in "~*+":
				atom.element = "H"

	def draw(self, name, scale):
		colors = {
			'1': (255,255,255),
			'Ar': (255, 0, 0),
			'2': (0, 255, 0),
			'3': (0, 0, 255),
			'S': (255, 255, 0),
			'O': (255, 0, 0),
			'N': (0, 0, 255),
			'Cl': (0, 255, 0),
			'Br': (180, 0, 0),
			'C': (128, 128, 128),
			'H': (220, 220, 220)
			}

		bounds = self.bounding_box()
		xres = int(scale * abs(bounds[0][0] - bounds[1][0]))
		yres = int(scale * abs(bounds[0][1] - bounds[1][1]))
		img = Image.new("RGB", (xres, yres))
		draw = ImageDraw.Draw(img)
		for bond in self.bonds:
			pts = sum([x.xyz[:2] for x in bond.atoms], tuple())
			pts = [(coord-bounds[0][i%2])*scale for i,coord in enumerate(pts)]

			draw.line(pts,fill=colors[bond.type])
			s = (scale * .25)
			for x in xrange(2):
				if bond.atoms[x].element not in "C":
					circle = (pts[x*2] - s,pts[x*2+1] - s, pts[x*2] + s, pts[x*2+1] + s)
					draw.ellipse(circle, fill=colors[bond.atoms[x].element])
		img.save(name+".png")

	def __getitem__(self, key):
		for x in self.bonds:
			if key in [y.element for y in x.atoms]:
				return x
		else:
			raise KeyError(key)

	@property
	def mol2(self):
		string = """@<TRIPOS>MOLECULE\nMolecule Name\n%d %d\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM""" %(len(self.atoms), len(self.bonds))
		string += "\n".join([x.mol2 for x in self.atoms] + 
						["@<TRIPOS>BOND", ] + 
						[x.mol2 for x in self.bonds])
		return string
	
	@property
	def gjf(self):
		string = "\n".join([x.gjf for x in self.atoms]) + "\n\n"
		bonds = self.bonds[:]
		for atom in self.atoms:
			s = str(atom.id) + " "
			for bond in tuple(bonds):
				if atom.id in [x.id for x in bond.atoms]:
					try:
						s += [str(x.id) + " " + (bond.type if bond.type != "Ar" else "1.5") + " " for x in bond.atoms if x.id != atom.id][0]
						bonds.remove(bond)
					except:
						pass
			string += s + "\n"
		return string

##############################################################################

def read_data(filename):
	try:
		f = open(os.path.join("data",filename), "r")
	except:
		try:
			f = open(os.path.join("data",filename.lower()), "r")
		except:
			raise Exception(3, "Bad Substituent Name: %s" %filename)
	atoms = []
	bonds = []
	state = 0
	for line in f:
		if line == "\n":
			state = 1
		elif state == 0:
			e,x,y,z = line.split()[-4:]
			atoms.append(Atom(x,y,z,e, atoms))
		elif state == 1:
			a1, a2, t = line.split()
			bonds.append(Bond((atoms[int(a1)-1], atoms[int(a2)-1]), t, bonds))
	f.close()
	return atoms, bonds

def merge(bond1, bond2, frag):
	#bond1 <= (bond2 from frag)
	if bond1.atoms[0].element[0] in "~*+":
		R1, C1 = bond1.atoms
	else:
		C1, R1 = bond1.atoms
	if bond2.atoms[0].element[0] in "~*+":
		R2, C2 = bond2.atoms
	else:
		C2, R2 = bond2.atoms

	R2x, R2y = R2.x, R2.y
	C1x, C1y = C1.x, C1.y
	theta = math.atan2(R1.y-C1.y, R1.x-C1.x) - math.atan2(C2.y-R2.y, C2.x-R2.x)
	frag.rotate(theta, (R2x, R2y), (C1x, C1y))
	
	if bond1.atoms[0].element[0] in "~*+":
		bond1.atoms = (C2, C1)
	else:
		bond1.atoms = (C1, C2)
	[x.parent.remove(x) for x in (bond2, R1, R2)]

def chain(frag, n):
	raise Exception(0, "Chains are not Implemented")

def stack(frag, x, y, z):
	frags = [frag]
	bb = frag.bounding_box()
	size = tuple(f-i for i, f in zip(bb[0], bb[1]))
	for i,axis in enumerate((x,y,z)):
		if axis <= 1:
			continue
		axisfrags = copy.deepcopy(frags)
		for num in xrange(1,axis):
			use = [0,0,0]
			use[i] = num*(2+size[i]) 
			for f in axisfrags:
				a = copy.deepcopy(f)
				a.displace(*use)
				frags.append(a)
	return Molecule(frags)

##############################################################################

class Output(object):
	def __init__(self, args):
		self.basis_ = args.basis.replace("/", "_").replace("-", "_").replace("(","").replace(")", "")
		self.dft = "TDDFT" if args.t else "DFT"
		self.errors = []
		self.mem = args.mem
		self.longname = args.longname | args.verbose
		self.error = args.error | args.verbose
		self.scale = args.d
		self.args = args

		self.n = '' if args.n <= 1 else "n%i" %(args.n)

		self.x = '' if args.x <= 1 else "x%i" %(args.x)
		self.y = '' if args.y <= 1 else "y%i" %(args.y)
		self.z = '' if args.z <= 1 else "z%i" %(args.z)
		
		if args.listfiles:
			this = [y for y in [self.read_listfile(x) for x in args.listfiles] if y]
			self.names = set(args.names + sum([x for x in this if x],[]))
		else:
			self.names = set(args.names)

		for name in self.names:
			try:
				molecule = self.build(*self.parse(name))
				if self.longname:
					nameparts = [x for x in [name, self.basis_, self.x, self.y, self.z, self.n, self.dft] if x]
				else:
					nameparts = [x for x in [name, self.x, self.y, self.z, self.n, self.dft] if x]
				filename = "_".join(nameparts)
				try:
					f = open(os.path.join(args.folder, filename+".gjf"), "w")
				except IOError:
					self.errors.append(("", "Bad Folder Name"))
					break
				self.write_file(molecule, f, filename)
				if self.scale:
					molecule.draw(filename, self.scale)
			except Exception as (num, message, ):
				self.errors.append((name, message))
			print name, "---- Done"
		
		if self.error:
			print "\n---- Errors ----"
			print "\n".join([" - ".join(x) for x in self.errors])

	def write_file(self, molecule, f, name):
		if self.args.n > 1:
			base = copy.deepcopy(molecule)
			molecule = chain(molecule, self.args.n)

		if any((self.args.x, self.args.y, self.args.z)):
			molecule = stack(molecule, self.args.x, self.args.y, self.args.z)
		starter = [
					"%%mem=%s"%self.mem,
					"%%chk=%s.chk" %name,
					"# opt %s geom=connectivity" %self.args.basis,
					"",
					name,
					"",
					"0 1",
					""
					]
		f.write("\n".join(starter))
		f.write(molecule.gjf)

	def read_listfile(self, filename):
		try:
			names = []
			f = open(filename, "r")
			for line in f:
				names.append(line.strip())
			return names
		except:
			self.errors.append((filename, "File Does Not Exist"))

	def parse(self, name):
		parts = name.split("_")
		core = None

		for part in parts:
			if part.upper() in CORES :
				core = part
				break
		if not core:
			raise Exception(1, "Bad Core Name")
		i = parts.index(core)
		left = parts[:i][0] if parts[:i] else None
		right = parts[i+1:]
		if len(right) > 1:
			middle = right[0]
			right = right[1]
		else:
			try:
				letter = right[0][0]
				if letter.lower() in DB and letter.lower() != letter:
					middle = letter
					right = right[0][1:]
				else:
					middle = None
					right = right[0]
			except:
				middle = None
		return core, left, middle, right

	def build(self, corename, leftname, middlename, rightname):
		core = Molecule(read_data(corename))
		structure = [x for x in [rightname, leftname] + [middlename] * 2 if x != None]
		if middlename:
			midx = structure.index(middlename)
		else:
			midx = len(structure)

		fragments = []
		for side in structure:
			this = []
			for i,char in enumerate(side):
				this.append(Molecule(read_data(char)))
				if char in OTHERS and i == (len(side)-1):
					this.append(Molecule(read_data(char)))
			fragments.append(this)

		for j, side in enumerate(fragments):
			this = [core]+side
			for i,part in enumerate(side):
				bondb = part.next_open()
				c = bondb.connection()
				#Used to fix X-groups
				if j >= midx and c == "~":
					c = "*" + c
				#Used to fix R-groups
				if i == (len(side)-1) and c == "*+":
					bonda = this[i-1].next_open(c)
				else:
					bonda = this[i].next_open(c)

				if bonda and bondb:
					merge(bonda, bondb, part)

		a = Molecule([core]+[item for sublist in fragments for item in sublist])
		a.close_ends()
		return a

if __name__ == '__main__':
	if len(sys.argv) > 1:
		out = Output(parser.parse_args(sys.argv[1:]))
	else:
		args = raw_input('Arguments: ')
		out = Output(parser.parse_args(args.strip().split()))
		paw = raw_input("<Press Enter>")
	#for letter in "abcdefgh":
	#	for x in xrange(1,5):
	#		out = Output(parser.parse_args(["24%s_TON_2J_24%s" %(letter, letter), '-o', '../molecules/', '-m', str(x)]))
	#		out.write_file()