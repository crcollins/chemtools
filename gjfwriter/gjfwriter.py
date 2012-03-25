#!/usr/bin/env python

import math
import os
import copy
import argparse
import sys

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

try:
	import Image
	import ImageDraw
	parser.add_argument('-d', type=int, action="store", default=0, help="Used to scale an output image. (0 by default, meaning no picture)")
except:
	pass

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
			#means the molecule was made from read_data()
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

			maxx = max(atom.x, maxx)
			maxy = max(atom.y, maxy)
			maxz = max(atom.z, maxz)
		return (minx, miny, minz), (maxx, maxy, maxz)

	def open_ends(self):
		openbonds = []
		for x in self.bonds:
			if any(True for atom in x.atoms if atom.element[0] in "+*~"):
				openbonds.append(x)
		return openbonds
	
	def next_open(self, conn="~*+"):
		#scans for the first available bond in order of importance.
		bonds = self.open_ends()
		for x in conn:
			for bond in bonds:
				if x in [atom.element[0] for atom in bond.atoms]:
					return bond
		try:
			#check the second bond type
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
		if not Image:
			return
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
			'H': (220, 220, 220),
			'Si': (128, 170, 128)
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
		#rotate to standard view
		img.rotate(-90).save(name+".png")

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

	def merge(self, bond1, bond2, frag):
		#bond1 <= (bond2 from frag)
		#find the part to change
		if bond1.atoms[0].element[0] in "~*+":
			R1, C1 = bond1.atoms
		elif bond1.atoms[1].element[0] in "~*+":
			C1, R1 = bond1.atoms
		else:
			raise Exception(5, "bad bond")
		if bond2.atoms[0].element[0] in "~*+":
			R2, C2 = bond2.atoms
		elif bond2.atoms[1].element[0] in "~*+":
			C2, R2 = bond2.atoms
		else:
			raise Exception(6, "bad bond")

		#saved to prevent overwriting them
		R2x, R2y = R2.x, R2.y
		C1x, C1y = C1.x, C1.y
		#angle of 1 - angle of 2 = angle to rotate
		theta = math.atan2(R1.y-C1.y, R1.x-C1.x) - math.atan2(C2.y-R2.y, C2.x-R2.x)
		frag.rotate(theta, (R2x, R2y), (C1x, C1y))
		
		if bond1.atoms[0].element[0] in "~*+":
			bond1.atoms = (C2, C1)
		else:
			bond1.atoms = (C1, C2)
		#remove the extension parts
		[x.parent.remove(x) for x in (bond2, R1, R2)]

	def chain(self, left, right, n):
		#raise Exception(0, "Chains are not Implemented")
		frags = [copy.deepcopy(self)]
		lidx, ridx = self.bonds.index(left), self.bonds.index(right)
		for i in xrange(n-1):
			a = copy.deepcopy(self)
			if i == 0:
				frags[i].merge(frags[i].bonds[ridx], a.bonds[lidx], a)
			else:
				frags[i].merge(frags[i].bonds[ridx-1], a.bonds[lidx], a)
			frags.append(a)
		return Molecule(frags)

	def stack(self, x, y, z):
		frags = [self]
		bb = self.bounding_box()
		size = tuple(maxv-minv for minv, maxv in zip(bb[0], bb[1]))
		for i,axis in enumerate((x,y,z)):
			#means there is nothing to stack
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

##############################################################################

class Output(object):
	def __init__(self, args):
		self.basis_ = args.basis.replace("/", "_").replace("-", "_").replace("(","").replace(")", "")
		self.dft = "TDDFT" if args.t else "DFT"
		self.errors = []
		self.mem = args.mem
		self.longname = args.longname | args.verbose
		self.error = args.error | args.verbose
		try:
			self.scale = args.d
		except:
			self.scale = 0
		self.args = args

		self.n = '' if args.n <= 1 else "n%i" %(args.n)

		self.x = '' if args.x <= 1 else "x%i" %(args.x)
		self.y = '' if args.y <= 1 else "y%i" %(args.y)
		self.z = '' if args.z <= 1 else "z%i" %(args.z)

		for name in self.get_names():
			try:
				molecule = self.build(*self.parse(name))
				filename = self.get_filename(name)
				try:
					f = open(os.path.join(args.folder, filename+".gjf"), "w")
				except IOError:
					self.errors.append(("", "Bad Folder Name"))
					break
				self.write_file(molecule, f, filename)
				if self.scale:
					molecule.draw(filename, self.scale)
			except Exception as message:
					self.errors.append(message)
			print name, "---- Done"
		

		if self.error:
			print "\n---- Errors ----"
			for x in self.errors:
				if type(x) == tuple:
					print " - ".join([str(x[0]), x[1]])
				else:
					print repr(x)

	def get_filename(self, name):
		if self.longname:
				nameparts = [x for x in [name, self.basis_, self.x, self.y, self.z, self.n, self.dft] if x]
		else:
			nameparts = [x for x in [name, self.x, self.y, self.z, self.n, self.dft] if x]
		return "_".join(nameparts)

	def get_names(self):
		if self.args.listfiles:
			that = [self.read_listfile(x) for x in self.args.listfiles] 
			this = [y for y in that if y]
			other = sum([x for x in this if x],[])
			names = set(self.args.names + other)
		else:
			names = set(self.args.names)
		return names

	def write_file(self, molecule, f, name):
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
		core = (Molecule(read_data(corename)), 0, corename, corename)
		struct = [x if x else "A" for x in [rightname, leftname] + [middlename] * 2]
		midx = 2

		fragments = []
		for side in struct:
			this = []
			for i, char in enumerate(side):
				if char in "10":
					raise Exception(4, "10 not implemented")
				this.append((Molecule(read_data(char)), i, char, side))
				if char.islower() and i == len(side)-1 :
					#checks to make sure the prev group was aryl for r-groups
					if this[i-1][2] not in '2389' and this[i-1][2].isdigit():
						this.append((Molecule(read_data(char)), i, char, side))
			fragments.append(this)
		ends = []
		#bond all of the fragments together
		for j, side in enumerate(fragments):
			this = [core]+side
			for (part, partidx, char, sidename) in side:
				bondb = part.next_open()
				c = bondb.connection()
				#enforces lowercase to be r-group
				if char.islower():
					c = "+"
				if j >= midx and c == "~":
					c = "*" + c
				#fixes x-groups connecting to ~
				if j >= midx and c == "*" and char.isupper():
					c = "~" + c
				bonda = this[partidx][0].next_open(c)
				if bonda and bondb:
					this[partidx][0].merge(bonda, bondb, part)
			ends.append(this[max([x[1] for x in side])][0].next_open('~'))

		#merge the fragments into single molecule
		out = [core[0]]
		for x in fragments:
			for y in x:
				out.append(y[0])
		a = Molecule(out)

		#multiplication of molecule/chain
		limit = all(ends[:2]) and all('~' in x.connection() for x in ends[:2])
		if self.args.n > 1 and limit:
			a = a.chain(ends[0], ends[1], self.args.n)
		if any((self.args.x, self.args.y, self.args.z)):
			a = a.stack(self.args.x, self.args.y, self.args.z)
		a.close_ends()
		return a

if __name__ == '__main__':
	if len(sys.argv) > 1:
		out = Output(parser.parse_args(sys.argv[1:]))
	else:
		args = raw_input('Arguments: ')
		out = Output(parser.parse_args(args.strip().split()))
		paw = raw_input("<Press Enter>")
