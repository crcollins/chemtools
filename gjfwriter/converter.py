import os

cores = ['CCC', 'CNN', 'CON', 'CSN', 'TCC', 'TNN', 'TON', 'TSN']
xrgroups = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l']
aryl = ['2', '3', '4', '5', '6', '7', '8', '9']
coreparts = ["*~0", "*~1", "~0", "~1"]
xrparts = ["*+0", "*+1"]
arylparts = ["~0", "~1", "+0", "+1"]

partslist = [coreparts, xrparts, arylparts]

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

	def __str__(self):
		return self.element + " %f %f %f" %(self.xyz)


class Bond(object):
	def __init__(self, atoms, type_, parent=None):
		self.parent = parent

		self.atoms = atoms
		self.type = type_

	@property
	def length(self):
		return sum((x-y)**2 for (x,y) in zip(*tuple(a.xyz for a in self.atoms))) ** .5

	@property
	def id(self):
		return self.parent.index(self)+1

	@property
	def mol2(self):
		return "%d %d %s" %(self.atoms[0].id, self.atoms[1].id, self.type)


def parse_mol2(filename):
	f = open(filename, "r")
	atoms = []
	bonds = []
	state = -2
	for line in f:
		if "@<TRIPOS>" in line:
			state += 1
		elif state == 0:
			x,y,z,e = line.split()[-4:]
			atoms.append(Atom(x,y,z,e, atoms))
		elif state == 1:
			a1, a2, t = line.split()[-3:]
			bonds.append(Bond((atoms[int(a1)-1], atoms[int(a2)-1]), t, bonds))
	f.close()
	return atoms, bonds

added = ["Sg", "Bh", "Hs", "Mt"]

for fname in os.listdir("mol2"):
	name = fname[:-5]
	for i,x in enumerate([cores, xrgroups, aryl]):
		if name in x:
			parts = partslist[i]

	atoms, bonds = parse_mol2("mol2/"+fname)
	f = open("converted/"+name, "w")

	for atom in atoms:
		if atom.element in added:
			atom.element = parts[added.index(atom.element)]
		f.write(str(atom)+"\n")
	f.write("\n")
	for bond in bonds:
		f.write(bond.mol2+"\n")
	f.close()
