import numpy as np
import matplotlib.pyplot as plt

def get_cms(coords, weights):

	mol_mass = np.sum(weights)
	x_sum = 0	#
	y_sum = 0	# weighted sums of each component
	z_sum = 0	#
	
	for z in zip(coords, weights):
		x_sum += z[0][0] * z[1]
		y_sum += z[0][1] * z[1]
		z_sum += z[0][2] * z[1]

	xcms = x_sum / mol_mass
	ycms = y_sum / mol_mass
	zcms = z_sum / mol_mass
	return xcms, ycms, zcms


def convert_coords(mol_coords, atom_weights, new_cms):
	"""
	Converts coordinates of molecules atoms according to coordinates 
	of a new center of masses.

	mol_coords: an array of 3-Tupel with coordinates of each atom in the molecule
	atom_weight>: an array of atom weights (in u) for each 3-Tupel in mol_coords
	new_cms: Coordinates (as a 3-Tupel) of a new center of masses
	"""
	x_cms, y_cms, z_cms = get_cms(mol_coords, atom_weights)	

	r_N = []	# array of radius-vektors of each atom relative to the cms

	for at in mol_coords:
		r_x = round(at[0] - x_cms,4)
		r_y = round(at[1] - y_cms,4)
		r_z = round(at[2] - z_cms,4)
		r_N.append((r_x, r_y, r_z))

	new_coords = []
	for r in r_N:
		x_new = new_cms[0] + r[0]
		y_new = new_cms[1] + r[1]
		z_new = new_cms[2] + r[2]
		new_coords.append((x_new, y_new, z_new))
	return new_coords

def main():
	data = np.genfromtxt(fname = "coords.txt", dtype = float)

	molecule = []
	for z in zip(data[:,0], data[:,1], data[:,2]):
		molecule.append(z)

	plt.scatter(data[:,0], data[:,1])
	new_mol = convert_coords(molecule, [12,12,12,12,12,12,12,12,12,12], (0,0,0))
	for at in new_mol:
		plt.scatter(at[0], at[1], color="black")

	new_mol = convert_coords(molecule, [12,12,12,12,12,12,12,12,12,12], (125,125,0))
	for at in new_mol:
		plt.scatter(at[0], at[1], color="black")


	# plt.ylim(260, 350)
	plt.show()

if __name__ == "__main__":
	main() 

