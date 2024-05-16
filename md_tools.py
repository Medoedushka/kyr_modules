import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import re
import sys

def numdensity3D(boxsize : tuple, particles : np.ndarray, h : int = 5, axis : int = 2) -> np.ndarray:
	"""
	Calculates a number density of particles inside a box with specified dimensions.

	Parameters
	----------
	boxsize : tuple
		A 3-tuple with components x, y, z that are specifing dimensions of a box.
	particles : ndarray
		NumPy ndarray of the shape (n, 3) with 3D coordinates of the n particles.
	h : int, default=5
		Specifies how datailed the number density function should be calculated.
	axis : {0, 1, 2}
		An axis along which the number density function will be calculated.

	Returns
	-------
	rho : ndarray
		A number denisty function with the shape (n, 2). Values at (n, 0) are corresponding to coordinates along the `axis`;
		values at (n, 1) corresponding to the number density function.
	"""
	
	# Validation of input parameters
	if (boxsize[0] == 0 or boxsize[1] == 0 or boxsize[2] == 0):
		raise ValueError("Calculation of 3D number density function requieres a three dimensional box! \
				   One of the submitted components is equal zero or missing.")
	if (h <= 0):
		raise ValueError("Wrong value of calculation step!")
	if (not (axis in {0, 1, 2})):
		raise ValueError("Improper axis value was submitted!")

	move = [0, 0, 0] # a vector that transform coordinatates of particles to internal coordinate system
	move[axis] = min(particles[:,axis])
	q_max = boxsize[axis] # maximum value of the internal coordinate q along which particles will be counted

	# search for axes that will be constant
	free_axis = [ax for ax in [0, 1, 2] if ax != axis]
	L1 = boxsize[free_axis[0]]
	L2 = boxsize[free_axis[1]]
	transformed_cms = particles - move 

	rho = np.ndarray(shape=(1,2))
	q = 0
	counter = 0
	while (q <= q_max):
		for cms in transformed_cms:
			if (cms[axis] >= q and cms[axis] <= (q + h)): counter += 1
		rho_q = counter / (L1 * L2 * h)
		rho = np.append(rho, [[q, rho_q]], axis=0)
		q += h
		counter = 0
	rho = np.delete(rho, 0, axis=0)
	return rho

def gmx_extract(index, path, rename):
    """
    Extracts data from GROMACS .edr-file and renames it.
    Returns averaged value of an extracted property.
    
    index   :      index of data file in the GROMACS gmx_mpi energy module
    path    :      a path to a .edr-file
    rename  :      new name for extracted data-file
    """

    
    COMMAND = f"echo {index} 0 | gmx_mpi energy -f {path} &> /dev/null" # &> /dev/null -> hides output of a command
    subprocess.run(COMMAND, shell=True, executable="/bin/bash") # TRY CATCH
    os.rename("energy.xvg", rename + ".xvg")

    data = np.genfromtxt(rename + ".xvg", skip_header=24, dtype=float)
    return np.average(data[:,1])

def get_cms(mol_coord, weights):
	mol_mass = np.sum(weights)

	x_sum = 0	#
	y_sum = 0	# weighted sum of each component
	z_sum = 0	#

	for mol in zip(mol_coord, weights):
		x_sum += mol[0][0] * mol[1]
		y_sum += mol[0][1] * mol[1]
		z_sum += mol[0][2] * mol[1]

	x_cms = x_sum / mol_mass
	y_cms = y_sum / mol_mass
	z_cms = z_sum / mol_mass
	return x_cms, y_cms, z_cms

def get_lastsim(fname : str):
	"""
	Searches in a folder <fname> a simulation part with a biggest number.

	fname : a name of a folder which has simulations sim1, sim2, etc.
	"""

	files = os.listdir(fname)
	r1 = re.compile(r'sim\d+\.gro')
	indexes = []
	for f in files:
		if (r1.match(f)):
			r2 = re.compile(r'[^\d+]')
			indexes.append(int(r2.sub('', f)))
	return max(indexes)
	
def prep_orthoboxy_sim(Dx, sx):
    sw = 0.31
    Dw = 2.3e-9
    tau_w = 0.5
    f = 2500
    tau_x = Dw/Dx * (sx/sw)**2 * tau_w # in ns
    dtau_x = tau_x / f * 1000 # in ps
    n = tau_x / 2e-6
    dn = dtau_x / 2e-3
    
    print(f"tau_x = {tau_x} ns")
    print(f"dtau_x = {dtau_x} ps")
    print(f"n = {n}")
    print(f"dn = {dn}")

# -------- TESTING --------
def main():
	print(sys.path)
if __name__ == "__main__":
	main()
