from __future__ import print_function
import numpy as np
import os
import re
import subprocess
import shutil
import pandas as pd
import matplotlib.pyplot as plt


class US_Analysis:
	def __init__(self, paraent_directory):
		AT = []
		os.chdir(paraent_directory)
		folders = [name for name in os.listdir('./') if os.path.isdir(os.path.join('./', name))]

		for folder in folders:
			if bool(re.search(r"[-+]?\d*\.\d+|\d+", folder)):
				at = float(re.search(r"[-+]?\d*\.\d+|\d+", folder).group())
				if os.path.isfile('Z' + str(at) + '/colvar'):
					AT.append(at)

		AT.sort()
		self.min_at = AT[0]
		self.max_at = AT[-1]
		
		self.data = []
		for at in AT:
			for line in open('Z' + str(at) + '/plumed1.dat'):
				if re.search("KAPPA=", line):
					sl = line.split()
					for s in sl:
						if re.search("KAPPA", s):
							spring_constant = float(s.split("=")[1])

			df = pd.read_csv('Z' + str(at) + '/colvar', delim_whitespace=True, header=None, usecols=[0,1], comment='#')
			self.data.append({"at": at, "spring_constant": spring_constant, "colvar": df.values})



		# ### CONSIDER PBC
		# for each in self.data:
		# 	colvar = each["colvar"]

		# 	for i in range(len(colvar)):
		# 		x = colvar[i][0]
		# 		y = colvar[i][1]
		# 		yf = each["at"]

		# 		ysize = 1
		# 		dy = y - yf
		# 		dy -= ysize * int(round(dy * ysize))
		# 		y = yf + dy

		# 		colvar[i][1] = y



	def __repr__(self):
		ats = []
		scs = []
		
		for each in self.data:
			ats.append(str(each["at"]))
			scs.append(str(each["spring_constant"]))

		at = ", ".join(ats)
		sc = ", ".join(scs)
		nwins = len(self.data)

		return "NUMBER OF WINDOWS: {}\n AT: {}\n SPRING_CONSTANTS: {}".format(nwins, at, sc)


	# def _get_histogram(self, data, min_value, max_value, nbins):
	# 	x = np.linspace(min_value, max_value, nbins)
	# 	y = np.zeros(nbins)

	# 	for value in data:
	# 		i = int((value - min_value)*nbins/(max_value-min_value))
	# 		y[i] += 1
	# 	y /= len(data)

	# 	return np.stack((x,y), axis = -1)


	def comp_prob(self, min_value, max_value, nbins):
		bins = np.linspace(min_value, max_value, nbins+1)
		for each in self.data:
			colvar = each["colvar"]
			read_start = int(len(colvar) * 0.2)
			read = colvar[read_start:]
			y, x = np.histogram(read[:,1], bins=bins)
			y = y/len(read)
			x = x[1:]
			each["prob"] = np.stack((x,y), axis=-1)
			# each["prob"][:,0] vs each["prob"][:,1]
		self.save_prob()


	def save_prob(self):
		for each in self.data:
			x = each["prob"][:,0]
			y = each["prob"][:,1]
			plt.plot(x, y, label=str(each["at"]))
		plt.legend()
		plt.savefig('prob.png')

	def plot_prob(self):
		for each in self.data:
			x = each["prob"][:,0]
			y = each["prob"][:,1]
			plt.plot(x, y, label=str(each["at"]))
		plt.legend()
		plt.show()

	def comp_re(self, min_value, max_value, nbins, nblocks=5, cut=0.08):
		bins = np.linspace(min_value, max_value, nbins+1)
		
		for each in self.data:
			colvar = each["colvar"]
			g = np.histogram(colvar[:,1], bins=bins)[0] / len(colvar)

			re = []

			for n in range(nblocks):
				start = n/nblocks
				end = (n+1)/nblocks
				read_start = int(len(colvar) * start)
				read_end = int(len(colvar) * end)
				read = colvar[read_start:read_end]

				f = np.histogram(read[:,1], bins = bins)[0] / len(read)

				f_g = np.zeros(nbins)
				for i in range(nbins):
					if g[i] == 0:
						f_g[i] = 1
					elif f[i] == 0:
						#f[i] = 1/len(read)
						f_g[i] = 1/len(read)/g[i]
					else:
						f_g[i] = f[i]/g[i]

				re.append(f @ np.log(f_g))

			each["re"] = re


		for each in self.data:
			for val in each["re"]:
				if val > cut:
					p = ""
					for val in each["re"]:
						p += "{: 2.3f}  ".format(val)
					print("AT {: 2.3f}: {}".format(each["at"], p))
					break


	def make_metadata(self, prefix='./'):
		file = open(prefix + 'metadata.dat', 'w')
		for each in self.data:
			file.write('{} {} {}\n'.format('Z' + str(each["at"]) + "/colvar", each["at"], each["spring_constant"]))
		file.close()


	def errors(self, nblocks=5):
		shutil.rmtree('errors', ignore_errors=True)
		os.makedirs('errors')

		for each in self.data:
			colvar = each["colvar"]
		
			for n in range(nblocks):
				start = n/nblocks
				end = (n+1)/nblocks
				read_start = int(len(colvar) * start)
				read_end = int(len(colvar) * end)
				read = colvar[read_start:read_end]

				file = open('errors/col{}_{}'.format(each["at"], n), 'w')
				
				for values in read:
					x  = values[0]
					y  = values[1]
					
					yf = each["at"]
					ysize = 1
					dy = y - yf
					dy -= ysize * int(round(dy * ysize))
					y = yf + dy

					file.write('{: 2.4f} {: 2.4f}\n'.format(x, y))
				file.close()


		for n in range(nblocks):
			file = open('errors/metadata{}.dat'.format(n), 'w')
			for each in self.data:
				file.write('col{}_{} {} {}\n'.format(each["at"], n, each["at"], each["spring_constant"]))
			file.close()


		file = open('errors/run.sh', 'w')
		for n in range(nblocks):
			file.write('wham {} {} 300 1e-10 310 0 metadata{}.dat pmf{}.dat\n'.format(self.min_at, self.max_at, n, n))
		file.close()



	def add(self, AT, spring_constant):
		file = open('add_simulation.sh', 'w')
		file.write('''#!/bin/bash
cp -r template Z{}
cp eq/step7.gro Z{}
cd Z{}
sed -i.bak 's/_AT/{}/' run_md1.sub
sed -i.bak 's/_AT/{}/' run_md2.sub
sed -i.bak 's/_KAPPA/{}/' plumed1.dat
sed -i.bak 's/_KAPPA/{}/' plumed2.dat
sed -i.bak 's/_AT/{}/' plumed1.dat
sed -i.bak 's/_AT/{}/' plumed2.dat
rm -rf *.bak
cd ..
'''.format(AT, AT, AT, AT, AT, spring_constant, spring_constant, AT, AT))
		file.close()

		subprocess.call(['bash', 'add_simulation.sh'])
		subprocess.call(['rm', '-rf', 'add_simulation.sh'])
		








