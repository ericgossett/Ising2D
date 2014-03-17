import matplotlib.pyplot as plt
from scipy import optimize
from numpy import *
import numpy as np
import csv 
k_b = float(1.3806503*pow(10,-23))
#==========================================================================================
def columnExtract(c1,filename):
	outfile = open(filename,'rU')
	data = [row[c1] for row in csv.DictReader(outfile)]
	outfile.close()
	return data

#==========================================================================================
def IsingPlot():
	
	fig1 = plt.figure(1)

	#`````````````````````````````SUBPLOT 1`````````````````````````````#
	subplot1 = fig1.add_subplot(211)
	filepath ="/Users/Eric/Github_Repos/Ising2D/Data.csv"
	
	x = columnExtract('T',filepath)
	x = asarray(x,dtype=float)
	y = columnExtract('M',filepath)
	y = asarray(y,dtype=float)
	y_abs = np.absolute(y)

	subplot1.plot(x,y_abs,'go')
	plt.grid(True)
	plt.axvline(x=2.269,color='r',ls='dashed')
	plt.ylabel(' $\mid <M> \mid\ per\ site\ \mu]$')
	plt.title('$Temperature\ vs.\ \mid <M> \mid\ per\ site$')	


	subplot3 = fig1.add_subplot(212)

	x = columnExtract('T',filepath)
	x = asarray(x,dtype=float)
	y = columnExtract('E',filepath)
	y = asarray(y,dtype=float)
	
	subplot3.plot(x,y,'yo')
	plt.grid(True)
	plt.axvline(x=2.269,color='r',ls='dashed')
	plt.xlabel('$Temperature\ [J/k_{B}]$')
	plt.ylabel('$< E >\ per\ site\ [ J ]$')
	plt.title('$Temperature\ vs.\ < E >\ per\ site$ ')
	
	

	fig3 = plt.figure(2)
	subplot4 = fig3.add_subplot(211)

	x = columnExtract('T',filepath)
	x = asarray(x,dtype=float)
	y = columnExtract('C',filepath)
	y = asarray(y,dtype=float)
	
	
	subplot4.plot(x,y,'ko')
	
	plt.grid(True)
	plt.axvline(x=2.269,color='r',ls='dashed')
	plt.ylabel('$C_{v}\ [J/K_{B}$')
	plt.title('$Specfic\ Heat\ C_{v}\ vs.\ Temperature $')

	subplot4= fig3.add_subplot(212)

	T = columnExtract('T',filepath)
	T = asarray(T,dtype=float)
	X = columnExtract('X',filepath)
	X = asarray(X,dtype=float)
	X_abs = np.absolute(X)
	
	subplot4.plot(T,X_abs,'bo')
	plt.grid(True)
	plt.axvline(x=2.269,color='r',ls='dashed')
	plt.xlabel('$Temperature\ [J/k_{B}]$')
	plt.ylabel('$\chi\ [\mu/k_{B}]$')
	plt.title('$Magnetic\ Suseptibility\ \chi\ vs.\ Temperature$')

	
	plt.show()

#==========================================================================================
def main():
	IsingPlot()

if __name__ == "__main__":
	main()