#!/bin/env python
import subprocess
import sys

def install(packages):
	subprocess.call([sys.executable, "-m", "pip", "install", package])

if __name__=='__main__':
	required_packages=['tensorly', 'argparse', 'joblib', 'matplotlib', 'lifelines', 'seaborn']
	strout=[]
	valid=1
	for package in required_packages:
		try:
			install(package)
		except:
			strout.append(package)
			valid=0

	if valid:
		print('\nAll required packages installed.')
		print('Ready to execute monti.py')
	else:
		print('ERROR')
		for e in strout:
			print("%s module not correctly installed."%(e))

	


