"""
Author: DozeDuck-Shangze
Date: 2023-Jan-10th
"""
import sys
import os
import getopt
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


args=sys.argv[1:]                                                                   # get arguments from usr's input; filename is sys.arg[0], so args start from [1:]
smile_data_lib = ''
target_molecule_smile = ''
rank_list_name = 'rank_list.dat'
# similarity_threshold = ''
number_output = 100
error_list = 'error_list.dat'


#optarg for the input example-c 2 -i PLDXPAL -t done.trg -p 50 -N 6 -n 80000
try:
   opts, args = getopt.getopt(args,"h:t:d:o:s:n:e:",["help","target_molecule_smile =",
                                    "smile_data_lib =", "rank_list_name =", "similarity_threshold =", "number_output =", "error_list ="])
except getopt.GetoptError:
   print ('rdkit_fingerprint.py -t <target_molecule_smile.smi> -d <smile_data_lib.dat> -o <rank_list_name.dat> -s <0.8> -n <1000> -e <error_smi.dat>')
   sys.exit(2)                                                                      # Exiting the program raises a SystemExit exception, 0 means normal exit, and others are abnormal exits.

for opt, arg in opts:                                                               # Generate several pairs of value, e.g: opr,arg = -i,PLDXPAL
   if opt == '-h':
      print ('rdkit_fingerprint.py -t <target_molecule_smile.smi> -d <smile_data_lib.dat> -o <rank_list_name.dat> -s <0.8> -n <1000> -e <error_smi.dat>')
      sys.exit()
   elif opt in ("-t", "--target_molecule_smile"):
      target_molecule_smile = str(arg)
   elif opt in ("-d", "--smile_data_set"):
      smile_data_lib = str(arg)
   elif opt in ("-o", "--rank_list_name"):
      rank_list_name = str(arg)
   elif opt in ("-s", "--similarity_threshold"):
      similarity_threshold = float(arg)
   elif opt in ("-n", "--number_output"):
      number_output = int(arg)
   elif opt in ("-e", "--error_list"):
      error_list = str(arg)


class similarity_calculator():
	smile_lib = []							# build a list contains all smile strings
	target_smile = "lalala"

	def __init__(self, smile_data_lib, target_smile_file, error_list, rank_list_name, number_output):
		self.target_smile_reader(target_smile_file)
		self.dataset_reader(smile_data_lib)
		self.build_dict(self.smile_lib, self.target_smile, error_list, rank_list_name, number_output)


	def target_smile_reader(self,filename):
		f = open(filename, "r")
		for line in f:
			self.target_smile = str(line)
		f.close()

	def dataset_reader(self,filename):
		f = open(filename, "r")
		data = f.readlines()                                                     # To make each line as a list so that can be cutted into slices
		for line in data:
			self.smile_lib.append(line[0:-1])                                    # To delete '\n'
		f.close()

	def tanimoto_calc(self,smi1, smi2):
		mol1 = Chem.MolFromSmiles(smi1)
		mol2 = Chem.MolFromSmiles(smi2)
		fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
		fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
		s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
		return s

	def build_dict(self, smile_lib, target_smile, error_list, rank_list_name, number_output):
		os.system('rm ' + error_list)
		# print("the dataset is: " + smile_lib[0])
		# print("the target smi is: " + target_smile)
        # similarity_dict ={}  # 声明字典
		similarity_dict = {smile_lib[0]:self.tanimoto_calc(str(target_smile), str(smile_lib[0]))}
		# print("Cooooooooooooooooll")
		# print("1575 element of dict: " + str(similarity_dict))
		a=1
		for i in range(len(smile_lib)):
			# print(a)
			try:
			    similarity_dict[smile_lib[i]] = self.tanimoto_calc(str(target_smile), str(smile_lib[i]))
			except:
			    os.system('echo "' +str(a)+ ' error smiles string: ' +  str(smile_lib[i]) + '" >> ' + error_list)
			a+=1
		# print(similarity_dict)
		with open(rank_list_name,"w") as rank_list:
			a = 1
			for i in sorted(similarity_dict.items(), key = lambda kv:(kv[1], kv[0]), reverse=True):
			    if(a <= number_output):
			        js = json.dumps(i, ensure_ascii=False)
			        rank_list.write(js)
			        rank_list.write('\n')
			    a += 1
            

a = similarity_calculator(smile_data_lib, target_molecule_smile, error_list, rank_list_name, number_output)

