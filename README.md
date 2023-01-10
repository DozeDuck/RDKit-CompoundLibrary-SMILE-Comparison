# RDKit-CompoundLibrary-SMILE-Comparison
# 用户提供分子SMILE库（如本文件）+ 目标化合成的SMILE（如本文件）
$ python rdkit_fingerprint.py -t gsk2830371.smi -d smile_test.dat -e error_list.dat -o rank_list.dat -n 100
-t: target compound's SMILE
-d: data set SMILE library
-e: the SMILEs which can't be used in the data set SMILE library
-o: output file, ranked based on Tanimoto similarity
-n: amount of SMILEs for output 
