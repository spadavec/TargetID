from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.externals import joblib
import csv
import pandas as pd
import numpy as np
import sys
from collections import OrderedDict
import requests
from rdkit import rdBase
import argparse



def read_csv(filename):
    smiles = []
    activities = []

    with open(filename) as ifile:
        data = csv.reader(ifile, delimiter=',')

        for line in data:
            smiles.append(line[0])
            activities.append(line[1])

    return smiles, activities

def write_csv(output, filename):
    outfile = filename[:-4] + '-out.csv'

    with open(outfile, "w") as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['SMILES', 'CHEMBL ID', 'Probability', 'EC50'])
        writer.writerows(output)

def main():
    collect = []

    parser = argparse.ArgumentParser(description='options parser for TargetID')
    parser.add_argument('--input', dest="input_filename", required=True)
    parser.add_argument('--threshold', dest="threshold", required=True)
    parser.add_argument('--output', dest="output_filename", required=True)
    args = parser.parse_args()

    filename = args.input_filename
    threshold = args.threshold
    output_filename = args.output_filename

    smiles, activities = read_csv(filename)

    morgan_nb = joblib.load('models_21/1uM/mNB_1uM_all.pkl')
    classes = list(morgan_nb.targets)

    print "Getting top predictions..."
    for i,x in enumerate(smiles):
        mol = Chem.MolFromSmiles(x)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=2048)
        res = np.zeros(len(fp), np.int32)
        DataStructs.ConvertToNumpyArray(fp,res)

        probas = list(morgan_nb.predict_proba(res.reshape(1,-1))[0])
        predictions = pd.DataFrame(list(zip(classes, probas)), columns=[ 'id', 'proba'])

        top10_pred = predictions.sort_values(by = 'proba', ascending = False).head(50)

        for index, row in top10_pred.iterrows():
            if float(row['proba']) >= float(threshold):
                result = [x, row['id'],row['proba'], activities[i]]
                collect.append(result)


    write_csv(collect,output_filename)


if __name__ == '__main__':
    main()
