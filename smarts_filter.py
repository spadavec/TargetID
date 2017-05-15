import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import argparse




def read_csv(filename):
    smiles = []
    activity = []

    with open(filename) as ifile:
        data = csv.reader(ifile, delimiter=',')

        for line in data:
            smiles.append(line[0])
            activity.append(line[1])

    return smiles, activity

def retain_mcs(smarts, smiles, activities):
    """ Finds all compounds in input file which contain designated SMARTS pattern """

    retained = []
    core_mols = [Chem.MolFromSmarts(x) for x in smarts]
    print ''
    print "Finding SMARTS matches..."
    for i,x in enumerate(smiles):
        mol = Chem.MolFromSmiles(x)
        for smarts in core_mols:
            if mol.HasSubstructMatch(smarts):
                retained.append([x, activities[i]])

    return retained

def write_csv(output, filename):
    with open(filename, "wb") as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['SMILES', 'Activities'])
        writer.writerows(output)


def main():
    parser = argparse.ArgumentParser(description='options parser for smarts_filter.py')
    parser.add_argument('--input', dest="input_filename", required=True)
    parser.add_argument('--output', dest="output_filename", required=True)
    parser.add_argument('--smarts', dest="smarts_pattern")
    args = parser.parse_args()
    filename = args.input_filename
    output_filename = args.output_filename
    smarts_pattern = args.smarts_pattern

    smiles, activities = read_csv(filename)
    if smarts_pattern:
        match = smarts_pattern
    else:
        match = ['[#6]-[#8]-[#6]1:[#6]:[#7]:[#6]:[#6]2:[#7]:1:[#6](:[#7]:[#7]:2)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1','[#7]-[#6](=O)-c1cncc2nnc(-c3ccccc3)n12']

    print ''
    print "Total of {} entry molecules ".format(len(smiles))
    matched = retain_mcs(match, smiles, activities)
    print ''
    print "Total of {} matching compounds ".format(len(matched))
    print ''

    write_csv(matched, output_filename)


if __name__ == '__main__':
    main()
