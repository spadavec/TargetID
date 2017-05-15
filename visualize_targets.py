import csv
import pandas
import operator
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import requests

def read_csv(filename):
    smiles = []
    chembl_id = []
    proba = []
    ec50 = []

    with open(filename) as ifile:
        data = csv.reader(ifile, delimiter=',')
        firstline = True

        for line in data:
            if firstline:
                firstline = False
                continue
            smiles.append(line[0])
            chembl_id.append(line[1])
            proba.append(line[2])
            ec50.append(line[3])


    return smiles, chembl_id, proba, ec50

def fetch_target_information(target):
    re = requests.get('https://www.ebi.ac.uk/chembl/api/data/target/{0}.json'.format(target))
    return (target, re.json()['organism'], re.json()['pref_name'])


def write_csv(output, filename):

    with open(filename, "w") as ofile:
        writer = csv.writer(ofile)
        writer.writerow(['CHEMBL_ID', 'Organism', 'Target Name', 'Target Count'])
        writer.writerows(output)


def main():
    output = []


    parser = argparse.ArgumentParser(description='options parser for visualize_targets.py')
    parser.add_argument('--input', dest="input_filename", required=True)
    parser.add_argument('--output', dest="output_filename", required=True)
    args = parser.parse_args()
    filename = args.input_filename
    output_filename = args.output_filename


    smiles, chembl_id, proba, ec50 = read_csv(filename)

    # Remove duplicates
    unique_chembl_id = set(chembl_id)

    print ''
    print "Total number of unique protein targets: {}".format(len(unique_chembl_id))


    # Generate dictionary for unique ids and their counts
    unique_chembl_id_d = {key: 0 for key in unique_chembl_id}
    for i,x in enumerate(chembl_id):
        unique_chembl_id_d[x] += 1

    sorted_chembl = sorted(unique_chembl_id_d.items(), key=operator.itemgetter(1), reverse=True)

    # Get the cutoff value so that only top 10 values are shown
    df = pandas.DataFrame(unique_chembl_id_d.items(), columns=['chembl_id', 'count'])
    dfList = df['count'].tolist()
    dfList.sort(reverse=True)
    cutoff = dfList[10]
    subset = df.loc[df['count'] > cutoff]

    subset_sorted = subset.sort_values(by='count', ascending=False)
    print ''
    print 'The top 5 most populated targets: '
    print subset_sorted.head(n=5)
    print ''

    unique_top_chembl_ids = subset_sorted['chembl_id'].tolist()


    print 'Grabbing information for unique targets...'
    for ids in unique_top_chembl_ids:
        target, organism, name = fetch_target_information(ids)
        #sprint "target, organism, name : {}, {}, {}".format(target, organism, name)
        temp = [ids, organism, name,unique_chembl_id_d[ids]]
        output.append(temp)

    print ''
    print "Writing output..."
    print ''

    write_csv(output, output_filename)

    print "Generating plot..."
    sns.set(style="white", context="talk")
    sns.barplot(subset_sorted['chembl_id'], subset['count'], palette="GnBu_d")
    plt.xticks(rotation=35)
    plt.tick_params(labelsize=10)
    sns.plt.title('Chembl ID Count')
    sns.plt.show()



if __name__ == '__main__':
    main()
