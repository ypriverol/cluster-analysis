__author__ = 'yperez'

import Bio
import sys, getopt
import csv
from contaminant import ClusterEntry


def annotateContaminant(contaminants, clusterEntry):

    for contaminant in contaminants:

        if clusterEntry.max_sequence in contaminant.seq:
            clusterEntry.max_sequence_contaminant.append(contaminant.id)
        if clusterEntry.second_max_sequence in contaminant.seq:
            clusterEntry.second_max_sequence_contaminant.append(contaminant.id)
        if clusterEntry.third_max_sequence in contaminant.seq:
            clusterEntry.third_max_sequence_contaminant.append(contaminant.id)

def main(argv):
    input_file = ''
    output_file = ''
    contaminant = ''

    # Get the scripts parameters
    try:
        opts, args = getopt.getopt(argv, "hi:o:c:", ["ifile=", "ofile=", "cfile"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
        elif opt in ("-c", "--cfile"):
            contaminant = arg

    # Read the contaminant database

    contaminants = []
    from Bio import SeqIO

    handle = open(contaminant, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id)
        contaminants.append(record)

    print(len(contaminants))


    # Read the clusters information.
    clusters = []
    i = 0
    with open(input_file) as clusterFile:
        reader = csv.reader(clusterFile, delimiter='\t')
        for row in reader:
            if i != 0:
                clusterEntry = ClusterEntry.ClusterEntry(row[0], row[1],
                                                         row[2], row[3],
                                                         row[4], row[5],
                                                         row[6], row[7],
                                                         row[8], row[9],
                                                         row[10], row[11],
                                                         row[12], row[13],
                                                         row[14], row[15],
                                                         row[16], row[17],
                                                         row[18], row[19],
                                                         row[20])
                annotateContaminant(contaminants, clusterEntry)
                clusters.append(clusterEntry)
            else:
                print(row)
            i = i + 1

    print(len(clusters))


if __name__ == "__main__":
    main(sys.argv[1:])


def usage():
    print('annotateContaminants.py -i <inputfile> -o <outputfile>')






