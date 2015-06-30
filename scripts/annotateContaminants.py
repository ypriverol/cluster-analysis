__author__ = 'yperez'

import Bio
import sys, getopt
import csv
from contaminant import ClusterEntry
from utils.ClusterUtils import NA


def annotateContaminant(contaminants, clusterEntry):

    for contaminant in contaminants:

        if clusterEntry.pep_seq != None and clusterEntry.pep_seq in contaminant.seq:
            clusterEntry.pep_seq_contaminant.append(contaminant.id)
        if clusterEntry.pep_seq_second != None and clusterEntry.pep_seq_second in contaminant.seq:
            clusterEntry.pep_seq_second_contaminant.append(contaminant.id)
        if clusterEntry.pep_seq_third != None and clusterEntry.pep_seq_third in contaminant.seq:
            clusterEntry.pep_seq_third_contaminant.append(contaminant.id)
        if clusterEntry.pep_seq_four != None and clusterEntry.pep_seq_four in contaminant.seq:
            clusterEntry.pep_seq_four_contaminant.append(contaminant.id)

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
        contaminants.append(record)


    # Read the clusters information.
    clusters = []
    i = 0
    with open(input_file) as clusterFile:
        reader = csv.reader(clusterFile, delimiter='\t')
        for row in reader:
            if i != 0:
                # col = 0
                # for rowvalue in row:
                #     print(header[col], " \t", rowvalue)
                #     col += 1

                clusterId            = NA(row[0])                         # cluster id
                precursor_mz         = NA(row[1])                      # precursor_mz
                precursor_highest_mz = NA(row[2])
                avg_charge           = NA(row[3])
                avg_highest_charge   = NA(row[4])
                max_charge           = NA(row[5])
                min_charge           = row[6]
                max_highest_charge   = row[7]
                min_highest_charge   = row[8]
                max_mz               = row[9]
                min_mz               = row[10]
                mz_range             = row[11]
                max_mz_highest       = row[12]
                min_mz_highest       = row[13]
                mz_range_highest     = row[14]
                num_spectra          = row[15]
                num_projects         = row[16]
                project              = row[17]
                num_project_highest  = row[18]
                project_highest      = row[19]
                multi_highest        = row[20]
                num_peptides         = row[21]
                num_psms             = row[22]
                num_species          = row[23]
                species              = row[24]
                num_species_highest  = row[25]
                species_highest      = row[26]
                max_ratio            = row[27]
                pep_seq              = NA(row[28])
                pep_count            = NA(row[29])
                pep_seq_second       = NA(row[30])
                pep_count_second     = NA(row[31])
                pep_seq_third        = NA(row[32])
                pep_count_third      = NA(row[33])
                pep_seq_four         = NA(row[34])
                pep_count_four       = NA(row[35])
                filename             = NA(row[36])

                clusterEntry = ClusterEntry.ClusterEntry(clusterId,precursor_mz, precursor_highest_mz, avg_charge,
                                                         avg_highest_charge, max_charge, min_charge , max_highest_charge,
                                                         min_highest_charge, max_mz,  min_mz , mz_range, max_mz_highest, min_mz_highest,
                                                         mz_range_highest, num_spectra, num_projects,  project, num_project_highest,  project_highest, multi_highest,
                                                         num_peptides,  num_psms, num_species, species, num_species_highest,  species_highest, max_ratio, pep_seq,
                                                         pep_count,  pep_seq_second, pep_count_second,   pep_seq_third,  pep_count_third, pep_seq_four, pep_count_four, filename)
                annotateContaminant(contaminants, clusterEntry)
                clusters.append(clusterEntry)
            else:
                print(row)
                header = row
            i = i + 1

    print(len(clusters))
    f = open(output_file, 'w')
    for clusterEntry in clusters:
        f.write(clusterEntry.string)



if __name__ == "__main__":
    main(sys.argv[1:])


def usage():
    print('annotateContaminants.py -i <inputfile> -o <outputfile>')






