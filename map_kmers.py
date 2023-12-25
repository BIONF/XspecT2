import statistics
import pickle
import os
from sklearn.cluster import KMeans
import numpy as np


def create_genome_kmer_list(genome_dir, k, genus):
    """Erstellt für jedes Genom eine Liste aller kmere"""

    # list of files in the directory
    files = sorted(os.listdir(genome_dir))

    # loop over all files in the directory
    for file in files:
        # path to the current genome
        genome_path = genome_dir + "/" + file

        # Create a list of all kmers in a genome/species
        kmer_list = []

        # open the file with the genome
        with open(genome_path, "r") as file:
            genome = file.readlines()

        # loop over all contigs in the genome
        for contig in genome:
            # loop over all kmers in the contig
            for i in range(len(contig) - k + 1):
                kmer = contig[i : i + k]
                kmer_list.append(kmer)

        # filename for the pickle file
        filename = genome_path.split("/")[-1] + "_kmer_list.txt"

        # saves the list of all kmers in the genome in a pickle file
        with open("kmer_positions/" + genus + "/" + filename, "wb") as file:
            pickle.dump(kmer_list, file)


def identify_split_reads(score, kmer_hits_single):
    """This method identifies if the read is split because of HGT"""

    # notes
    """ 
    Identify split read regions from the kmer profiles, one side should cluster Ones and the other side should cluster Zeros
    for the respective  species.
    1.) Collect the kmer profiles, transforms the kmer_hit vector into a species_hit vector --> shows which positions in the read are covered
    2.) Identify part of kmer profiles with a high number of ones --> shows where the read is covered
    3.) Divide the kmer profiles into two clusters --> One with Zeros and one with Ones
        --> Maximize the number of ones in the first cluster and the number of zeros in the second cluster
    """

    # initialize variables
    kmer_profiles = [] * len(score)
    split_regions = []
    index_result = max(range(len(score)), key=score.__getitem__)

    # Collect the kmer profiles, transforms the kmer_hit vector into a species_hit vector
    for kmer_hits in kmer_hits_single:
        for i in range(len(kmer_hits)):
            kmer_profiles[i].append(kmer_hits[i])

    og_species = None
    novel_species = None
    clusters = cluster_zeros_ones(kmer_profiles[index_result])
    split_regions.append((clusters, i))
    # which part of the read belongs to the original species, 0 =  left, 1 = right
    if sum(clusters[0]) > sum(clusters[1]):
        og_species = 0
        novel_species = 1
    else:
        og_species = 1
        novel_species = 0

    for i, kmer_profile in enumerate(kmer_profiles):
        if i == index_result:
            continue
        clusters = cluster_zeros_ones(kmer_profile)
        # 0.3 and 0.6 are arbitrary values, they should be tested etc.
        # Find the complemantary species of the split read --> HGT Donor
        if (
            sum(clusters[og_species]) / len(clusters[og_species]) < 0.3
            and sum(clusters[novel_species]) / len(clusters[novel_species]) >= 0.6
        ):
            split_regions.append((clusters, i))
            break

    # split_regions = [([cluster_0, cluster_1], og_index), ([cluster_0, cluster_1], novel_index)]
    return split_regions


def cluster_zeros_ones(input_list, threshold=None):
    """This method divides a list of zeros and ones into two lists,
    maximizing occurences of zeros and ones in each list"""

    # min. length of a cluster
    if threshold == None:
        threshold = len(input_list) * 0.1

    # copy the input list
    input_list_copy = input_list[:]

    # convert zeros to -1 so that they function as penalty
    input_list[:] = [-1 if x == 0 else 1 for x in input_list]
    cluster_score = 0
    # calculate the score for the cluster for each possible split
    for i in range(threshold, len(input_list) - threshold):
        # goal is to maximize the score --> highest score is the best split, contains the most ones and least zeros (-1)
        score = max(sum(input_list[:i]), sum(input_list[i:]))
        if score > cluster_score:
            cluster_score = score
            split_index = i

    # split the input list into two clusters
    cluster_0 = input_list_copy[:split_index]
    cluster_1 = input_list_copy[split_index:]

    return [cluster_0, cluster_1]


def map_kmers(kmers, predictions, genus):
    """Map kmers to their respective genome"""

    # check for pure- or split-read
    if len(predictions) == 1:
        prediction = predictions[0]

        # get the kmer list for the genome
        # TODO: decide filename
        filename = "Test"
        with open(
            "kmer_positions/" + genus + "/" + filename + "_kmer_list.txt", "rb"
        ) as fp:
            kmer_list = pickle.load(fp)

        # map kmers to the genome
        kmer_positions = []
        for kmer in kmers:
            # get all indices of the kmer in the kmer list
            indices = [i for i, x in enumerate(kmer_list) if x == kmer]
            kmer_positions += indices

        # sort the list of kmer positions in ascending order
        # is this biologically correct? Its neccecary for not unqiue kmers but the true position of the kmer is lost
        kmer_positions.sort()

        # 1. Cluster iteration: cluster kmers that lie directly next to each other
        clusters = []
        for i, position in enumerate(kmer_positions):
            # skip already clustered kmers
            if position in cluster:
                continue
            cluster = []
            cluster.append(position)
            # add all kmers that lie directly next to each other to the cluster
            for j in range(i + 1, len(kmer_positions)):
                if kmer_positions[j] == kmer_positions[j - 1] + 1:
                    cluster.append(kmer_positions[j])
                else:
                    clusters.append(cluster)
                    break

        # 2. Cluster iteration: cluster clusters that are at max 10 kmers apart
        final_clusters = []
        # threshold for the distance between two clusters
        distance_threshold = 5
        for i, cluster in enumerate(clusters):
            if cluster[i] != final_clusters[-1]:
                final_clusters += cluster[i]
            if (clusters[i + 1][0] - clusters[i][-1]) < distance_threshold:
                final_clusters += clusters[i + 1]


# TODO:
# rename function to map_kmers and split into two functions -> second one to cluster kmers
def cluster_kmers(kmer_list, kmer_dict):
    """Map kmers to their respective genome"""
    clusters = {}
    contig_median_list = []
    # Schleife über alle kmere in kmer_list
    for i in range(len(kmer_list)):
        kmer = kmer_list[i]
        # Überprüfen, ob das kmer in kmer_dict vorhanden ist
        kmer_info = kmer_dict.get(kmer)
        if kmer_info is None:
            continue
        # Holt die Contig-ID und kmer-Position aus kmer_dict
        kmer_id = kmer_dict[kmer][1]
        kmer_pos = kmer_dict[kmer][0]
        # Füge das kmer dem entsprechenden Contig in das Dictionary hinzu
        # Ein Contig ist ein cluster
        if kmer_id not in clusters:
            clusters[kmer_id] = []
        clusters[kmer_id].append((kmer, kmer_pos))
    # Schleife über alle Contigs im Dictionary clusters
    for contig in clusters:
        contig_list = clusters[contig]
        contig_len = len(contig_list)
        if contig_len < 2:
            # print("Zu wenig kmere im Contig!")
            continue
        # Sortieren der kmere in der Liste nach der Position
        sorted_contig_list = sorted(contig_list, key=lambda x: x[1])
        distances = []
        # Schleife über die sortierte Liste von kmere im Contig
        for i in range(1, len(sorted_contig_list)):
            kmer_pos = sorted_contig_list[i][1]
            prev_kmer_pos = sorted_contig_list[i - 1][1]
            # Berechnung der Distanz zum vorherigen kmer
            distance = kmer_pos - prev_kmer_pos
            distances.append(distance)
        # Berechnung der Median-Entfernung für das aktuelle Contig
        # print(distances)
        median_distance = statistics.median(distances)
        # print(median_distance)
        contig_median_list.append((contig, median_distance, contig_len))
    # Summe der Contigs in contig_median_list
    num_contigs = len(contig_median_list)
    # Liste aller Contig-Größen
    contig_lengths = [x[2] for x in contig_median_list]
    # Liste aller Median-Entfernungen
    median_distances = [x[1] for x in contig_median_list]
    # Median der Median-Entfernungen berechnen
    if len(median_distances) > 0:
        median_of_medians = statistics.median(median_distances)
    else:
        median_of_medians = None
    # Ergebnisliste erstellen
    result = [num_contigs, median_of_medians, contig_lengths]
    return result


def main():
    # Verzeichnis mit den Genomen
    genome_dir = "path/to/genomes"
    # create_genome_kmer_list(genome_dir, 21, "Acinetobacter")


if __name__ == "__main__":
    main()
