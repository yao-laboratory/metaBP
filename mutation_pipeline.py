import argparse
import os
import mutation_id
import numpy as np
import pandas as pd
import datetime

def pre_process(output_fp, read_file_1, read_file_2):
    bbmap_folder = "bbmap/bbmap"
    adapters = bbmap_folder+"/resources/adapters.fa"
    phiX_adapters=bbmap_folder+"/resources/phix174_ill.ref.fa.gz"
    filename1 = (read_file_1.split("/")[-1]).split(".")
    filename2 = (read_file_2.split("/")[-1]).split(".")
    temp1 = output_fp+'/'+filename1[0]+"temp1.fq"
    temp2 = output_fp+'/'+filename2[0]+"temp2.fq"

    output1 = output_fp+"/sequencing_reads/"+filename1[0]+"_filtered.fastq"
    output2= output_fp+"/sequencing_reads/"+filename2[0]+"_filtered.fastq"

    if os.path.exists(output_fp+"/sequencing_reads")==False:
       os.system("mkdir "+output_fp+"/sequencing_reads")

    if os.path.exists(output1)==False:
        command1 = "{}/bbduk.sh -Xmx7g in1={} in2={} out1={} out2={} minlen=10 qtrim=rl trimq=20 ktrim=r k=21 mink=11 ref={} hdist=1 threads=2 tbo tpe".format(bbmap_folder, read_file_1, read_file_2, temp1, temp2, adapters)
        os.system(command1)

        command2 = "{}/bbduk.sh in1={} in2={} out1={} out2={} ref={} hdist=1 k=21 threads=2".format(bbmap_folder, temp1, temp2, output1, output2, phiX_adapters)
        os.system(command2)
    
    if os.path.exists(temp1)==True:
       os.system("rm "+temp1+" "+temp2)
    
    return output1, output2

def isolate_clusters(cluster_fp):
    cluster_file = open(cluster_fp, "r+")
    file_data = []
    for i in cluster_file:
        file_data.append(i)
    cluster_file.close()
    clusters = {}
    for i in range(len(file_data)):
        cluster = []
        linedata = file_data[i].split(" ")
        hasChar = linedata[0].startswith(">")
        # find the identity sequence for a cluster
        if (len(linedata) == 1) and (hasChar == True):
            clusName = linedata[0]
            in_cluster = True
            counter = 0
            while in_cluster:
                i = i + 1
                next = file_data[i]
                next_split = next.split(" ")
                if i == (len(file_data) - 1):
                    cluster.append(next)
                    in_cluster = False
                    clusters[clusName] = cluster
                elif (len(next_split) == 1) and (next.startswith(">")):
                    in_cluster = False
                    i = i - 1
                    clusters[clusName] = cluster
                else:
                    cluster.append(next)
    return clusters

def isolate_short_prots(clusters_fp, output_fp):
    clusters = isolate_clusters(clusters_fp)
    fragments = {}

    # find clusters with short and long proteins
    mixed_keys = []

    for i in clusters.keys():
        cluster = clusters[i]
        for j in cluster:
            if len(j) > 100 and (j.startswith(">") == False):
                mixed_keys.append(i)

    long_prots = {}
    mixed_keys = list(set(mixed_keys))
    for mk in mixed_keys:
        cluster = clusters[mk]
        # remove the mixed cluster from the list of final clusters
        clusters.pop(mk)
        for j in cluster:
            ind = cluster.index(j)
            if (j.startswith(">") == False) and (len(j) > 100):
                long_prots[cluster[ind - 1]] = j
                cluster.pop(ind - 1)
                cluster.pop(ind - 1)
            elif (j.startswith(">") == False) and (len(j) < 100):
                fragments[cluster[ind - 1]] = j
                cluster.pop(ind - 1)
                cluster.pop(ind - 1)

    # now rewrite the clusters file, along with the fragments file and long protein file
    fragments_fp = output_fp+"/fragments.fasta"
    long_prots_fp = output_fp+"/long_proteins.fasta"
    short_prot_clusters_fp = output_fp+"/clusters_short_prots.fasta"
    short_prots_fp = output_fp+"/short_prots.fasta"
    fragment_file = open(fragments_fp, "w+")
    long_prots_file = open(long_prots_fp, "w+")
    cluster_file = open(short_prot_clusters_fp, "w+")
    short_prots = open(short_prots_fp, "w+")

    for i in fragments.keys():
        fragment_file.write(i)
        fragment_file.write(fragments[i])

    for i in long_prots.keys():
        long_prots_file.write(i)
        long_prots_file.write(long_prots[i])

    for i in clusters.keys():
        if len(clusters[i]) >= 2:
            cluster_file.write(i)
            for j in clusters[i]:
                cluster_file.write(j)
                short_prots.write(j)

    fragment_file.close()
    long_prots_file.close()
    cluster_file.close()
    short_prots.close()

    return short_prot_clusters_fp

def scrub_file_data(file_data):
    scrubbed_file_data = []
    for i in range(len(file_data) - 1):
        i_data = file_data[i].split(" ")
        i_next = file_data[i + 1].split(" ")
        if i == (len(file_data) - 1):
            k = i_data[0].replace("\n", "")
            scrubbed_file_data.append(k)
            scrubbed_file_data.append(file_data[i + 1])
            scrubbed_file_data.append(file_data[i + 2])
            scrubbed_file_data.append(file_data[i + 3])
        elif (len(i_data) == 1) and (len(i_next) == 1):
            k = file_data[i].replace("\n", "")
            scrubbed_file_data.append(k)
        else:
            scrubbed_file_data.append(
                file_data[i]
            )  # this works but the sequence is still split
    sequence_info = []
    consolidated_sequences = []
    sequence_info.append(scrubbed_file_data[0])
    for j in range(len(scrubbed_file_data) - 1):
        sequence = ""
        linedata = scrubbed_file_data[j].split(" ")
        hasChar = linedata[0].startswith(">")
        # find the identifier for a sequence
        if hasChar == True:
            in_sequence = True
            counter = 0
            while in_sequence:
                j = j + 1
                next = scrubbed_file_data[j]
                next_split = next.split(" ")
                if j == len(scrubbed_file_data) - 1:
                    in_sequence = False
                if next.startswith(">"):
                    in_sequence = False
                    sequence_info.append(scrubbed_file_data[j])
                    j = j - 1
                else:
                    sequence += next
            consolidated_sequences.append(sequence)
    sequences = []
    for i in range(len(consolidated_sequences)):
        sequences.append(sequence_info[i])
        sequences.append(consolidated_sequences[i])
    return sequences

def make_prot_df(output_fp, short_prots_fp):
    prot_file = open(short_prots_fp, "r+")

    file_data = []
    for i in prot_file:
        file_data.append(i)
    prot_file.close()
    proteins = {}
    for i in range(len(file_data)):
        if file_data[i].startswith(">"):
            info = file_data[i].split(" ")
            prot_id = info[0]
            proteins[prot_id] = file_data[i + 1]

    seqs = []
    for k in proteins.keys():
        seqs.append(proteins[k])

    acc = ["x"] * len(proteins)
    prot_df = pd.DataFrame()
    prot_df["accession"] = acc
    prot_df["identifier"] = proteins.keys()
    prot_df["protein"] = seqs

    csv_fp = output_fp+"tmp/sequences.csv"
    prot_df.to_csv(csv_fp)
    return csv_fp

def format_clusters(clusters_fp, output_fp):
    new_clus_fp = output_fp+"/clusters_all_seqs.fasta"
    new_clus = open(new_clus_fp, "w+")
    clus_df = pd.read_csv(clusters_fp)
    clus_df = clus_df[
        ["seqId", "identifier", "protein", "length", "repseq", "groupNum"]
    ]
    clusters = {}
    clus_df = clus_df.sort_values(by="groupNum", ascending=True)

    groups = clus_df["groupNum"].reset_index()
    rep_seqs = clus_df["repseq"].reset_index()
    rep_seqs = rep_seqs[rep_seqs["repseq"] == True]

    repseq = []
    for i in rep_seqs["index"]:
        repseq.append([clus_df["identifier"][i], clus_df["groupNum"][i]])

    for i in repseq:
        group = i[1]
        clusters[group] = clus_df[clus_df["groupNum"] == group].sort_values(
            by="repseq", ascending=False
        )

    for i in clusters.keys():
        cluster = clusters[i]
        rep = cluster[cluster["repseq"] == True]
        # change these to write later
        clus_str = rep["identifier"].item() + "\n"
        new_clus.write(clus_str)
        for index, row in cluster.iterrows():
            seq_info = "" + row["identifier"] + " len:" + str(row["length"]) + "\n"
            prot = row["protein"]
            new_clus.write(seq_info)
            new_clus.write(prot)

    new_clus.close()
    return new_clus_fp

def mutation_pipeline(
    read_file_1,
    read_file_2,
    seq_file_path,
    output_fp,
    clust_type,
    mutation_id_window=10,
):
    if os.path.exists(output_fp)==False:
        os.system("mkdir "+output_fp)
    if os.path.exists(output_fp+"/tmp")==False:
        os.system("mkdir "+output_fp+"/tmp")

    prots_fp = ""
    if (read_file_1 is not None) and (read_file_2 is not None):
        # filter data
        print(datetime.datetime.now(), ": Filtering sequencing read data")
        filtered_read1, filtered_read2 = pre_process(output_fp, read_file_1, read_file_2)
        print(datetime.datetime.now(), ": Filtering completed")
        # run plass
        plass_string = (
            "plass/bin/plass assemble "
            + filtered_read1
            + " "
            + filtered_read2
            + " " + output_fp + "/assembly.fas "+output_fp+"/tmp --min-seq-id 0.8 >/dev/null 2>&1"
        )
        print(datetime.datetime.now(), ": Running PLASS")
        print(plass_string)
        os.system(plass_string)
        print(datetime.datetime.now(), ": PLASS assembly completed")
        prots_fp = output_fp + "/assembly.fas"
    elif seq_file_path is not None:
        prots_fp = seq_file_path
    else:
        print("Error locating input files, please try again.")

    clusters_fp = ""
    if clust_type == 0:
        # run linclust
        linclust_str = "mmseqs/bin/mmseqs easy-linclust " +prots_fp+ " " +output_fp+"/clusters "+output_fp+"/tmp --min-seq-id 0.9 --cov-mode 1 -c 0.5 --cluster-mode 2 >/dev/null 2>&1"
        print(datetime.datetime.now(), ": Running Linclust")
        os.system(linclust_str)
        print(datetime.datetime.now(), ": Protein clustering complete")
        clusters_fp = output_fp+"/clusters_all_seqs.fasta"
    elif clust_type == 1:
        print(datetime.datetime.now(), ": Running RBiotools Linclust")
        # create csv file of protein sequences to pass in to RBiotools
        seq_csv = make_prot_df(output_fp, prots_fp)
        # then call R code
        os.system("Rscript get_clusters.r >/dev/null 2>&1")
        # identify mutations in clusters
        clusters_fp = format_clusters(output_fp+"/tmp/clusters.csv", output_fp)
        print(datetime.datetime.now(), ": Protein clustering complete")

    else:
        print("Error choosing correct version of Linclust, please try again.")
    
    # isolate clusters with short proteins
    print(datetime.datetime.now(), ": Isolating short proteins")
    clusters_short_prots = isolate_short_prots(clusters_fp, output_fp)
    print(datetime.datetime.now(), ": Short protein isolation complete")
    mutations_output_fp = output_fp+"/mutations.txt"
    print(datetime.datetime.now(), ": Identifying mutations")
    mutation_id.id_mutations(clusters_short_prots, mutations_output_fp, mutation_id_window, output_fp)
    print(datetime.datetime.now(), ": Mutation identification complete")
    # insert ML model steps here
    
def main():
    parser = argparse.ArgumentParser(
        prog="protein_based_annotation_pipeline",
        description="this method executes the mutation identification and annotation pipeline.",
    )
    subparser = parser.add_subparsers(
        dest="method",
        help="enter the the filepaths for the sequence read files and the name for the output file.",
    )

    mutation_id_parser = subparser.add_parser(
        "call_mutations",
        help="this method takes in sequence read data and returns mutation and annotation information for the sequences.",
    )
    mutation_id_parser.add_argument(
        "-i1",
        dest="read_file_1",
        type=str,
        help="the file path of the first sequence read file (.fastq.gz) to use in paired-end assembly",
    )
    mutation_id_parser.add_argument(
        "-i2",
        dest="read_file_2",
        type=str,
        help="the file path of the second sequence read file (.fastq.gz) to use in paired-end assembly",
    )
    mutation_id_parser.add_argument(
        "-s",
        dest="seq_file_path",
        type=str,
        help="the file path of the .fasta file of sequences if plass assembly is to be bypassed",
    )
    mutation_id_parser.add_argument(
        "-o", dest="output_file_path", type=str, help="the file path for the output directory"
    )
    mutation_id_parser.add_argument(
        "-clust", dest="clust_type", type=int, help="the version of linclust used for clustering. 0 for the original Linclust, 1 for the RBiotools Linclust"
    )
    
    mutation_id_parser.add_argument(
        "-mw",
        dest="mutation_id_window",
        type=int,
        help="the window on either side of a detected mutation within which to check that the similarity threshold (90%) is met. Default=10",
    )

    arguments = parser.parse_args()

    if arguments.method == "call_mutations":
        read_file_1 = arguments.read_file_1
        read_file_2 = arguments.read_file_2
        seq_file_path = arguments.seq_file_path
        output_file_path = arguments.output_file_path
        clust_type = arguments.clust_type
        mutation_id_window = arguments.mutation_id_window
        # call function
        mutation_pipeline(
            read_file_1,
            read_file_2,
            seq_file_path,
            output_file_path,
            clust_type,
            mutation_id_window,
        )
    else:
        print("Incorrect input, please check parameters and try again")
        
if __name__ == "__main__":
    main()
