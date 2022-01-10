import os
import mutation_id
import argparse
import pandas as pd
import numpy as np

def isolate_short_prots(prots_fp):
    # function to isolate proteins that are from 1-100 amino acids in line (peptides)
    sequences_file = open(prots_fp, "r+")
    file_data = []
    
    for i in sequences_file:
        file_data.append(i)
        
    # consolidate sequences on multiple lines
    sequences = scrub_file_data(file_data)

    len_lim = 101
    short_sequences = []
    
    for i in range(len(sequences)):
        check_start = sequences[i].startswith(">")
        seq_len = len(sequences[i])
        if (check_start == False and seq_len < len_lim):
            short_sequences.append(sequences[i-1])
            short_sequences.append(sequences[i])
    
    short_seq_fp = "short_assembly.fas"
    short_seq_file = open(short_seq_fp, "w+")
    
    for i in short_sequences:
        short_seq_file.write(i)

    return short_seq_fp

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
            scrubbed_file_data.append(file_data[i])  # this works but the sequence is still split
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

def make_prot_df(short_prots_fp):
    prot_file = open(short_prots_fp, "r+")
    
    file_data = []
    for i in prot_file:
        file_data.append(i)
    
    proteins = {}
    for i in range(len(file_data)):
        if file_data[i].startswith(">"):
            info = file_data[i].split(" ")
            prot_id = info[0]
            proteins[prot_id] = file_data[i+1]
    
    seqs = []
    for k in proteins.keys():
        seqs.append(proteins[k])
   
    acc = ["x"] * len(proteins)
    prot_df = pd.DataFrame()
    prot_df["accession"] = acc
    prot_df["identifier"] = proteins.keys()
    prot_df["protein"] = seqs
    
    csv_fp = "sequences.csv"
    prot_df.to_csv(csv_fp)
    return csv_fp

def format_clusters(clusters_fp):
    new_clus_fp = "clusters_all_seqs.fasta"
    new_clus = open(new_clus_fp, "w+")
    clus_df = pd.read_csv(clusters_fp)
    clus_df = clus_df[["seqId", "identifier", "protein", "length", "repseq", "groupNum"]]
    clusters = {}
    clus_df = clus_df.sort_values(by = "groupNum", ascending=True)
    
    groups = clus_df["groupNum"].reset_index()
    rep_seqs = clus_df["repseq"].reset_index()
    rep_seqs = rep_seqs[rep_seqs["repseq"] == True]
    
    repseq = []
    for i in rep_seqs["index"]:
        repseq.append([clus_df["identifier"][i], clus_df["groupNum"][i]])
    
    for i in repseq:
        group = i[1]
        clusters[group] = clus_df[clus_df["groupNum"] == group].sort_values(by = "repseq", ascending = False)
    
    
    for i in clusters.keys():
        cluster = clusters[i]
        rep = cluster[cluster["repseq"] == True]
        #change these to write later
        clus_str = rep["identifier"].item() + "\n"
        new_clus.write(clus_str)
        for index, row in cluster.iterrows():
            seq_info = ""+ row["identifier"]+ " len:"+str(row["length"]) + "\n"
            prot = row["protein"]
            new_clus.write(seq_info)
            new_clus.write(prot)
        
    return new_clus_fp
    
    
def annotation_pipeline(read_file_1, read_file_2, seq_file_path, output_fp, plass_bypass_flag=False, mutation_id_window=10):
    prots_fp = ""
    if plass_bypass_flag == False:
        # run plass
        plass_string = "plass/bin/plass assemble " +read_file_1+ " " +read_file_2+ " assembly.fas ./tmp --min-seq-id 0.8 >/dev/null 2>&1"
        print("Running Plass")
        os.system(plass_string)
        prots_fp = "assembly.fas"
    elif seq_file_path is not None:
        prots_fp = seq_file_path
    
    # isolate short proteins here
    short_prots_fp = isolate_short_prots(prots_fp)
    # create csv file of protein sequences to pass in to RBiotools
    seq_csv = make_prot_df(short_prots_fp)
    # then call R code
    os.system("Rscript get_clusters.r >/dev/null 2>&1")
    print("Running Linclust")
    
    # identify mutations in clusters
    clusters_fp = format_clusters("clusters.csv")
    output_fp = "mutations.txt"
    print("Identifying mutations")
    mutation_id.id_mutations(clusters_fp, output_fp, mutation_id_window)
    
    #insert ML model steps here
    return

def main():
    parser = argparse.ArgumentParser(prog="protein_based_annotation_pipeline",description="this method executes the mutation identification and annotation pipeline.")
    subparser = parser.add_subparsers(dest="method",help="enter the the filepaths for the sequence read files and the name for the output file.")
    
    mutation_id_parser = subparser.add_parser("call_mutations",help="this method takes in sequence read data and returns mutation and annotation information for the sequences.")
    mutation_id_parser.add_argument("-i1", dest="read_file_1", type=str, help="the file path of the first sequence read file (.fastq.gz) to use in paired-end assembly")
    mutation_id_parser.add_argument("-i2", dest="read_file_2", type=str, help="the file path of the second sequence read file (.fastq.gz) to use in paired-end assembly")
    mutation_id_parser.add_argument("-s", dest="seq_file_path", type=str, help="the file path of the .fasta file of sequences if plass is to be bypassed")
    mutation_id_parser.add_argument("-o", dest="output_file_path", type=str, help="file path for the output file")
    mutation_id_parser.add_argument("-pbf",dest="plass_bypass_flag", type=int, help="0 if you would like to assemble your sequences using plass, 1 if you have pre-assembled sequences")
    mutation_id_parser.add_argument("-mw",dest="mutation_id_window", type=int, help="the window on either side of a detected mutation within which to check that the similarity threshold (90%) is met. Default=10")
    
    arguments = parser.parse_args()
    
    if arguments.method == "call_mutations":
        read_file_1 = arguments.read_file_1
        read_file_2 = arguments.read_file_2
        seq_file_path = arguments.seq_file_path
        output_file_path = arguments.output_file_path
        pbf = arguments.plass_bypass_flag
        plass_bypass_flag = False
        if pbf == 1:
            plass_bypass_flag = True
        mutation_id_window = arguments.mutation_id_window
        #call function
        annotation_pipeline(read_file_1, read_file_2, seq_file_path, output_file_path, plass_bypass_flag, mutation_id_window)
        print(read_file_1, read_file_2, seq_file_path, output_file_path, plass_bypass_flag, mutation_id_window) 
    else: 
        print("Incorrect input, please check parameters and try again")
        
if __name__ == "__main__":
    main()

