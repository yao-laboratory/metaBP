import os

def id_mutations(input_fp, output_fp, sim_threshold, input_folder):
    # open both files
    inputfile = open(input_fp, "r+")
    outputfile = open(output_fp, "w+")
    inputfile.seek(0, 0)
    file_data = []
    # loop through input file
    for i in inputfile:
        file_data.append(i)
    inputfile.close()
    for i in range(len(file_data)):
        cluster = []
        linedata = file_data[i].split(" ")
        hasChar = linedata[0].startswith(">")
        # find the identity sequence for a cluster
        if (len(linedata) == 1) and (hasChar == True):
            # have it identify all the sequences in the cluster and put them into a 2D array with the extra information
            # then have anonther function handle only the clusters sent in and send back a list with all the extra information to write into the file
            clusName = linedata[0]
            in_cluster = True
            counter = 0
            while in_cluster:
                i = i + 1
                next = file_data[i]
                next_split = next.split(" ")
                if i == (len(file_data)-1):
                    cluster.append(next)
                    in_cluster = False
                    clus_no_stars = remove_stars(cluster)
                    aligned_cluster = align_cluster(clus_no_stars, input_folder)
                    find_mutations(clusName, aligned_cluster, sim_threshold, outputfile)
                elif (len(next_split) == 1) and (next.startswith(">")):
                    in_cluster = False
                    i = i - 1
                    clus_no_stars = remove_stars(cluster)
                    aligned_cluster = align_cluster(clus_no_stars, input_folder)
                    find_mutations(clusName, aligned_cluster, sim_threshold, outputfile)
                else:
                    cluster.append(next)
    outputfile.close()
    print("Done")
    
    
def remove_stars(cluster):
    new_cluster = []
    for i in cluster:
        j = i.replace("*", "")
        new_cluster.append(j)
    return new_cluster

def align_cluster(cluster, input_folder):
    if len(cluster) < 3:
        return cluster
    cluster_file_path = input_folder+"/cluster.fas"
    alignment_file_path = input_folder+"/aligned.fas"
    temp_cluster_file = open(cluster_file_path, "w+")
    for i in cluster:
        temp_cluster_file.write(i)
    temp_cluster_file.close()
    os.system("clustalo/clustalo-1.2.4-Ubuntu-x86_64 -i "+input_folder+"/cluster.fas -o "+input_folder+"/aligned.fas -v >/dev/null 2>&1")
    alignment_file = open(alignment_file_path, "r+")
    file_data = []
    
    for j in alignment_file:
        file_data.append(j)
    # insert get alignmened clusters function here
    aligned_cluster = get_aligned_sequences(file_data)
    alignment_file.seek(0, 0)
    alignment_file.truncate()
    alignment_file.close()
    os.remove(cluster_file_path)
    os.remove(alignment_file_path)
    return aligned_cluster

def get_aligned_sequences(file_data):
    scrubbed_file_data = []
    for i in range(len(file_data) - 1):
        i_data = file_data[i].split(" ")
        i_next = file_data[i + 1].split(" ")
        if i == (len(file_data) - 2):
            k = i_data[0].replace("\n", "")
            scrubbed_file_data.append(k)
            scrubbed_file_data.append(file_data[i + 1])
        elif ((len(i_data) == 1) and (len(i_next) == 1)) and (i_next[0].startswith(">") == False):
            k = file_data[i].replace("\n", "")
            scrubbed_file_data.append(k)
        else:
            scrubbed_file_data.append(
                file_data[i]
            )  # this works but the sequence is still split
    sequence_info = []
    aligned_sequences = []
    sequence_info.append(scrubbed_file_data[0])
    for j in range(len(scrubbed_file_data) - 1):
        sequence = ""
        linedata = scrubbed_file_data[j].split(" ")
        hasChar = linedata[0].startswith(">")
        hasChar_next = scrubbed_file_data[j+1].startswith(">")
        # find the identifier for a sequence
        if (hasChar_next == False) and (hasChar == True):
            in_sequence = True
            counter = 0
            while in_sequence:
                j = j + 1
                next = scrubbed_file_data[j]
                next_split = next.split(" ")
                if j == len(scrubbed_file_data) - 1:
                    in_sequence = False
                if (next.startswith(">")):
                    in_sequence = False
                    sequence_info.append(scrubbed_file_data[j])
                    j = j - 1
                else:
                    sequence += next
            aligned_sequences.append(sequence)
    aligned_clusters = []
    for i in range(len(aligned_sequences)):
        aligned_clusters.append(sequence_info[i])
        aligned_clusters.append(aligned_sequences[i])
    return aligned_clusters

def find_mutations(clusName, cluster, sim_threshold, outputfile):
    if len(cluster) < 3:
        outputfile.write(clusName)
        for i in cluster:
            outputfile.write(i)
        return
    # cluster is just a list of strings
    main_seq_info = cluster[0]
    main_seq = cluster[1]
    # streamline this
    del cluster[0]
    del cluster[0]
    # first separate the data into information and sequences (have the same indices), then compare them
    seq_info = []
    sequences = []
    for x in range(len(cluster)):
        if x % 2 == 0:
            seq_info.append(cluster[x])
        else:
            sequences.append(cluster[x])
    mutations = []
    main_seq_information = {""}
    for i in sequences:
        additional_info = ""
        if main_seq != i:
            main_list = list(main_seq)
            seq_list = list(i)
            mutation_info = ""
            loop_range = 0
            if len(main_list) < len(seq_list):
                loop_range = len(main_list)
            else:
                loop_range = len(seq_list)
            for x in range(loop_range):
                check_list = check_mutation_window(
                    main_list, seq_list, sim_threshold, loop_range, x
                )
                check_left = check_list[0]
                check_right = check_list[1]
                if check_left and check_right:
                    if main_list[x] != seq_list[x]:
                        if (main_list[x] != "\n") and (seq_list[x] != "\n"):
                            mutation_info += "" + str(x + 1) + seq_list[x] + ","
                            main_seq_information.add(
                                "" + str(x + 1) + main_seq[x] + ","
                            )
            additional_info += mutation_info
        mutations.append(additional_info[:-1])
    # write the new information to the output file
    # first write the cluster name and main sequence:
    outputfile.write(clusName)
    main_seq_info = main_seq_info.replace("\n", "")
    main_seq_info += " "
    main_seq_information = list(main_seq_information)
    main_seq_information.sort()
    main_seq_mutations = ""
    for i in main_seq_information:
        main_seq_info += i
        main_seq_mutations += i
    main_seq_info += "\n"
    outputfile.write(main_seq_info)
    outputfile.write(main_seq)
    for j in range(len(sequences)):
        if mutations[j] == "":
            new_info = seq_info[j].replace("\n", "") + " " + main_seq_mutations + "\n"
        else:
            new_info = seq_info[j].replace("\n", "") + " " + mutations[j] + "\n"
        outputfile.write(new_info)
        outputfile.write(sequences[j])
        new_info = ""
        
def check_mutation_window(main_list, seq_list, sim_threshold, loop_range, x):
    left_check = False
    right_check = False
    left_match = []
    right_match = []
    # in range on the left
    if x == 0:
        left_check = True
        num_loop = 0
    elif x - sim_threshold > -1:
        num_loop = sim_threshold + 1
    else:
        num_loop = sim_threshold - x
    for i in range(2, num_loop):
        if main_list[x - i] == seq_list[x - i]:
            left_match.append(True)
        else:
            left_match.append(False)
    # in range on the right
    if x + sim_threshold < (loop_range - 1):
        num_loop = sim_threshold + 1
    else:
        num_loop = loop_range - x
    for i in range(2, num_loop):
        if main_list[x + i] == seq_list[x + i]:
            right_match.append(True)
        else:
            right_match.append(False)
    left_thresh = int(0.9 * (len(left_match)))

    count_tl = 0
    for j in left_match:
        if j == True:
            count_tl += 1
    if count_tl >= left_thresh:
        left_check = True
    right_thresh = int(0.9 * (len(right_match)))

    count_tr = 0
    for j in right_match:
        if j == True:
            count_tr += 1
    if count_tr >= right_thresh:
        right_check = True
    return [left_check, right_check]

