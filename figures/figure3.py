import matplotlib.pyplot as plt
import annotate_proteins

data_df = pd.read_csv(out_dir+"metabp_ecfreq_mice.csv", index_col="EC")
vals = data_df["Total"].iloc[0:11]
labels = data_df.index[0:11]

plt.rcdefaults()
fig, ax = plt.subplots(figsize=(15,10))
y_pos = np.arange(len(labels))
ax.barh(y_pos, vals, align='center', color="darkgreen", alpha=0.5)
ax.set_yticks(y_pos, labels=labels)
ax.invert_yaxis()
ax.set_xlabel("Number of Sequences", fontsize = 12)
ax.set_ylabel("EC", fontsize = 12)
ax.set_title("MetaBP Top 10 EC", fontsize=16)
plt.savefig("meta_ec_freq_bar.pdf",dpi=800)
plt.show()

data_df = pd.read_csv(out_dir+"eggnog_ecfreq_mice.csv", index_col="EC")
vals = data_df["Total"].iloc[0:11]
labels = data_df.index[0:11]

plt.rcdefaults()
fig, ax = plt.subplots(figsize=(15,10))
y_pos = np.arange(len(labels))
ax.barh(y_pos, vals, align='center', color="teal", alpha=0.5)
ax.set_yticks(y_pos, labels=labels)
ax.invert_yaxis()
ax.set_xlabel("Number of Sequences", fontsize = 12)
ax.set_ylabel("EC", fontsize = 12)
ax.set_title("Eggnog Top 10 EC", fontsize=16)
plt.savefig("egg_ec_freq_bar.pdf",dpi=800)
plt.show()

data_df = pd.read_csv(out_dir+"eggnog_ogfreq_mice.csv", index_col="OGs")
vals = data_df["Total"].iloc[0:11]
l = data_df.index[0:11]
labels = []
for i in l:
    labels.append(i.split("|")[-1])
    
plt.rcdefaults()
fig, ax = plt.subplots(figsize=(15,10))
y_pos = np.arange(len(labels))
ax.barh(y_pos, vals, align='center', color="teal", alpha=0.5)
ax.set_yticks(y_pos, labels=labels)
ax.invert_yaxis()
ax.set_xlabel("Number of Sequences", fontsize = 12)
ax.set_ylabel("OGs", fontsize = 12)
ax.set_title("Eggnog Top 10 OGs", fontsize=16)
plt.savefig("egg_og_freq_bar.pdf",dpi=800, bbox_inches="tight")
plt.show()

db_dir = "/work/yaolab/shared/2022_small_peptide/mitra/db_protid_split/"
db_fps = [db_dir+"prot_ids1_annot.txt", db_dir+"prot_ids2_annot.txt", db_dir+"prot_ids3_annot.txt", db_dir+"prot_ids4_annot.txt", db_dir+"prot_ids5_annot.txt"]
db = annotate_proteins.load_db(db_fps)

data_df = pd.read_csv(out_dir+"metabp_taxfreq_mice.csv", index_col="tax_id")
vals = data_df["Total"].iloc[0:20].values
l = data_df.index[0:20].values
labels = []
for i in l:
    lab = db[db["tax_id"]==i]["species"].values[0]
    labels.append(lab)
act_vals = []
act_labs = []
for i in range(len(labels)):
    if "metagenome" not in labels[i]:
        act_vals.append(vals[i])
        act_labs.append(labels[i])
        

