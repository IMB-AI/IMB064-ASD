import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["pdf.fonttype"] = 42

ROOT_DIR = "set up ROOT directory"

### configure
path2output = ROOT_DIR + "Analysis/IMB064_BTBR_metagenome/result/screening_GEM_glu_gaba_proportion/"
path2gem = ROOT_DIR + "GEMreconstruction/result"
gemtool = "gapseq"
tooldate = "20231030"
medium = "gut"
#targets = ["IMB19", "IMB002", "IMB012", "IMB013", "IMB005", "IMB020", "IMB011", "IMB023", "IMB017",
#           "IMB019", "IMB046", "IMB048", "IMB049", "IMB052", "IMB045", "IMB015", "IMB016", "IMB064"]
targets = [ "IMB064", "MM4-1A", "PS128", "DR7"]

healthy_mets = ["Taurine", "5-Aminopentanoate", "GABA", "L-Glutamine", "Spermine", "Putrescine", "(S)-2-Aminobutanoate",
                "Citrulline", "Taurodeoxycholate", "Betaine", "Taurolithocholate", "L-Serine", "1-Methylhistidine", "L-Valine",
                "Cholate", "O-Butanoylcarnitine", "Carnitine", "L-Asparagine", "Succinate",
                "L-Histidine", "O-Acetylcarnitine", "Taurochenodeoxycholate", "L-Isoleucine", "Arachidonate", "8,11,14-Icosatrienoate",
                "L-Threonine", "Docosahexaenoic acid", "Taurocholate", "Serotonin"]
healthy_mets_target = ["GABA"]
# "C16 sphingomyelin"
disease_mets = ["L-Glutamate", "Histamine", "L-Aspartate", "Homocysteine", "Gynesine", "2-Hydroxyglutarate",
                "indole-3-propionate", "Glycolithocholate", "Phenethylamine", "Dopamine"]
disease_mets_target = ["L-Glutamate"]
#
# "Lecithin"


### integrate metabolites profile of individual bacteria
mets = [ pd.read_csv(f"{path2gem}/{strain}/{gemtool}/{medium}/metabolites_produced_uptaken.txt", sep="\t")
        for strain in targets ]
for df, dataset in zip(mets, targets): 
        df["dataset"] = dataset
        df["metabolites"] = [ rxn.split("-e0 ")[0] for rxn in df["rxn.name"] ]

mets = [ df.loc[ df["metabolites"].isin(healthy_mets + disease_mets), :] for df in mets ]
mets_df = pd.concat(mets, axis=0)
mets_df["abs.mtf.flux"] = abs(mets_df["mtf.flux"])
# output = mets_df.groupby("dataset")["abs.mtf.flux"].sum().to_frame().sort_values(by="abs.mtf.flux", ascending=False)


### good & bad flux calculation
#bacteria_rank = pd.DataFrame(output.index)
bacteria_rank = pd.DataFrame(targets, columns=["dataset"])
for met in healthy_mets + disease_mets:
        df = mets_df.loc[mets_df["metabolites"] == met, ["dataset", "metabolites", "mtf.flux"]]
        df = pd.pivot_table(df, index="dataset", columns="metabolites", values="mtf.flux").reset_index()
        bacteria_rank = pd.merge(bacteria_rank, df, on="dataset", how="left")
bacteria_rank = bacteria_rank.fillna(0).set_index("dataset")


good_mets = bacteria_rank.loc[:, bacteria_rank.columns.isin(healthy_mets)]
good_mets_sum = good_mets.sum(axis=1).to_frame().rename(columns={0: "good_flux"})
bad_mets = bacteria_rank.loc[:, bacteria_rank.columns.isin(disease_mets)]
bad_mets_sum = bad_mets.sum(axis=1).to_frame().rename(columns={0: "bad_flux"})


total_mets = pd.concat([good_mets_sum,bad_mets_sum], axis=1)
#print(total_mets.sort_values(by="good_flux", ascending=False))
#print(total_mets.sort_values(by="bad_flux", ascending=True))
total_mets["total_effective_flux"] = [ gf + (-1)*bf for gf, bf in zip(total_mets["good_flux"], total_mets["bad_flux"]) ]
total_mets = total_mets.sort_values(by="total_effective_flux", ascending=False)



### plot - good metabolites
fig, ax = plt.subplots(figsize=(10, 7))
pos_bottom = np.zeros(len(good_mets)); neg_bottom = np.zeros(len(good_mets))
for met in good_mets.columns:
        values = good_mets[met].values
        
        positives = np.where(values > 0, values, 0)
        negatives = np.where(values < 0, values, 0)
        
        if met in healthy_mets_target:
                color = "crimson"
        else:
                color = "lightgray"

        ax.bar(np.arange(len(good_mets.index)), positives, 0.7, bottom=pos_bottom, color=color, label=met if met in healthy_mets_target else None)
        ax.bar(np.arange(len(good_mets.index)), negatives, 0.7, bottom=neg_bottom, color=color)

        pos_bottom += positives
        neg_bottom += negatives

ax.set_xticks(np.arange(len(good_mets)))
ax.set_xticklabels(good_mets.index, rotation=45)
plt.ylabel("flux of beneficial metabolites (mmol∗gDW−1)")
ax.axhline(0, color='black', linewidth=0.8)
ax.legend(title="Highlighted")

plt.tight_layout()
plt.savefig(path2output + "/gaba_flux_proportion.png", format="png")
plt.savefig(path2output + "/gaba_flux_proportion.pdf", format="pdf")


### plot - bad metabolites
fig, ax = plt.subplots(figsize=(10, 7))
pos_bottom = np.zeros(len(bad_mets)); neg_bottom = np.zeros(len(bad_mets))
for met in bad_mets.columns:
        values = bad_mets[met].values
        
        positives = np.where(values > 0, values, 0)
        negatives = np.where(values < 0, values, 0)
        
        if met in disease_mets_target:
                color = "crimson"
        else:
                color = "lightgray"

        ax.bar(np.arange(len(bad_mets.index)), positives, 0.7, bottom=pos_bottom, color=color, label=met if met in disease_mets_target else None)
        ax.bar(np.arange(len(bad_mets.index)), negatives, 0.7, bottom=neg_bottom, color=color)

        pos_bottom += positives
        neg_bottom += negatives

ax.set_xticks(np.arange(len(bad_mets)))
ax.set_xticklabels(bad_mets.index, rotation=45)
plt.ylabel("flux of not beneficial metabolites (mmol∗gDW−1)")
ax.axhline(0, color='black', linewidth=0.8)
ax.legend(title="Highlighted")

plt.tight_layout()
plt.savefig(path2output + "/glu_flux_proportion.png", format="png")
plt.savefig(path2output + "/glu_flux_proportion.pdf", format="pdf")




### plot
# fig, ax = plt.subplots(figsize=(10,7))
# bottom = np.zeros(len(targets))
# for met, flux in bacteria_rank.to_dict(orient="list").items():
#         p = ax.bar(bacteria_rank.index, abs(np.array(flux)), label=met, bottom=bottom)
#         bottom += abs(np.array(flux))

# ax.legend(loc="upper right")
# plt.xticks(rotation=45)
# plt.savefig(path2output + "/ASD_absolute_flux_gabaglu.png", format="png")

# fig, ax = plt.subplots(figsize=(10,7))
# plt.bar(total_mets.index, total_mets["total_effective_flux"].tolist())
# plt.xticks(rotation=45)
# plt.ylabel("total metabolic flux beneficial for ASD (mmol∗gDW−1)")
# plt.savefig(path2output + "/ASD_total_effective_flux_gabaglu.png", format="png")
# plt.savefig(path2output + "/ASD_total_effective_flux_gabaglu.pdf", format="pdf")

# fig, ax = plt.subplots(figsize=(10,7))
# plt.bar(total_mets.index, total_mets["good_flux"].tolist())
# plt.xticks(rotation=45)
# plt.ylabel("flux of GABA (mmol∗gDW−1)")
# plt.savefig(path2output + "/gaba_flux.png", format="png")
# plt.savefig(path2output + "/gaba_flux.pdf", format="pdf")

# fig, ax = plt.subplots(figsize=(10,7))
# plt.bar(total_mets.index, total_mets["bad_flux"].tolist())
# plt.xticks(rotation=45)
# plt.ylabel("flux of glutamate (mmol∗gDW−1)")
# plt.savefig(path2output + "/glu_flux.png", format="png")
# plt.savefig(path2output + "/glu_flux.pdf", format="pdf")

