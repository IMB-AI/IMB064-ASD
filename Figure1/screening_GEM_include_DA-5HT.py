import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["pdf.fonttype"] = 42

ROOT_DIR = "set up root directory"

### configure
path2output = ROOT_DIR + "Analysis/IMB064_BTBR_metagenome/result/screening_GEM/"
path2gem = ROOT_DIR + "GEMreconstruction/result"
gemtool = "gapseq"
tooldate = "20231030"
medium = "gut"

targets = ["IMB064", "MM4-1A", "PS128", "DR7"]

healthy_mets = ["Taurine", "5-Aminopentanoate", "GABA", "L-Glutamine", "Spermine", "Putrescine", "(S)-2-Aminobutanoate",
                "Citrulline", "Taurodeoxycholate", "Betaine", "Taurolithocholate", "L-Serine", "1-Methylhistidine", "L-Valine",
                "Cholate", "O-Butanoylcarnitine", "Carnitine", "L-Asparagine", "Succinate",
                "L-Histidine", "O-Acetylcarnitine", "Taurochenodeoxycholate", "L-Isoleucine", "Arachidonate", "8,11,14-Icosatrienoate",
                "L-Threonine", "Docosahexaenoic acid", "Taurocholate", "Serotonin"]

disease_mets = ["L-Glutamate", "Histamine", "L-Aspartate", "Homocysteine", "Gynesine", "2-Hydroxyglutarate",
                "indole-3-propionate", "Glycolithocholate", "Phenethylamine", "Dopamine"]
# "Lecithin"

#''' 20260903 JHLEE metabolites info
healthy_mets_df = pd.DataFrame(healthy_mets)
healthy_mets_df["type"] = [ "beneficial for ASD" ] * len(healthy_mets_df)

disease_mets_df = pd.DataFrame(disease_mets)
disease_mets_df["type"] = [ "harmful for ASD" ] * len(disease_mets_df)

mets_info = pd.concat([healthy_mets_df, disease_mets_df]).rename(columns={0: "metabolites"})
mets_info.to_excel(path2output + "beneficial_harmful_mets_ASD.xlsx", index=False)
#'''


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

### plot
fig, ax = plt.subplots(figsize=(10,7))
bottom = np.zeros(len(targets))
for met, flux in bacteria_rank.to_dict(orient="list").items():
        p = ax.bar(bacteria_rank.index, abs(np.array(flux)), label=met, bottom=bottom)
        bottom += abs(np.array(flux))

ax.legend(loc="upper right")
plt.xticks(rotation=45)
plt.savefig(path2output + "/ASD_absolute_flux.png", format="png")

fig, ax = plt.subplots(figsize=(10,7))
plt.bar(total_mets.index, total_mets["total_effective_flux"].tolist())
plt.xticks(rotation=45)
plt.ylabel("total metabolic flux beneficial for ASD (mmol∗gDW−1)")
plt.savefig(path2output + "/ASD_total_effective_flux.png", format="png")
plt.savefig(path2output + "/ASD_total_effective_flux.pdf", format="pdf")

fig, ax = plt.subplots(figsize=(10,7))
plt.bar(total_mets.index, total_mets["good_flux"].tolist())
plt.xticks(rotation=45)
plt.ylabel("produced metabolic flux beneficial for ASD (mmol∗gDW−1)")
plt.savefig(path2output + "/ASD_good_flux.png", format="png")
plt.savefig(path2output + "/ASD_good_flux.pdf", format="pdf")

fig, ax = plt.subplots(figsize=(10,7))
plt.bar(total_mets.index, total_mets["bad_flux"].tolist())
plt.xticks(rotation=45)
plt.ylabel("uptaken metabolic flux harmful for ASD (mmol∗gDW−1)")
plt.savefig(path2output + "/ASD_bad_flux.png", format="png")
plt.savefig(path2output + "/ASD_bad_flux.pdf", format="pdf")

