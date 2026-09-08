import os
import pandas as pd


# Set Group info
g1_name = 'PBS'
g1_sample = ['MG1-1','MG1-2','MG1-3','MG1-4','MG1-5','MG1-6']
g2_name = 'IMB064'
g2_sample = ['MG2-1','MG2-2','MG2-3','MG2-4','MG2-5','MG2-6']
versus = f'{g1_name} vs {g2_name}'

# Read ASV result
exp_result_path = 'Target_ASV_result_with_relab.tsv'
asv_result = pd.read_csv(exp_result_path, index_col=0, sep='\t')

# Get picrust2 microbe ls
pic_microbe = list(set(asv_result["PICRUSt2_old"]))
pic_microbe_id = [m.split(" | ")[0] for m in pic_microbe]

# Check Gapseq output
fn_gapseq_res = pd.read_csv("ALL_gapseq_res.tsv", sep="\t")

# make sample, pic classificaion, relabun dict
test = asv_result[["PICRUSt2_old"]+g1_sample+g2_sample].groupby("PICRUSt2_old").sum()
test = test.reset_index().melt(id_vars="PICRUSt2_old",var_name="Sample",value_name="RelAbun")
test_dict = {m:{} for m in set(test["Sample"])}
for i in range(len(test)):
    test_dict[test.loc[i,"Sample"]][test.loc[i,"PICRUSt2_old"]] = test.loc[i,"RelAbun"]

# calculate total score (gapseq flux * relabun)
fn_gapseq_res["ex_name"] = fn_gapseq_res["ex"] + " | " + fn_gapseq_res["rxn.name"]
test2 = fn_gapseq_res[["ex_name", "mtf.flux", "microbe"]].copy()
for k in test_dict.keys():
    test2[k] = test2.apply(lambda x: test_dict.get(k, {}).get(x["microbe"], 0)*x["mtf.flux"], axis=1)

# melt for draw plot
test3 = test2[["ex_name","microbe"]+g1_sample+g2_sample].copy()
test3 = test3.melt(id_vars=["ex_name","microbe"],var_name="Sample",value_name="RelAb*Flux")
test3["Group"] = [g1_name if s in g1_sample else g2_name for s in test3["Sample"]]
test3.to_csv("ALL_gapseq_res_with_relab.tsv", sep="\t", index=False)