import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set Group info
g1_name = 'PBS'
g1_sample = ['MG1-1','MG1-2','MG1-3','MG1-4','MG1-5','MG1-6']
g2_name = 'IMB064'
g2_sample = ['MG2-1','MG2-2','MG2-3','MG2-4','MG2-5','MG2-6']
versus = f'{g1_name} vs {g2_name}'

# Read silva and gg2 result
exp_result_path = 'lefse_input.xlsx'
gg2_result = pd.read_excel(exp_result_path, index_col=0, sheet_name='GG2 result')

# Set variable
txn_db = "GG2"
txn_lvl = "Genus"
txn_result = gg2_result.copy()

# Set main table
target_df = txn_result[txn_result['Level']==txn_lvl]
target_df = target_df[target_df[f'p-val(kruskall-wallis)({versus})']!='-']
target_df[f'p-val(kruskall-wallis)({versus})'] = target_df[f'p-val(kruskall-wallis)({versus})'].astype(float)
target_df = target_df[target_df[f'p-val(kruskall-wallis)({versus})']<0.05]
target_df = target_df[target_df[f'Log LDA({versus})']>2]
target_df.loc[target_df[f'Class with highest mean({versus})'] == "PBS", f'Log LDA({versus})'] *= -1

# taxa naming change
target_df['Label'] = [';_'.join([t.split('__')[-1] for t in txn.split(';')[1:]]) for txn in target_df.index]

# Prepare input data for drawing plot
lefse_df = target_df.iloc[:-1,:][['Label',
                                  f'Class with highest mean({versus})',
                                  f'Log LDA({versus})']].copy()
lefse_df[f'Log LDA({versus})'] = lefse_df[f'Log LDA({versus})'].astype(float)
lefse_df = lefse_df.sort_values(f"Log LDA({versus})", ascending=False)

# Draw plot
f, ax = plt.subplots(figsize=(4,8))

sns.set_theme(style="whitegrid")
sns.barplot(data=lefse_df, x=f"Log LDA({versus})", y="Label", hue= f"Class with highest mean({versus})",
              orient="h", palette=['tab:red','tab:blue'], linewidth=1, edgecolor="w",
              hue_order=[g1_name, g2_name], ax=ax) #'tab:orange','tab:blue','tab:gray'
ax.xaxis.grid(False)
ax.yaxis.grid(True)

plt.title(f"LEfSe\n({versus})", fontsize=12, fontweight='bold')
plt.ylabel(f"{txn_db} ({txn_lvl} level)")
plt.xlabel("Log LDA")
plt.legend(loc='lower right')

plt.savefig(f'{txn_db}_{txn_lvl}_v2.pdf', bbox_inches='tight')
plt.show()