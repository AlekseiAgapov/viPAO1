#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap


matrix = pd.read_csv("present_max50pct_10kb_mq20_primary_matrix.tsv", sep="\t")
matrix = matrix.set_index("sample")
matrix = matrix.apply(pd.to_numeric)

width = max(8, 0.34 * matrix.shape[1])
height = max(3, 0.45 * matrix.shape[0] + 1.2)

plt.figure(figsize=(width, height))
sns.heatmap(
    matrix,
    cmap=ListedColormap(["white", "#2ca25f"]),
    linewidths=0.5,
    linecolor="black",
    cbar=False,
    square=True,
)

plt.xlabel("Phage ID")
plt.ylabel("")
plt.title("Phage presence/absence")
plt.tight_layout()
plt.savefig("phage_presence_absence.pdf", bbox_inches="tight")
