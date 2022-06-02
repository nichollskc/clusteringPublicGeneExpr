import numpy as np
import pandas as pd

ncols = 5

for sig in ["ifn.csv", "nk_dipp.csv", "cd4_activation_green_30.csv",
            "cd4_activation_yellow_30.csv", "cd4_activation_black_30.csv", "exhaustion_down_wherry.csv"]:
    print(sig)

    genes = pd.read_csv("data/pathways/processed/" + sig)["SYMBOL"]
    print(len(genes))

    desired_length = int(np.ceil(len(genes) / ncols)) * ncols
    n_blanks = desired_length - len(genes)
    genes_grid = np.reshape(list(genes) + [""] * n_blanks, (int(desired_length / ncols), ncols))

    genes_df = pd.DataFrame(genes_grid)

    print(genes_df.to_latex(index=False))


