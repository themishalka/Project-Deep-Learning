import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import itertools
import random
import os

import matplotlib as mpl
mpl.rcParams["image.cmap"] = "viridis"



### ------------- Defining and sampling parameter combinations ------------

def sample_combinations(search_space, sample_ratio=0.6, seed=22):
    """
    Sample a subset of hyperparameter combinations from a full grid.

    Parameters
    ----------
    search_space : dict
        A dictionary where each key is a hyperparameter name and each value
        is a list of possible values

    sample_ratio : float, optional (default=0.6)
        Fraction of all possible combinations to randomly sample.
        Must be in the range (0, 1].

    seed : int or None, optional (default=22)
        Random seed for reproducible sampling. If None, sampling is not seeded.

    Returns
    -------
    list of dict
        A list of sampled combinations, where each element is a dictionary
        mapping hyperparameter names to a chosen value.
    """

    if seed is not None:
        random.seed(seed)

    # generate all combinations
    all_combinations = [
        dict(zip(search_space.keys(), values))
        for values in itertools.product(*search_space.values())
    ]

    if not 0 < sample_ratio <= 1:
        raise ValueError("sample_ratio must be in the range (0, 1].")

    n_samples = max(1, int(sample_ratio * len(all_combinations)))

    # if sample_ratio == 1, just return everything
    if n_samples >= len(all_combinations):
        return all_combinations

    return random.sample(all_combinations, k=n_samples)


### ------------- Exploring the hyperparameter tuning output --------------

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def explore_tuning(search_space, 
                   output_tuning, 
                   stage_name="Stage 1", 
                   param_name_map=None):
    """
    Explore and visualize tuning results for a given search stage.

    Parameters
    ----------
    search_space : dict
        Dictionary defining the hyperparameters that were tuned.
        Keys are hyperparameter names and values are the candidate values
        explored during the search. 

    output_tuning : list of dict
        List of run results, where each element contains at least:
        - "params": dict of hyperparameter values used for that run
        - "val_loss": float validation loss of the run

    stage_name : str, optional (default="Stage 1")
        Label used in figure titles to distinguish exploration phases.

    param_name_map : str, optional (default=None)
        Explicit name of parameters for clearer visualisation

    Returns
    -------
    pandas.DataFrame
        A tidy DataFrame containing hyperparameters and validation loss,
        which can be reused for further custom analysis.
    """
    tuned_params = list(search_space.keys())

    results = pd.DataFrame([
        {
            **{k: r["params"].get(k, np.nan) for k in tuned_params},
            "mse": r.get("mse", np.nan),
            "mae": r.get("mae", np.nan),
            "rmse": r.get("rmse", np.nan),
        }
        for r in output_tuning
    ])

    print("=== Metric Summary ===")
    print(results[["mse", "mae", "rmse"]].describe(), "\n")
    
    # --- Scatter: MSE vs MAE ---
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=results,
        x="mse", y="mae",
        palette="viridis", s=70
    )
    plt.xlabel("MSE (Validation)")
    plt.ylabel("MAE (Validation)")
    #plt.title(f"MSE vs MAE across hyperparameter runs — {stage_name}")
    plt.tight_layout()
    plt.show()

    # --- Correlation: hyperparameters vs MSE ---
    numeric_cols = results.select_dtypes(include=[np.number])
    numeric_cols = numeric_cols.drop(columns=[col for col in ['mae', 'rmse'] if col in numeric_cols.columns])


    if "mse" not in numeric_cols.columns:
        print("No MSE values found — cannot compute correlations.")
        return results

    param_cols = [c for c in numeric_cols.columns if c != "mse"]
    if not param_cols:
        print("No numeric hyperparameters to correlate with MSE.")
        return results

    corr = numeric_cols.corr()["mse"].drop("mse")

    # Apply readable label mapping
    if param_name_map is not None:
        corr.index = [param_name_map.get(x, x) for x in corr.index]

    plt.figure(figsize=(8, 6))
    sns.barplot(x=corr.index, y=corr.values)
    plt.xticks(rotation=40, ha="right")
    plt.ylabel("Correlation with MSE")
    plt.xlabel("Hyperparameter")
    #plt.title(f"Hyperparameters vs MSE — {stage_name}")
    plt.tight_layout()
    plt.show()

    return results

    


### ---------------------------- Saving model ------------------------------

def save_model(
    run_result,
    save_dir,
    file_name
):
    """
    Save model and hyperparameters

    Parameters
    ----------
    run_result : dict
        Dictionary from tuning run, expected keys:
        - "model" : trained Keras model
        - "params": dict of hyperparameters
    save_dir : str
        Directory where artifacts will be written.
    file_name : str
        Base file name (no extension).
    """

    os.makedirs(save_dir, exist_ok=True)

    model_path  = os.path.join(save_dir, f"{file_name}.keras")
    params_path = os.path.join(save_dir, f"{file_name}_params.npy")

    model = run_result["model"]

    # ---- Save model + params ----
    model.save(model_path)
    np.save(params_path, run_result["params"], allow_pickle=True)

    print(f"Saved model to {model_path}")


# ----------------- Creating label map for visualisation -------------------

def create_map(X, meta, ref_col, mapping_col):
    tumor_type_map = meta.set_index(ref_col)[mapping_col]
    tumor_for_test = X.index.map(tumor_type_map)
    tumor_for_test = tumor_for_test.fillna("NA") 
    tumor_map = tumor_for_test.astype("category")
    colors_tumor = tumor_map.codes

    return tumor_map, colors_tumor


# ---------------------------- PCA projection -------------------------------

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def plot_pca(Z, label_map, n_components=2, title="PCA of Latent Space",
             cmap_name="viridis", s=40, random_state=22):

    pca = PCA(n_components=n_components, random_state=random_state)
    Z_2d = pca.fit_transform(Z)

    # Label handling
    if hasattr(label_map, "cat"):
        categories = label_map.cat.categories
        color_codes = label_map.cat.codes.values
    else:
        categories, inv = np.unique(label_map, return_inverse=True)
        color_codes = inv

    cmap = getattr(plt.cm, cmap_name)
    norm = matplotlib.colors.Normalize(vmin=color_codes.min(),
                                       vmax=color_codes.max())

    # Plot
    plt.figure(figsize=(8,5))
    plt.scatter(Z_2d[:,0], Z_2d[:,1],
                c=color_codes, cmap=cmap, norm=norm, s=s)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)

    # Legend
    handles = [
        plt.Line2D([], [], marker="o", linestyle="",
                   color=cmap(norm(i)), label=str(cat))
        for i, cat in enumerate(categories)
    ]
    plt.legend(handles=handles, title="Label", loc="upper right")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.show()

    print("Explained variance ratio:", pca.explained_variance_ratio_)
    return Z_2d, pca



# ---------------------------- UMAP projection ------------------------------

from umap.umap_ import UMAP
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def plot_umap(
    Z,
    label_map,
    n_neighbors=20,
    min_dist=0.3,
    title="UMAP of Latent Space",
    cmap_name="viridis",
    s=40,
    random_state=22,
):
    """
    Z : array-like, shape (n_samples, n_latent)
        Latent space embeddings.
    label_map : pandas.Categorical or Series
        Labels for coloring (e.g. tumor_map, grade_map).
    """

    # Run UMAP
    reducer = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=random_state,
    )
    Z_umap = reducer.fit_transform(Z)

    # Convert labels → color codes
    if hasattr(label_map, "cat"):
        categories = label_map.cat.categories
        color_codes = label_map.cat.codes.values
    else:
        # fallback for plain arrays / series
        categories, inv = np.unique(label_map, return_inverse=True)
        color_codes = inv

    # Shared colormap + normalization
    cmap = getattr(plt.cm, cmap_name)
    norm = matplotlib.colors.Normalize(
        vmin=color_codes.min(),
        vmax=color_codes.max()
    )

    # ---- Plot ----
    plt.figure(figsize=(8, 5))
    plt.scatter(
        Z_umap[:, 0],
        Z_umap[:, 1],
        c=color_codes,
        cmap=cmap,
        norm=norm,
        s=s,
    )

    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title(title)

    # Legend using exact category→code mapping
    handles = [
        plt.Line2D(
            [], [], marker="o", linestyle="",
            color=cmap(norm(i)),
            label=str(cat),
        )
        for i, cat in enumerate(categories)
    ]

    plt.legend(handles=handles, title="Label", loc="upper right")
    plt.tight_layout()
    plt.show()

    return Z_umap


# ------------------------- t-SNE projection ----------------------------

from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def plot_tsne(Z, label_map, perplexity=30, n_iter=1000,
              learning_rate="auto", init="pca",
              title="t-SNE of Latent Space",
              cmap_name="viridis", s=40, random_state=22):

    tsne = TSNE(n_components=2,
                perplexity=perplexity,
                max_iter=n_iter,
                learning_rate=learning_rate,
                init=init,
                random_state=random_state)

    Z_2d = tsne.fit_transform(Z)

    # Label handling
    if hasattr(label_map, "cat"):
        categories = label_map.cat.categories
        color_codes = label_map.cat.codes.values
    else:
        categories, inv = np.unique(label_map, return_inverse=True)
        color_codes = inv

    cmap = getattr(plt.cm, cmap_name)
    norm = matplotlib.colors.Normalize(vmin=color_codes.min(),
                                       vmax=color_codes.max())

    # Plot
    plt.figure(figsize=(8,5))
    plt.scatter(Z_2d[:,0], Z_2d[:,1],
                c=color_codes, cmap=cmap, norm=norm, s=s)

    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.title(title)

    # Legend
    handles = [
        plt.Line2D([], [], marker="o", linestyle="",
                   color=cmap(norm(i)), label=str(cat))
        for i, cat in enumerate(categories)
    ]
    plt.legend(handles=handles, title="Label", loc="upper right")
    plt.tight_layout()
    plt.show()

    return Z_2d, tsne

