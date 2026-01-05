import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import itertools
import random
import os


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

def explore_tuning(search_space, output_tuning, stage_name="Stage 1"):
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

    Returns
    -------
    pandas.DataFrame
        A tidy DataFrame containing hyperparameters and validation loss,
        which can be reused for further custom analysis.
    """

    tuned_params = list(search_space.keys())

    results = pd.DataFrame([
        {k: r["params"].get(k, np.nan) for k in tuned_params}
        | {"val_loss": r.get("val_loss", np.nan)}
        for r in output_tuning
    ])

    # --- Summary statistics ---
    print("=== Validation Loss Summary ===")
    print(results[["val_loss"]].describe(), "\n")

    # --- Distribution of validation loss ---
    plt.figure(figsize=(7, 5))
    sns.histplot(results["val_loss"], bins=20, kde=True)
    plt.xlabel("Validation Loss")
    plt.ylabel("Count")
    plt.title(f"Distribution of Validation Losses — {stage_name}")
    plt.show()

    # --- Correlation between hyperparameters and validation loss ---
    # only compute correlations for numeric columns
    numeric_cols = results.select_dtypes(include=[np.number])

    if "val_loss" in numeric_cols.columns and len(numeric_cols.columns) > 1:
        param_corr = numeric_cols.corr()["val_loss"].drop("val_loss")

        plt.figure(figsize=(8, 5))
        sns.barplot(x=param_corr.index, y=param_corr.values)
        plt.xticks(rotation=45)
        plt.ylabel("Correlation with Validation Loss")
        plt.title(f"Hyperparameters vs Validation Loss — {stage_name}")
        plt.show()
    else:
        print("Not enough numeric hyperparameters to compute correlations.")

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


