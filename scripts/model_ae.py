import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import mean_squared_error, mean_absolute_error
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from tensorflow.keras.optimizers import Adam, RMSprop
from tensorflow.keras.utils import to_categorical
import tensorflow as tf

# --------------------------------------------------------------------
# -------------------------- Hyperparameters -------------------------
# --------------------------------------------------------------------

params = {
    # Encoder
    "enc_units1": 512,
    "enc_units2": 0,
    "enc_units3": 0,
    "drop1": 0.2,

    # Latent space
    "latent_dim": 32,

    # Decoder 
    "dec_units1": 128,
    "dec_units2": 0,
    "dec_units3": 0,

    # Shared
    "seed": 22,
    "activation": "relu",
    "optimizer": "adam",
    "lr": 1e-3,
    "epochs": 200,
    "batch_size": 16,
    "val_split": 0.2
} # this grid can handle asymmetric encoder/decoder


np.random.seed(params["seed"])
tf.random.set_seed(params["seed"])

# --------------------------------------------------------------------
# ---- Data loading and preprocessing --------------------------------
# --------------------------------------------------------------------

def load_data(path="../data/merged_expression.csv", test_size=0.2, seed=params["seed"]):
    X = pd.read_csv(path, index_col=0)
    y = X.index

    # train test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, shuffle=True, random_state=seed
    )

    # scaling
    X_train = StandardScaler().fit_transform(X_train)
    X_test = StandardScaler().fit_transform(X_test)

    return X_train, X_test, y_train, y_test

# --------------------------------------------------------------------
# ---------------------- Building the model --------------------------
# --------------------------------------------------------------------

def build_ae(X_train, params):
    # Input layer
    inputs = Input(shape=(X_train.shape[1],), name="input_layer")
    x = inputs

    # -------- Encoder --------
    x = Dense(params["enc_units1"], activation=params["activation"])(x)
    if params["drop1"] > 0:
        x = Dropout(params["drop1"])(x)
    x = BatchNormalization()(x)

    if params["enc_units2"] > 0:
        x = Dense(params["enc_units2"], activation=params["activation"])(x)
        if params["drop1"] > 0:
            x = Dropout(params["drop1"])(x)
        x = BatchNormalization()(x)

    if params["enc_units3"] > 0:
        x = Dense(params["enc_units3"], activation=params["activation"])(x)
        if params["drop1"] > 0:
            x = Dropout(params["drop1"])(x)
        x = BatchNormalization()(x)

    # Latent layer (bottleneck)
    latent = Dense(params["latent_dim"], activation="linear", name="latent")(x)
    x = latent

    # -------- Decoder --------
    if params["dec_units1"] > 0:
        x = Dense(params["dec_units1"], activation=params["activation"])(x)

    if params["dec_units2"] > 0:
        x = Dense(params["dec_units2"], activation=params["activation"])(x)

    if params["dec_units3"] > 0:
        x = Dense(params["dec_units3"], activation=params["activation"])(x)

    # Output layer — reconstruct input
    outputs = Dense(X_train.shape[1], activation="linear")(x)

    # Create model
    model = Model(inputs=inputs, outputs=outputs, name="autoencoder")

    # -------- Optimizer --------
    if params["optimizer"].lower() == "rmsprop":
        opt = RMSprop(learning_rate=params["lr"])
    else:
        opt = Adam(learning_rate=params["lr"])

    model.compile(optimizer=opt, loss="mse")

    return model


# --------------------------------------------------------------------
# -------------------- Callbacks -------------------------------------
# --------------------------------------------------------------------

def get_callbacks(params):
    return [
        EarlyStopping(monitor="val_loss", patience=15, restore_best_weights=True),
        ReduceLROnPlateau(monitor="val_loss", factor=0.5, patience=5, min_lr=1e-5)
    ]

# --------------------------------------------------------------------
# ---- Training ------------------------------------------------------
# --------------------------------------------------------------------

def train_ae(model, X_train, params):
    
    callbacks = get_callbacks(params)

    history = model.fit(
        x=X_train,
        y=X_train,
        epochs=params["epochs"],
        batch_size=params["batch_size"],
        validation_split=params["val_split"],
        callbacks=callbacks,
        verbose=0
    )
    return history

# --------------------------------------------------------------------
# ---- Evaluation ----------------------------------------------------
# --------------------------------------------------------------------

def evaluate_ae(model, X_test):
    X_rec = model.predict(X_test, verbose=0)
    mse = mean_squared_error(X_test, X_rec)
    mae = mean_absolute_error(X_test, X_rec)
    rmse = np.sqrt(mse)

    metrics = {
        "mse": mse,
        "mae": mae,
        "rmse": rmse
    }

    print(
        f"Reconstruction metrics — "
        f"MSE={mse:.4f}, MAE={mae:.4f}, RMSE={rmse:.4f}"
    )

    return metrics

# --------------------------------------------------------------------
# -------------------- Main execution -------------------------------
# --------------------------------------------------------------------

if __name__ == "__main__":
    X_train, X_test, y_train, y_test = load_data()
    model = build_ae(X_train, params)
    history = train_ae(model, X_train, params)
    evaluate_ae(model, X_test)
