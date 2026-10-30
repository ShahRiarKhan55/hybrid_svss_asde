from lib import *
import time
from sklearn.compose import TransformedTargetRegressor as TTR
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

CV = lambda: KFold(5, shuffle=True, random_state=0)
GRIDS = {"svr": {"svr__C": [0.1, 1, 10, 100], "svr__gamma": [0.01, 0.1, 1, 10], "svr__epsilon": [0.01, 0.1, 0.3]},
         "rf": {"randomforestregressor__min_samples_leaf": [1, 3, 5, 10], "randomforestregressor__max_features": [1.0, 0.5]},
         "mlp": {"mlpregressor__hidden_layer_sizes": [(16, 16), (32, 32), (64, 64)], "mlpregressor__alpha": [1e-3, 1e-1, 10.0]}}
EST = {"svr": SVR(kernel="rbf"), "rf": RandomForestRegressor(n_estimators=200, random_state=0), "mlp": MLPRegressor(activation="tanh", max_iter=1000, random_state=0)}


def gpr_pred(Xtr, ytr, Xte, kind):
    if kind == "add": m = fit(Xtr, ytr, "p_add")
    else:  # one-lengthscale RBF, same normalizer, bounds and optimiser settings
        m = GPy.models.GPRegression(Xtr, ytr.reshape(-1, 1), GPy.kern.RBF(2), normalizer=True)
        m.kern.lengthscale.constrain_bounded(*LS, warning=False); m.kern.variance.constrain_bounded(*VAR, warning=False); m.Gaussian_noise.variance.constrain_bounded(*NOISE, warning=False)
        m.optimize(max_iters=200, messages=False)
    mu, var = m.predict(Xte); return mu.ravel(), var.ravel()


def sk_pred(Xtr, ytr, Xte, kind):  # scaled inputs, standardised target; hyperparameters by 5-fold CV (MAE) on the training sites only
    pipe = TTR(make_pipeline(StandardScaler(), EST[kind]), transformer=StandardScaler())
    gs = GridSearchCV(pipe, {f"regressor__{k}": v for k, v in GRIDS[kind].items()}, cv=CV(), scoring="neg_mean_absolute_error", n_jobs=1).fit(Xtr, ytr)
    return gs.predict(Xte), {k.split("__", 1)[1]: v for k, v in gs.best_params_.items()}


def dgp_pred(Xtr, ytr, Xte, seed):  # layers [1,3,2], 50 inducing points, 200 iterations (predictor_comparison/run.py deepgp_fit)
    np.random.seed(seed)
    tr = pd.DataFrame(dict(X=Xtr[:, 0], Y=Xtr[:, 1], Value=ytr)); te = pd.DataFrame(dict(X=Xte[:, 0], Y=Xte[:, 1], Value=0.0))
    return run.deepgp_fit(tr, te, ind=50, iters=200)
