
# checking if there is issues with ipcw
import numpy as np
import pandas as pd
import h5py
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import requests

from pycox.evaluation import EvalSurv

# Load mb
url = "https://github.com/jaredleekatzman/DeepSurv/raw/refs/heads/master/experiments/data/metabric/metabric_IHC4_clinical_train_test.h5"
file_name = "metabric_IHC4_clinical_train_test.h5"

with open(file_name, 'wb') as f:
    f.write(requests.get(url).content)

# Read train and test
with h5py.File(file_name, 'r') as f:
    x_train = f['train']['x'][()]
    e_train = f['train']['e'][()]
    t_train = f['train']['t'][()]
    x_test = f['test']['x'][()]
    e_test = f['test']['e'][()]
    t_test = f['test']['t'][()]


print(x_train.shape)
print(e_train.shape) 
print(t_train.shape) 

df_train = pd.DataFrame(x_train)
df_train['time'] = t_train + 1e-8
df_train['status'] = e_train

df_test = pd.DataFrame(x_test)
df_test['time'] = t_test + 1e-8 
df_test['status'] = e_test

# Fit cox
cox = CoxPHFitter()
cox.fit(df_train, duration_col='time', event_col='status')

baseline_surv = cox.baseline_survival_
times_grid = baseline_surv.index.values

# Get risk
risk_scores = cox.predict_partial_hazard(df_test)

# Predict on the test sample
surv_probs = np.array([baseline_surv.values.flatten() ** rs for rs in risk_scores]).T
surv_df = pd.DataFrame(surv_probs, index=times_grid)

# Conver to right type
times_test = df_test['time'].values
events_test = df_test['status'].values.astype(int)

# Check with and without the censoring
ev_no_ipcw = EvalSurv(surv_df, times_test, events_test, censor_surv=None)
ev_ipcw = EvalSurv(surv_df, times_test, events_test, censor_surv="km")

print("C-index without IPCW:", ev_no_ipcw.concordance_td())
print("C-index with IPCW (KM):", ev_ipcw.concordance_td())
