#%% Libraries

import os
import pandas as pd
import matplotlib.pyplot as plt


#%% Function for writing namelist file

def write_namelist(ALBEDO, CANMOD, CONDCT, DENSTY, EXCHNG, HYDROL):
    
  nlst = f"""
&nam_grid
  NNx = 1,
  NNy = 1,
  NNsmax = 3,
  NNsoil = 4,
/
&nam_layers
  DDzsnow = 0.1, 0.2, 0.4,
  DDzsoil = 0.1, 0.2, 0.4, 0.8,
/
&nam_driving
  zzT = 10,
  zzU = 10,
  met_file = 'data/input_SLF_5WJ.txt',
  out_file = 'data/output_SLF_5WJ.txt',
/
&nam_modconf
  NALBEDO = {ALBEDO},
  NCANMOD = {CANMOD},
  NCONDCT = {CONDCT},
  NDENSTY = {DENSTY},
  NEXCHNG = {EXCHNG},
  NHYDROL = {HYDROL},
  NSNFRAC = 3,
  NRADSBG = 0,
  NZOFFST = 0,
  NOSHDTN = 1,
  LHN_ON  = .FALSE.,
  LFOR_HN = .TRUE.,
/
&nam_modpert
  LZ0PERT = .FALSE.,
/
&nam_modtile
  CTILE = 'open',
  rtthresh = 0.1,
/
&nam_results
   CLIST_DIAG_RESULTS = 'rotc', 'hsnt', 'swet', 'slqt', 'swtb', 'swtd', 'lwtr', 'romc', 'sbsc',
   CLIST_STATE_RESULTS = 'tsfe', 'scfe',
/
&nam_location
  fsky_terr = 0.9483,
  slopemu = 0.0122,
  xi = 0,
  Ld = 1,
  lat = 46.8296,
  lon = 9.8092,
  dem = 2540,
  pmultf = 1,
  fveg = 0,
  hcan = 0,
  lai = 0,
  vfhp = 1, 
  fves = 0,
/
  """
  
  nlst_file = open("nlst_tmp.nam", "w")
  nlst_file.write(nlst)
  nlst_file.close()


#%% Function for running the model

def run_fsm(ALBEDO=2, CANMOD=0, CONDCT=1, DENSTY=3, EXCHNG=1, HYDROL=2):

  write_namelist(ALBEDO, CANMOD, CONDCT, DENSTY, EXCHNG, HYDROL)

  status = os.system('FSM_TXT.exe nlst_tmp.nam')

  os.remove('nlst_tmp.nam')

  if status != 0:
    print("Error")
    return None

  names = ['year','month','day','hour','Ds','fsnow','SWE','Tsrf','Nsnow']
  df = pd.read_csv('data/output_SLF_5WJ.txt', delim_whitespace=True, names=names)

  return df


#%% Function for plotting results

def plot(df, label, field='SWE', ylabel='Snow water equivalent (mm)'):

  plt.plot(df[field], label=label)
  plt.ylabel(ylabel)


#%% Run model for different albedo options

df = run_fsm(ALBEDO=0)
plot(df, 'ALBEDO=0')
df = run_fsm(ALBEDO=1)
plot(df, 'ALBEDO=1')
df = run_fsm(ALBEDO=2)
plot(df, 'ALBEDO=2')
plt.legend()


#%% Run model for different canopy options

df = run_fsm(CANMOD=0)
plot(df, 'CANMOD=0')
df = run_fsm(CANMOD=1)
plot(df, 'CANMOD=1')
plt.legend()


#%% Run model for different conductivity options

df = run_fsm(CONDCT=0)
plot(df, 'CONDCT=0')
df = run_fsm(CONDCT=1)
plot(df, 'CONDCT=1')
plt.legend()


#%% Run model for different density options

df = run_fsm(DENSTY=0)
plot(df, 'DENSTY=0')
df = run_fsm(DENSTY=1)
plot(df, 'DENSTY=1')
df = run_fsm(DENSTY=2)
plot(df, 'DENSTY=2')
df = run_fsm(DENSTY=3)
plot(df, 'DENSTY=3')
plt.legend()


#%% Run model for different density options

df = run_fsm(EXCHNG=0)
plot(df, 'EXCHNG=0')
df = run_fsm(EXCHNG=1)
plot(df, 'EXCHNG=1')
df = run_fsm(EXCHNG=2)
plot(df, 'EXCHNG=2')
plt.legend()


#%% Run model for different hydraulics options

df = run_fsm(HYDROL=0)
plot(df, 'HYDROL=0')
df = run_fsm(HYDROL=1)
plot(df, 'HYDROL=1')
df = run_fsm(HYDROL=2)
plot(df, 'HYDROL=2')
plt.legend()
