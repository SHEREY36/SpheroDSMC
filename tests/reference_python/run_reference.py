import sys, yaml, copy, numpy as np, json, time
sys.path[:0]=["/home/muhammed/Documents/Thesis/DSMC_V2/contracts/python","/home/muhammed/Documents/Thesis/DSMC_V2/Coll_Models_v2/src","/home/muhammed/Documents/Thesis/DSMC_V2/DSMC_0D_v2/src"]
import os; os.chdir("/home/muhammed/Documents/Thesis/DSMC_V2")
from dsmc_v2.simulation import run_simulation
cfg=yaml.safe_load(open("DSMC_0D_v2/config/encounter_unit_model_v1.yaml"))
seed=int(sys.argv[1]); out=sys.argv[2]
N=10000; AR=2.0; rho=0.763944
Vp=np.pi/6+np.pi*0.25*(AR-1); m=rho*Vp; phi=0.01; V=N*Vp/phi; L=V**(1/3)
cfg["particle"].update({"AR":AR,"radius":0.5,"mass":float(m)})
cfg["system"].update({"alpha":0.8,"kTt":1.0,"kTr":1.0,"phi":phi,"domain":[L,L,L],"particle_count":N})
cfg["time"].update({"dt":0.16,"dtau":0.05,"t_end":1e9,"tau_end":10.0})
cfg.setdefault("simulation",{})["exact_initial_temperatures"]=True
t0=time.time()
diag=run_simulation(cfg, seed, out)
print("wall",time.time()-t0, {k:diag[k] for k in ("collisions","cpp","encounter_clock_mean_ratio")})
