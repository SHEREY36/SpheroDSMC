#!/bin/bash
# Serial wall time for the first contact per particle (tau_c = 1), from
# T_tr = T_rot = 1, for the DEM box sizes (64 d, phi = 0.01) and N = 250000.
# Writes runs/bench/summary.txt: AR alpha N wall_seconds_to_tau_c1 encounters
# Run on an otherwise idle machine. ARS and ALPHAS may be set in the environment.
ARS=${ARS:-"1.1 1.5 2.0 2.5 3.0"}
ALPHAS=${ALPHAS:-"0.50 0.80 0.95"}
mkdir -p runs/bench
echo "# AR alpha N wall_s_to_tau_c1 encounters" > runs/bench/summary.txt
declare -A NDEM=([1.1]=4354 [1.25]=3642 [1.5]=2861 [2.0]=2003 [2.5]=1541 [3.0]=1252)
for AR in $ARS; do
  for a in $ALPHAS; do
    [ -f models/model_AR${AR}_a${a}.npz ] || continue
    for N in ${NDEM[$AR]} 250000; do
      p=runs/bench/AR${AR}_a${a}_N$N
      cat > $p.nml <<NML
&run
  model_file = 'models/model_AR${AR}_a${a}.npz'
  output_prefix = '$p'
  n_particles = $N
  seed = 7
  tau_end = 1.0
/
NML
      ./bin/dsmc0d $p.nml > /dev/null
      w=$(grep "wall seconds to tau_c = 1" $p.log | awk '{print $NF}')
      e=$(grep "^encounters:" $p.log | awk '{print $NF}')
      echo "$AR $a $N $w $e" | tee -a runs/bench/summary.txt
    done
  done
done
