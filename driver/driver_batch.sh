#!/bin/bash
#SBATCH --job-name=nemo4-imhotep
#SBATCH --time=06:00:00                          # Time limit hrs:min:sec
#SBATCH --constraint="[P100*12&CL48*37]"    # Select 12 GPU nodes and 36 CascadeLake 48-core nodes
#SBATCH --nodes=49                               # 48 total nodes: 12 GPU and 36 CL48 (NEMO) 1 CL48 (XIOS)

__conda_setup="$('/home/users/shao/lus/shao/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/users/shao/lus/shao/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/users/shao/lus/shao/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/users/shao/lus/shao/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate smartsim
module load smartredis/gcc-9.3.0/0.3.0 smartsim-deps/gcc-9.3.0/0.4.0 smartsim-redis/gcc-9.3.0/0.4.0
cd /home/users/shao/dev/hpe/smartsim-nemo/driver/
python driver_eORCA025.py clustered --num_nodes 36 --tasks_per_node 48 --node_features="[CL48|CL56]"
