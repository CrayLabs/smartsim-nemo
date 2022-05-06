# SmartSim-NEMO

## Purpose 
This repository contains a driver script and SLURM batch job that can
run an eORCA_025 simulation provided by the IMHOTEP project. At a high
level this example demonstrates the following:
- Running SmartSim workflows within a SLURM batch script, launching
  the simulation and orchestrator within the same allocation
- Initializing an orchestrator in 'clustered' mode
- How to specify different resources in the same allocation for each
  component to run on using Slurm 'constraints'
- Configuring SmartSim to launch an MPMD workload, the NEMO executable
  and separate I/O-server (XIOS) processes
- Using SmartSim to modify Fortran namelists that configure the model
  prior to execution

Note: This example does require access to HPC resources due to the
high resolution of the ocean simulation. On a modest HPC platform,
it can be run for development purposes on O(300) CPUs, but for
production purposes O(3000) cores can provide reasonable simulation
throughput.

## Limited availablity
The NEMO4 code does not currently have broad distribution rights and the
and eORCA_025 configuration has not yet been published. For access,
please contact
- Andrew Shao: `Andrew.Shao at hpe.com`
- Julien Le Sommer: `julien.lesommer at univ-grenoble-alpes.fr`

## Quickstart
Assuming that the user has compiled the NEMO and XIOS binaries, has
downloaded input files and atmospheric forcing, here are the steps
to launch the simulation:
1. Modify `driver/driver_batch.sh` to configure for the user/site
   including account information, how to activate the SmartSim
   environment, the number of requested resources etc.
2. Make sure that any configurable options to the SmartSim driver
   are consistent with resources requested in the batch job (see
   Usage section below)

## Usage
A large number of runtime arguments can be passed to the SmartSim
driver, which can be be accessed by running:

```
python driver/driver_eORCA025.py clustered --help
```

which outputs the following

```
NAME
    driver_eORCA025.py clustered - Run a nemo OM4_025 simulation with a cluster of databases used for machine-learning inference

SYNOPSIS
    driver_eORCA025.py clustered <flags>

DESCRIPTION
    Run a nemo OM4_025 simulation with a cluster of databases used for machine-learning inference

FLAGS
    --nemo_num_nodes=NEMO_NUM_NODES
        Default: 48
        number of nodes allocated to nemo
    --nemo_tasks_per_node=NEMO_TASKS_PER_NODE
        Default: 48
        how many MPI ranks to be run per node
    --nemo_node_features=NEMO_NODE_FEATURES
        Default: '[CL48|CL56]'
        features of the nodes that the nemo should run on
    --xios_num_nodes=XIOS_NUM_NODES
        Default: 1
        number of nodes allocated to XIOS
    --xios_tasks_per_node=XIOS_TASKS_PER_NODE
        Default: 8
        how many MPI ranks to be run per node
    --xios_node_features=XIOS_NODE_FEATURES
        Default: '[CL48|CL56]'
        features of the nodes that the XIOS should run on
    --nemo_cfg_path=NEMO_CFG_PATH
        Default: '/lus/cls01029/shao/dev/m2lines/NEMO_4.0.6_IMHOTEP/cfgs/eORC...
        full path to the nemo configuration directory full path to the nemo configuration directory
    --nemo_forcing_path=NEMO_FORCING_PATH
        Default: '/lus/scratch/shao/model_inputs/nemo4/FORCING_ATMOSPHERIQUE/...
        full path to the nemo forcing directory full path to the nemo forcing directory
    --nemo_input_path=NEMO_INPUT_PATH
        Default: '/lus/scratch/shao/model_inputs/nemo4/ORCA025.L75/eORCA025.L...
    --xios_exe=XIOS_EXE
        Default: '/lus/cls01029/shao/dev/m2lines/xios-2.5/bin/xios_server.exe'
        path to the XIOS executable
    --job_number=JOB_NUMBER
        Default: 1
        The job number (nn_no) of the nemo experiment
    --first_time_step=FIRST_TIME_STEP
        Default: 1
        The first time step in the segment
    --last_time_step=LAST_TIME_STEP
        Default: 72
        The last time step in the segment
    --restart_flag=RESTART_FLAG
        Default: '.false.'
        If ".true.", then this experiment is a restart. Note in Fortran notation
    --restart_directory=RESTART_DIRECTORY
        Default: './'
    --output_feedback_files=OUTPUT_FEEDBACK_FILES
        Default: ''
        Path to output model feedback files used for data assimilation
    --iceberg_output_directory=ICEBERG_OUTPUT_DIRECTORY
        Default: './'
        Path to output iceberg trajectory files
    --orchestrator_port=ORCHESTRATOR_PORT
        Default: 6780
        port that the database will listen on
    --orchestrator_interface=ORCHESTRATOR_INTERFACE
        Default: 'ipogif0'
        network interface bound to the database
    --orchestrator_nodes=ORCHESTRATOR_NODES
        Default: 3
        number of orchestrator nodes to use
    --orchestrator_nemo_node_features=ORCHESTRATOR_NEMO_NODE_FEATURES
        Default: 'P100'
        node features requested for the orchestrator nodes
```        

