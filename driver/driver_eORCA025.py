"""
Example of a SmartSim driver that can be used to run a NEMO4-based
eORCA025 model. Currently, this only configures and runs a simulation.
Note a database is spun up, but it is not used by NEMO4 at this time.
"""

from glob import glob
from smartsim import Experiment
from smartsim.database import Orchestrator
from smartsim.log import log_to_file
import warnings

def create_distributed_orchestrator(
    exp,
    orchestrator_port,
    orchestrator_interface,
    orchestrator_nodes,
    orchestrator_node_features,
    walltime
    ):

    orchestrator = exp.create_database(
        port = orchestrator_port,
        interface = orchestrator_interface,
        db_nodes = orchestrator_nodes,
        time=walltime,
        threads_per_queue=2,
        batch=True)

    if orchestrator_node_features:
        orchestrator.set_batch_arg("constraint", orchestrator_node_features)
    return orchestrator

def create_NEMO_model(
        experiment,
        walltime,
        num_nodes,
        tasks_per_node,
        NEMO_cfg_path,
        NEMO_forcing_path
    ):

    NEMO_run_settings = experiment.create_run_settings(
        NEMO_cfg_path,
        run_args = { 'time': walltime }
    )

    NEMO_run_settings.set_tasks_per_node(tasks_per_node)
    NEMO_run_settings.set_tasks(num_nodes*tasks_per_node)

    NEMO_model = experiment.create_model(
        f"{NEMO_cfg_path}/BLD/bin/nemo.exe",
        run_settings   = NEMO_run_settings,
    )

    NEMO_model.attach_generator_files(
        to_configure=[
            f"{NEMO_cfg_path}/EXPREF/namelist.eORCA025.L75-IMHOTEP02",
            f"{NEMO_cfg_path}/EXPREF/namelist-ice.eORCA025.L75-IMHOTEP02"
        ]
        to_copy=f"{NEMO_cfg_path}/EXPREF/*.xml",
        to_symlink=f"{NEMO_forcing_path}/*"
    )

    return NEMO_model

def configure_NEMO_model(
    model,
    job_number,
    first_time_step,
    last_time_step,
    restart_flag,
    restart_directory,
    output_feedback_files,
    iceberg_output_directory,
    ):

    NEMO_config_options = {
        "NN_NO": job_number,
        "NIT000": first_time_step,
        "NITEND": last_time_step,
        "RESTART": restart_flag,
        "CN_DIRRST": restart_directory,
        "CN_DIAOBS": output_feedback_files,
        "CN_DIRICB": iceberg_output_directory
    }

    model.params = NEMO_config_options
    model.register_incoming_entity(model)

    return model


def NEMO_clustered_driver(
    walltime="02:00:00",
    num_nodes=25,
    tasks_per_node=45,
    NEMO_cfg_path="/home/users/shao/dev/m2lines/NEMO_4.0.6_IMHOTEP/cfgs/eORCA025.L75-IMHOTEP02",
    NEMO_forcing_path="/lus/scratch/shao/model_inputs/NEMO4/ORCA025.L75/eORCA025.L75-IMHOTEP02-I",
    job_number = 1,
    first_time_step = 1,
    last_time_step = 72,
    restart_flag = ".false.",
    restart_directory = "",
    output_feedback_files = "",
    iceberg_output_directory = "./",
    orchestrator_port=6780,
    orchestrator_interface="ipogif0",
    orchestrator_nodes=3,
    orchestrator_node_features='P100',
    configure_only=False
    ):
    """Run a NEMO OM4_025 simulation with a cluster of databases used for
    machine-learning inference

    :param walltime: how long to allocate for the run, "hh:mm:ss"
    :type walltime: str, optional
    :param num_nodes: number of nodes allocated to each model member
    :type num_nodes: int, optional
    :param tasks_per_node: how many MPI ranks to be run per node
    :type tasks_per_node: int, optional
    :param NEMO_cfg_path: full path to the NEMO configuration directory
    :type NEMO_cfg_path: str, optional
    :param NEMO_forcing_path: full path to the NEMO forcing directory
    :type NEMO_forcing_path: str, optional
    :param job_number: The job number (nn_no) of the NEMO experiment
    :type job_number: int, optional
    :param first_time_step: The first time step in the segment
    :type first_time_step: int, optional
    :param last_time_step: The last time step in the segment
    :type last_time_step: int, optional
    :param restart_flag: If ".true.", then this experiment is a restart. Note in Fortran notation
    :type restart_flag: str, optional
    :param output_feedback_files: Path to output model feedback files used for data assimilation
    :type output_feedback_files: str, optional
    :param iceberg_output_directory: Path to output iceberg trajectory files
    :type iceberg_output_directory: str, optional
    :param model_node_features: (Slurm-only) Constraints/features for the
                                    node
    :type model_node_features: str, optional
    :param orchestrator_port: port that the database will listen on
    :type orchestrator_port: int, optional
    :param orchestrator_interface: network interface bound to the database
    :type orchestrator_interface: str, optional
    :param orchestrator_nodes: number of orchestrator nodes to use
    :type orchestrator_nodes: int, optional
    :param orchestrator_node_features: (Slurm-only) node features requested for
                                       the orchestrator nodes
    :type orchestrator_node_features: str, optional
    :param configure_only: If True, only configure the experiment and return
                           the orchestrator and experiment objects
    :type configure_only: bool, optional
    """

    experiment = Experiment("AI-EKE-NEMO", launcher="auto")
    NEMO_model = create_NEMO_model(
        experiment,
        walltime,
        num_nodes,
        tasks_per_node,
        NEMO_cfg_path,
        NEMO_forcing_path
     )
    configure_NEMO_model(
        NEMO_model,
        job_number,
        first_time_step,
        last_time_step,
        restart_flag,
        restart_directory,
        output_feedback_files,
        iceberg_output_directory,
    )
    orchestrator = create_distributed_orchestrator(
        experiment,
        orchestrator_port,
        orchestrator_interface,
        orchestrator_nodes,
        orchestrator_node_features,
        walltime
    )

    experiment.generate( NEMO_model, orchestrator, overwrite=True )
    if configure_only:
        return experiment, NEMO_model, orchestrator
    else:
        experiment.start(NEMO_model, orchestrator, summary=True)
        experiment.stop(orchestrator)

if __name__ == "__main__":
    import fire
    log_to_file("./NEMO_driver.log", log_level='info')
    fire.Fire({
        "clustered":NEMO_clustered_driver
    })
