#
# mlips.py
#
# This module provides Machine Learning Interatomic Potential (MLIPS) calculators
# and related functionality for the auto-kappa package.
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os
import logging
import numpy as np
from pathlib import Path

logger = logging.getLogger(__name__)


class MLIPSCalculatorFactory:
    """
    Supported models:
        - eSEN: EquiformerV2-based universal potential
        - MACE: Multi-Atomic Cluster Expansion potential
    """
    
    # Default model paths and configurations
    DEFAULT_CONFIGS = {
        'esen': {
            'model_dir': '/home/wangzeyu/home/research/ML-IPS/matbench-discovery/models/eSEN/checkpoints',
            'model_name': 'esen_30m_oam.pt',
            'use_gpu': True,
        },
        'mace': {
            'model_name': 'mace-mpa-0-medium',
            'model_url_template': 'https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/{model_name}.model',
            'use_gpu': True,
            'enable_cueq': False,
        }
    }
    
    @classmethod
    def create_calculator(cls, model_name, directory=None, **kwargs):
        """Create a MLIPS calculator based on the model name.
        
        Args
        ----
        model_name : str
            Name of the MLIPS model ("esen", "mace")
        directory : str, optional
            Output directory for calculations
        **kwargs : dict
            Additional model-specific arguments that override defaults
            
        Returns
        -------
        calculator : MLIPS calculator object
            The initialized calculator
            
        Raises
        ------
        ValueError
            If the model name is not supported
        """
        model_name_lower = model_name.lower()
        
        if model_name_lower == "esen":
            return cls._create_esen_calculator(directory=directory, **kwargs)
        elif model_name_lower == "mace":
            return cls._create_mace_calculator(directory=directory, **kwargs)
        else:
            raise ValueError(
                f"Unsupported MLIPS model: {model_name}. "
                f"Supported models: {list(cls.DEFAULT_CONFIGS.keys())}"
            )
    
    @classmethod
    def _create_esen_calculator(cls, directory=None, **kwargs):
        """Create eSEN (EquiformerV2) calculator.
        
        Args
        ----
        directory : str, optional
            Output directory
        **kwargs : dict
            Override default configurations:
            - model_dir: Directory containing the model checkpoint
            - model_name: Name of the checkpoint file
            - use_gpu: Whether to use GPU (default: True)
            - seed: Random seed for reproducibility
        
        Returns
        -------
        calc : OCPCalculator
            Initialized eSEN calculator
        """
        try:
            from fairchem.core import OCPCalculator
        except ImportError:
            raise ImportError(
                "fairchem package is required for eSEN calculator. "
                "Install it with: pip install fairchem-core"
            )
        
        # Merge default config with user overrides
        config = cls.DEFAULT_CONFIGS['esen'].copy()
        config.update(kwargs)
        
        # Build model checkpoint path
        model_dir = config.get('model_dir')
        model_name = config.get('model_name')
        model_ckpt = str(Path(model_dir) / model_name)
        
        # Check if model exists
        if not os.path.exists(model_ckpt):
            raise FileNotFoundError(
                f"eSEN model checkpoint not found: {model_ckpt}\n"
                f"Please download the model or update the path."
            )
        
        # Create calculator
        use_gpu = config.get('use_gpu', True)
        seed = config.get('seed', 0)
        
        logger.info(f"\n Creating eSEN calculator from: {model_ckpt}")
        logger.info(f" Using GPU: {use_gpu}")
        
        calc = OCPCalculator(
            checkpoint_path=model_ckpt, 
            cpu=not use_gpu,
            seed=seed
        )
        
        # Disable mixed precision for stability
        calc.trainer.scaler = None
        
        # Set directory if provided
        if directory is not None:
            calc.directory = directory
            os.makedirs(directory, exist_ok=True)
        
        return calc
    
    @classmethod
    def _create_mace_calculator(cls, directory=None, **kwargs):
        """Create MACE calculator.
        
        Args
        ----
        directory : str, optional
            Output directory
        **kwargs : dict
            Override default configurations:
            - model_name: Name of the MACE model
            - model_url_template: URL template for downloading model
            - use_gpu: Whether to use GPU (default: True)
            - enable_cueq: Enable charge equilibration (default: False)
        
        Returns
        -------
        calc : MACECalculator
            Initialized MACE calculator
        """
        try:
            from mace.calculators import mace_mp
            import torch
        except ImportError:
            raise ImportError(
                "mace-torch package is required for MACE calculator. "
                "Install it with: pip install mace-torch"
            )
        
        # Merge default config with user overrides
        config = cls.DEFAULT_CONFIGS['mace'].copy()
        config.update(kwargs)
        
        # Determine device
        use_gpu = config.get('use_gpu', True)
        device = "cuda" if (use_gpu and torch.cuda.is_available()) else "cpu"
        
        # Build model checkpoint URL/path
        model_name = config.get('model_name')
        
        # Check if model_name is a local path or needs to be downloaded
        if os.path.exists(model_name):
            checkpoint = model_name
            logger.info(f"\n Loading MACE model from local path: {checkpoint}")
        else:
            # Use URL template
            url_template = config.get('model_url_template')
            checkpoint = url_template.format(model_name=model_name)
            logger.info(f"\n Downloading MACE model: {model_name}")
            logger.info(f" From URL: {checkpoint}")
        
        logger.info(f" Using device: {device}")
        
        # Create calculator
        enable_cueq = config.get('enable_cueq', False)
        calc = mace_mp(
            model=checkpoint,
            device=device,
            enable_cueq=enable_cueq
        )
        
        # Set directory if provided
        if directory is not None:
            calc.directory = directory
            os.makedirs(directory, exist_ok=True)
        
        return calc
    
    @classmethod
    def create_mlips_calculator(cls, model_name, directory=None, **kwargs):
        """Legacy method name for backward compatibility."""
        return cls.create_calculator(model_name, directory=directory, **kwargs)
    
    @classmethod
    def update_model_config(cls, model_name, **kwargs):
        """Update default configuration for a specific model.
        Args
        ----
        model_name : str
            Name of the model to update
        **kwargs : dict
            Configuration parameters to update
        """
        model_name_lower = model_name.lower()
        if model_name_lower not in cls.DEFAULT_CONFIGS:
            raise ValueError(f"Unknown model: {model_name}")
        
        cls.DEFAULT_CONFIGS[model_name_lower].update(kwargs)
        logger.info(f"Updated default config for {model_name}: {kwargs}")


def get_mlips_optimizer(atoms, calc_type="esen", fmax=0.01, logfile=None, trajectory=None):
    """Get the appropriate optimizer for MLIPS calculations.
    
    Different MLIPS models may work better with different optimizers.
    This function selects the optimal optimizer based on the calculator type.
    
    Args
    ----
    atoms : ASE Atoms
        The structure to optimize
    calc_type : str
        Type of MLIPS calculator ("esen", "mace")
    fmax : float
        Force convergence criterion in eV/Å
    logfile : str, optional
        Path to log file
    trajectory : str, optional
        Path to trajectory file
    
    Returns
    -------
    optimizer : ASE optimizer object
        Configured optimizer instance
    """
    from ase.optimize import BFGS, LBFGS, FIRE
    
    calc_type_lower = calc_type.lower()
    
    # Select optimizer based on calculator type
    if calc_type_lower == "esen":
        # FIRE works well for eSEN
        optimizer = FIRE(atoms, logfile=logfile, trajectory=trajectory)
        logger.info(f" Using FIRE optimizer for {calc_type}")
    elif calc_type_lower == "mace":
        # LBFGS is efficient for MACE
        optimizer = LBFGS(atoms, logfile=logfile, trajectory=trajectory)
        logger.info(f" Using LBFGS optimizer for {calc_type}")
    else:
        # Default to BFGS
        optimizer = BFGS(atoms, logfile=logfile, trajectory=trajectory)
        logger.info(f" Using BFGS optimizer (default) for {calc_type}")
    
    return optimizer


def run_mlips_relaxation(calculator, atoms, calc_type="esen", fmax=0.01, max_steps=500,
                        directory=None, save_trajectory=True):
    """Run structure relaxation using MLIPS calculator.
    
    Args
    ----
    calculator : MLIPS calculator
        The MLIPS calculator instance
    atoms : ASE Atoms
        Structure to relax
    calc_type : str
        Type of calculator ("esen", "mace")
    fmax : float
        Force convergence criterion in eV/Å
    max_steps : int
        Maximum optimization steps
    directory : str, optional
        Output directory (uses calculator.directory if not specified)
    save_trajectory : bool
        Whether to save optimization trajectory
    
    Returns
    -------
    converged : bool
        Whether optimization converged
    """
    # Set directory if provided
    if directory is not None:
        calculator.directory = directory
    
    # Call the main function
    result = run_calculation(calculator, atoms, calc_type=calc_type, 
                           fmax=fmax, max_steps=max_steps)
    
    # Convert return value: 0 means success (converged), -1 means error (not converged)
    return result == 0


def write_vasprun_xml(atoms, forces, energy, stress=None, filename="vasprun.xml"):
    """Write a minimal vasprun.xml file for compatibility with VASP workflows.
    
    This is useful when using MLIPS calculators but need output in VASP format
    for downstream analysis tools.
    
    Args
    ----
    atoms : ASE Atoms
        The structure
    forces : np.ndarray
        Atomic forces in eV/Å
    energy : float
        Total energy in eV
    stress : np.ndarray, optional
        Stress tensor in eV/Å³
    filename : str
        Output filename
    """
    from ase.io import write
    
    # Create a copy to avoid modifying original
    atoms_copy = atoms.copy()
    
    # Set results
    atoms_copy.calc = None  # Clear any existing calculator
    
    # Create a simple SinglePointCalculator
    from ase.calculators.singlepoint import SinglePointCalculator
    
    if stress is not None:
        calc = SinglePointCalculator(
            atoms_copy,
            energy=energy,
            forces=forces,
            stress=stress
        )
    else:
        calc = SinglePointCalculator(
            atoms_copy,
            energy=energy,
            forces=forces
        )
    
    atoms_copy.calc = calc
    
    # Write to vasprun.xml format
    try:
        write(filename, atoms_copy, format='vasp-xml')
        logger.info(f" Written VASP-compatible output to {filename}")
    except Exception as e:
        logger.warning(f" Could not write vasprun.xml: {e}")


def check_mlips_requirements(model_name):
    """Check if required packages are installed for the specified MLIPS model.
    
    Args
    ----
    model_name : str
        Name of the MLIPS model
    
    Returns
    -------
    available : bool
        Whether the model can be used
    message : str
        Information or error message
    """
    model_name_lower = model_name.lower()
    
    if model_name_lower == "esen":
        try:
            import fairchem.core
            return True, f"eSEN is available (fairchem-core version: {fairchem.core.__version__})"
        except ImportError:
            return False, "eSEN requires fairchem-core. Install with: pip install fairchem-core"
    
    elif model_name_lower == "mace":
        try:
            import mace
            import torch
            cuda_available = torch.cuda.is_available()
            return True, f"MACE is available. CUDA: {cuda_available}"
        except ImportError:
            return False, "MACE requires mace-torch. Install with: pip install mace-torch"
    
    else:
        return False, f"Unknown MLIPS model: {model_name}"


# Backward compatibility: Keep the old function names as aliases
def create_mlips_calculator(model_name, directory=None, **kwargs):
    """Legacy function for creating MLIPS calculators.
    
    This is kept for backward compatibility. New code should use
    MLIPSCalculatorFactory.create_calculator() instead.
    """
    logger.warning(
        "create_mlips_calculator() is deprecated. "
        "Use MLIPSCalculatorFactory.create_calculator() instead."
    )
    return MLIPSCalculatorFactory.create_calculator(model_name, directory=directory, **kwargs)

def run_calculation(calc, atoms, 
                   calc_type="esen", fmax=0.01, max_steps=500):
    """ Run a calculation with either VASP or MLIPS
    
    Args
    -----
    calc : calculator object
        VASP calculator or MLIPS calculator (e.g., OCPCalculator, MACE)

    atoms : ASE Atoms object
        The structure to optimize
        
    calc_type : str
        Type of MLIPS calculator ("esen", "mace", etc.)
        
    fmax : float
        Force convergence criterion in eV/Å
        
    max_steps : int
        Maximum optimization steps
    
    """
    # MLIPS calculation using ASE
    atoms.calc = calc
    try:
        # Setup output directory and logging
        if hasattr(calc, 'directory') and calc.directory:
            import os
            os.makedirs(calc.directory, exist_ok=True)
            trajectory_file = f"{calc.directory}/opt.traj"
            logfile = f"{calc.directory}/opt.log"
        else:
            # Use current directory if no calculator directory specified
            trajectory_file = "opt.traj"
            logfile = "opt.log"
            
        # Test initial force calculation for debugging
        logger.info(f"Testing initial force calculation with {calc_type}...")
        try:
            initial_energy = atoms.get_potential_energy()
            initial_forces = atoms.get_forces()
            max_initial_force = np.max(np.linalg.norm(initial_forces, axis=1))
            logger.info(f"Initial energy: {initial_energy:.6f} eV")
            logger.info(f"Initial max force: {max_initial_force:.6f} eV/Å")
            
            if max_initial_force < 1e-5:
                logger.warning("⚠️  Initial forces are nearly zero - structure may already be optimized!")
                logger.info("Proceeding with optimization anyway...")
            
        except Exception as e:
            logger.error(f"❌ Error in initial force calculation: {e}")
            
        # Use the specialized optimizer selection function with proper ASE parameters
        opt = get_mlips_optimizer(atoms, calc_type=calc_type, fmax=fmax,
                                logfile=logfile, trajectory=trajectory_file)

        logger.info(f"Starting MLIPS optimization: fmax={fmax}, max_steps={max_steps}")

        # Run optimization - trajectory and logging handled automatically by optimizer
        opt.run(fmax=fmax, steps=max_steps)
        
        # Check final result
        final_forces = atoms.get_forces()
        max_final_force = np.max(np.linalg.norm(final_forces, axis=1))
        logger.info(f"Optimization completed. Final max force: {max_final_force:.6f} eV/Å")
        
        # Save final structure - determine output path
        if hasattr(calc, 'directory') and calc.directory:
            output_path = f"{calc.directory}/CONTCAR"
        else:
            output_path = "CONTCAR"
            
        from ase.io import write
        write(output_path, atoms, format='vasp')
        logger.info(f"MLIPS optimization completed. Final structure saved to {output_path}")
            
    except Exception as e:
        msg = f"\n Error in MLIPS calculation: {e}"
        logger.error(msg)
        return -1
        
    return 0

