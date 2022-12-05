# Model Predictive Control of a Satellite

This repository contains the code developed in the context of the PhD thesis entitled "Hybrid and Distributed Control and Guidance for Cooperative On-Orbit Servicing".

The code uses MATLAB and ACADOS.

## Code architecture
```
satellite_mpc/
│
├── README.md
├── setup_acados.m
│
├── scripts/
│   ├── closed_loop_mpc_free_flyer_simple.m
│   └── ...
│
├── src/
│   ├── free_flyer_simple/
│   │    ├── free_flyer_simple.m
│   │    └── ...
│   └── ...
│
└── utils/
    ├── check_acados_requirements.m
    └── ...
```

The ```scripts/``` directory contains the scripts that run the MPC simulations. The ```src/``` directory contains the models written in the ACADOS sintax. The ```utils/``` directory contains functions that are used throughout the code.

## Downloading and running the code

Clone this repository to the parent directory of your ACADOS installation.

To run the code, open MATLAB and run the script ```setup_acados.m```. This sets up the ACADOS environment variables and adds the repository's directories and subdirectories to the MATLAB path.

Then, choose which script you wish to simulate from the ```scripts/```folder and run it in MATLAB.