# HPC Performance Metrics Tools

This repository contains scripts and tools to estimate code running time, speedup, and efficiency on High-Performance Computing (HPC) systems. 

## Contents

- `performance_metrics.py`: Python script to calculate and visualize SCF time, speedup, and efficiency across various CPU configurations.
- `omp_parallelism.py`: Python script to analyze the impact of OpenMP (OMP) parallelism on calculation time.

## Requirements

- Python 3.x
- pandas
- matplotlib
- numpy

## Usage

### `performance_metrics.py`

This script calculates and visualizes SCF time, speedup, and efficiency based on CPU configurations specified in a CSV file.

#### Instructions:

1. Ensure you have Python 3.x installed on your system.
2. Install the required dependencies using `pip install -r requirements.txt`.
3. Place your data in a CSV file (e.g., `222_111kp_CPU_scf.csv`) organized with columns: `Number of CPUs`, `SCF Time`, and `Number of Nodes`.
4. Run the script using `python performance_metrics.py`.

### `omp_parallelism.py`

This script analyzes the impact of OpenMP parallelism on calculation time and visualizes the results in a bar plot.

#### Instructions:

1. Ensure you have Python 3.x installed on your system.
2. Install the required dependencies using `pip install -r requirements.txt`.
3. Run the script using `python omp_parallelism.py`.

## Contributors

- [Marco Fronzi](https://github.com/marcofronzi)

## License

This project is licensed under the MIT License
