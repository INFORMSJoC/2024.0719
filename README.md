[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Asymptotically Tight MILP Approximations for a Non-convex QCP

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported in the paper 
[Asymptotically Tight MILP Approximations for a Non-convex QCP](https://doi.org/10.1287/ijoc.2024.0719) by Shiyi Jiang, Jianqiang Cheng, Kai Pan, and Boshi Yang. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0719

https://doi.org/10.1287/ijoc.2024.0719.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Jiang2025,
  author =        {Jiang, Shiyi AND Cheng, Jianqiang AND Pan, Kai AND Yang, Boshi},
  publisher =     {INFORMS Journal on Computing},
  title =         {Asymptotically Tight MILP Approximations for a Non-convex QCP},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0719.cd},
  url =           {https://github.com/INFORMSJoC/2024.0719},
  note =          {Available for download at \url{https://doi.org/10.1287/ijoc/2024.0719.cd}, \url{https://github.com/INFORMSJoC/2024.0719}},
}  
```

## Description

The goal of this repository is to share data and results related to the non-convex quadratically constrained program (QCP), solved using our asymptotically tight mixed-integer linear programming (MILP) approximations.



## Building

We develop two asymptotically tight mixed-integer linear programming approximations for a non-convex quadratically constrained program. We demonstrate the effectiveness of our approaches via numerical experiments by applying our proposed approximations to randomly generated instances and two industrial applications.

## Prerequisites
- Python 3.8 or higher
- Required Libraries: Gurobi, NumPy, SciPy

## Structure

The source code is available in the [src](src) directory, which includes the following components:
- `ABprimaryProblem.py`: Defines the `qcqp_optimize` function, which solves the A,B decomposition of the original QCP problem. 
- `accuracy.py`: Defines the `print_info` function, which calculates the number of variables and constrains of different approximations given an approximation quality parameter $\epsilon$.
- `adaptiveRefineAlgo.py`: Defines the `qcqp_optimize` function, which uses an adaptive refined algorithm to find a feasible solution of the QCP problem.
- `approximate_sub.py`: Defines the `v_total_max` function, which calculates the upper bound of variable $v$.
- `bental_poly_appro_y.py`: Defines the `qcqp_optimize` function, which solves the first approximation where the SOC constrains are approximated using the first polyhedral approximation.
- `bental_poly_appro_z.py`: Defines the `qcqp_optimize` function, which solves the first approximation where the complement of SOC constrains are approximated.
- `bental_poly_appro.py`: Defines the `qcqp_optimize` function, which solves the first MILP approximation. 
- `bental_polyhedral_1.py`: Defines some functions, which are used to calculate coefficient matrices of the first polyhedral approximation. 
- `bental_polyhedral_2.py`: Defines some functions, which give another method to calculate coefficient matrices of the first polyhedral approximation. 
- `binary.py`: Gives a case study on the pure binary QCP problem. 
- `case1.py`, `case2.py`, `case3.py`, and `case4.py`: Give four cases on the small-scale QCP problems. 
- `directSolve.py`: Defines the `qcqp_optimize` function, which solves the original QCP problem. 
- `generateData.py`, `generateInteger.py`, and `generateNonconvex.py`: Define some functions, which are used to randomly generate data. 
- `linearApproximation.py`: Defines the `qcqp_optimize` function, which solves the second MILP approximation.
- `lineObj.py`: Defines the `qcqp_optimize` function, which solves the QCQP problems with linear objective. 
- `lowerboundLazy.py`: Defines the `qcqp_optimize` function, which solves the approximations with lazy constraints. 
- `sdpRelax.py`: Defines the function to solve the SDP relaxation of the original QCP problem. 
- `lpRelax.py`: Defines the function to solve the LP relaxation of the original QCP problem.  
- `socpRelax.py`: Defines the function to solve the SOCP relaxation of the original QCP problem. 
- `YZprimaryProblem.py`: Defines the `qcqp_optimize` function, which solves the Y,Z decomposition of the original QCP problem. 


## Data 

The instances are randomly generated according to some rules. Please see the [data](data) directory to view the data. 
This directory includes five folders: "Random," "SensitivityAnalyses," "JDE," "TTRS," and "Unitbox." 

The "Random" folder includes all randomly generated instances considered in Section 4.3. Each instance (denoted by "N_M1_M2_E_I") represents a QCP problem with N variables, M1 convex constraints, and M2 non-convex constraints. Note that E is the number of nonnegative eigenvalues in the matrix and I is the instance id. 

The "SensitivityAnalyses" folder includes all instances used for sensitivity analyses considered in Section 4.4. Each instance (denoted by "N_M1_M2_E_I") represents a QCP problem with N variables, M1 convex constraints, and M2 non-convex constraints. Note that E is the number of nonnegative eigenvalues in the matrix and I is the instance id. 

The "JDE" and "TTRS" folders include all instances of the joint decision and estimation (JDE) problem and the two-trust-region subproblem (TTRS) considered in Section 4.5, respectively. Each instance (denoted by "N_I") represents a corresponding problem with N variables. Note that I is the instance id.

The "Unitbox" folder includes all randomly generated unitbox instances considered in Appendix F. Each instance (denoted by "N_M_D_I") represents a unitbox QCP problem with N variables and M constraints. Note that D is the density (%) of the matrix and I is the instance id.


## Results

Please see the [results](results) directory to view the log information of Gurobi. 
This directory includes five folders: "Random," "SensitivityAnalyses," "JDE," "TTRS," and "Unitbox." 

The "Random" folder records the log information of solving randomly generated instances considered in Section 4.3. This folder includes six folders: "Approximations," "Approximations2," "Gurobi," "LP_Relax," "SDP_Relax," and "SOCP_Relax." In each folder, each file (named "N_M1_M2_E") records the log information of solving the corresponding approximation or original problem with N variables, M1 convex constraints, M2 non-convex constraints. Note that E is the number of nonnegative eigenvalues in the matrix. 


The "SensitivityAnalyses" folder records the log information of solving randomly generated instances used for sensitivity analyses considered in Section 4.4. This folder includes three folders: "Approximations," "Approximations2," and "Gurobi." In each folder, each file (named "N_M1_M2_E_A") records the log information of solving the corresponding approximation or original problem with N variables, M1 convex constraints, and M2 non-convex constraints. Note that E is the number of nonnegative eigenvalues in the matrix and A represents the parameter of the approximation quality ($\times 10^{-3}$), which is optional. 

The "JDE" and "TTRS" folders record the log information of solving the joint decision and estimation (JDE) problem and the two-trust-region subproblem (TTRS) considered in Section 4.5. Both of them include three folders: "Approximations," "Approximations2," and "Gurobi." In each folder, each file (named "N") records the log information of solving the corresponding approximation or original problem with N variables. 

The "Unitbox" folder records the log information of solving randomly generated unitbox instances considered in Appendix F. This folder includes six folders: "Approximations," "Approximations2," "Gurobi," "LP_Relax," "SDP_Relax," and "SOCP_Relax." In each folder, each file (named "N_M_D") records the log information of solving the corresponding approximation or original problem with N variables and M constraints. Note that D is the density (%) of the matrix.


## Run


### Individual Execution

For users who wish to run codes on a personal computer, we recommend using the following command to execute the program:

```bash
python directSolve.py <instance_path> <epsilon> >> <output_path>
```

Arguments:

- `<instance_path>`: Specifies the file path for the instance to be solved.
- `<epsilon>`: Specifies the approximation quality parameter.
- `<output_paht>`: Specifies the output file path.


### HPC Job Submission

For users who intend to execute codes on HPC servers, we recommend utilizing the following job script as a reference.

```bash
#!/bin/bash

#SBATCH --partition=h07q1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32gb
#SBATCH -J "benlb400397"
#SBATCH --output=./log/benlb400397.log



EXE_FILE="./bental_poly_appro.py"
FILE_NAME_PREFIX="100_3_2_97_"
EPSILON=0.001
OUTPUT_FOLDER="./QCQPData_results/Approximations/"


export SLURM_EXPORT_ENV=ALL
export GRB_LICENSE_FILE="./gurobi.lic"
module load Anaconda3/2023.03-1
source activate SHPC-env
pip3 install numpy gurobipy


for((i=0; i!=5; i++)) 
do
    python3 $EXE_FILE $FILE_NAME_PREFIX$i $EPSILON >> $OUTPUT_FOLDER$FILE_NAME_PREFIX
done

conda deactivate
module unload Anaconda3/2023.03-1 
```

## Support

For support in using the data and code, submit an [issue](https://github.com/INFORMSJoC/2024.0719/issues/new).
