# Simulated Moving Bed Simulation with CADET

This repository contains example simulations of Simulated Moving Bed (SMB) processes using **CADET-Process** and **CADET-RDM**. A four-zone SMB system with eight columns and a five-zone SMB system with five columns is examined.

This example reproduces part of the case study from:

* *"Efficient numerical simulation of simulated moving bed chromatography with a single-column solver"*
  Qiao-Le He, Samuel Leweke, Eric von Lieres
  *Computers & Chemical Engineering* (2018); 111:183-198.
  [doi:10.1016/j.compchemeng.2017.12.022.](https://www.sciencedirect.com/science/article/pii/S0098135417304520)

---

## Authors

* Katharina Paul
* Ronald Jäpel
* Hannah Lanzrath
* Johannes Schmölder

---

## Running the Example Simulation

1. Clone this repository.
2. Set up the environment using the `environment.yml` file.
3. Run the simulation:

   ```bash
   python main.py
   ```

The results will be stored in the `src` folder inside the `output` directory.

> **Note**: Running `cadet-rdm` requires [**Git LFS**](https://git-lfs.com/), which needs to be installed separately.
>
> * **Ubuntu/Debian**:
>
>   ```bash
>   sudo apt-get install git-lfs
>   git lfs install
>   ```
>
> * **macOS** (with Homebrew):
>
>   ```bash
>   brew install git-lfs
>   git lfs install
>   ```
>
> * **Windows**:
>   Download and install from [https://git-lfs.com](https://git-lfs.com)

---

## Output Repository

The output data for this case study can be found here:
[https://github.com/cadet/RDM-Example-Simulated-Moving-Bed-Output](https://github.com/cadet/RDM-Example-Simulated-Moving-Bed-Output)