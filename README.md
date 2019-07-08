## Introduction
The MATLAB code in gsim folder implements the simulations and plots the figures described in the paper "Location-free Spectrum Cartography" by Yves Teganya,  Daniel Romero, Luis Miguel Lopez-Ramos, and Baltasar Beferull-Lozano.

##### Requirements: Matlab 9.1.0(2016b) or later
## Guidelines
The first time you want to run a simulation after downloading the code, enter the folder gsim_locF and execute
```gsimStartup``` on the MATLAB prompt.
After that, you will be able to execute any simulation in the aforementioned paper by invoking the function in ```gsim.m```. that will have been created in the folder **gsim/**. 

The experiments reproducing different figures in the paper are organized in methods located in the file ```Experiments/LocFCartogrExperiments.m```. The comments before each method indicate which figure on the paper it generates.

You are now all set. For example, to run experiment 401 with 100 iterations, type ```gsim (0, 401, 100)```. To just display the results of the last execution of experiment 401 (stored in **Experiments/LocFCartogrExperiments_data**), type ```gsim(1, 401)```. 

For more information about the simulation evironment, please see <https://github.com/fachu000/GSim-base>.

## Citation
If you find our code helpful in your resarch or work, please cite our paper.
```
@article{teganya2019tlocation,
  title={Location-free Spectrum Cartography},
  author={Teganya, Yves and Romero, Daniel and Lopez-Ramos, Luis Miguel and Beferull-Lozano, Baltasar},
  journal={IEEE Transactions on Signal Processing},
  volume={67},
  number={15},
  pages={4013--4026},
  year={2019},
  publisher={IEEE}
}
```
