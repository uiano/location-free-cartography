## location-free cartography
The MATLAB code in gsim folder implements the simulations and plots the figures described in the paper "Location-free Spectrum Cartography" by Yves Teganya,  Daniel Romero, Luis Miguel Lopez-Ramos, and Baltasar Beferull-Lozano.

##### Requirements: Matlab 9.1.0(2016b) or later
## Guidelines
The first time you want to run a simulation after downloading the code, execute
```gsimStartup```.
After that, you will be able to execute any simulation in the aforementioned paper by invoking the function in ```gsim.m```. that will have been created in the folder **gsim/**.

The experiments reproducing different figures are organized in functions located in the file ```LocFCartogrExperiments.m```.
For example, to run experiment 401 with 100 iterations, type ```gsim (0, 401, 100)```. To just display the results of the last execution of experiment 4001, type ```gsim(1, 401)```. 
## Citation
If you find our code helpful in your resarch or work, please cite our paper.
```
@article{teganya2019tlocation,
  title={Location-free Spectrum Cartography},
  author={Teganya, Yves and Romero, Daniel and Lopez-Ramos, Luis Miguel and Beferull-Lozano, Baltasar},
  journal={IEEE Transactions on Signal Processing},
  volume={xxxx},
  number={x},
  pages={xxxx--xxxx},
  year={2019},
  publisher={IEEE}
}
```
