# GCMC-par

Work-in-progress project implementing multi-thread Grand Canonical Monte Carlo simulations powered by a neural network potential trained on Tensor Flow via Panna [1]. This is a partial implementation of the serial code, submitted as a final project for CS 267 (Applications of Parallel Computers) in Spring 2023. The current version is applied to a MOF system, with CO$`_2`$ as a guest molecule. The implementation so far loads a MOF unit cell in VASP format, computes the neighbor list including ghost particles due to periodic boundary conditions, loads neural network parameters, computes Behler-Parinello descriptors for the local chemical environments upon loaded parameters, and executes the neural network to calculate the potential energy prediction for the system. 


[1] Pellegrini, Franco, et al. "PANNA 2.0: Efficient neural network interatomic potentials and new architectures." arXiv preprint arXiv:2305.11805 (2023).
