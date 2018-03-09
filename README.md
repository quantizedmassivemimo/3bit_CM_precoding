# Simulator for 3-bit Constant-Modulus Precoding in VLSI with C3PO
(c) 2018 Christoph Studer, Oscar Castañeda, and Sven Jacobsson
e-mail: studer@cornell, oc66@cornell.edu, & sven.jacobsson@ericsson.com

More information about our research can be found at [http://vip.ece.cornell.edu] and [https://sites.google.com/site/durisi].

### Important information 

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein, and Christoph Studer, "VLSI Design of a 3-Bit Constant-Modulus Precoder for Massive MU-MIMO," IEEE International Symposium on Circuits and Systems (ISCAS), to appear in 2018

and clearly mention this in our paper.  

### How to start a simulation:

Simply run

```sh
precoder_sim
```

which starts a simulation in a 256 BS antennas, 16 users massive MIMO systems with 16-QAM modulation using ZF and MRT precoding (both infinite precision and 3-bit CM quantization), 2-bit (4-phase) C2PO, as well as the 3-bit C3PO algorithm introduced in the paper.

The simulator runs with predefined parameters. You can provide your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with other parameters, please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Version history
* Version 0.1: studer@cornell.edu - initial version for GitHub release
