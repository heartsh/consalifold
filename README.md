# NeoAliFold Program, Program for Maximum-expected-accuracy Estimations of RNA Consensus Secondary Structures 
This project provides the NeoAliFold Program, a program for the maximum-expected-accuracy estimations of RNA consensus secondary structures.

# Installation
This project has been written in mainly Rust, a systems programming language.
So first, you need to install the Rust compiler, package manager, and standard library. 
Visit [the Rust website](https://www.rust-lang.org) to see more about this language.
You can install these 3 components with 1 line as follows:
```bash
$ curl https://sh.rustup.rs -sSf | sh
```
The above installation is done by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), so you can easily switch a compiler to use. 
Also you need to install the [Centroid package](https://github.com/satoken/centroid-rna-package) to predict pairing probabilities by the RNAalipfold algorithm.
Now you can install the PhyloAliFold program and its dependent, the PhyloProb program, as follows: 
```bash
$ cargo install phyloprob # You input the probabilities computed by this program to "phylofold"
$ cargo install phyloalifold
```
Check if this program has been installed properly as follows:
```bash
$ phyloprob
$ phyloalifold
```
After the test, the figures shown in the paper of the PhyloFold program can be reproduced:
```bash
$ cd src
$ ./run_all.py # Install python packages required to the reproduction. Saved figures will appear at the "../assets/images" directory.
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).
