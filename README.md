# Prediction Tool for RNA Consensus Secondary Structures to Consider RNA Pairwise Structural Alignments
# Installation
This project is written mainly in Rust, a systems programming language.
You need to install Rust components, i.e., rustc (the Rust compiler), cargo (the Rust package manager), and the Rust standard library.
Visit [the Rust website](https://www.rust-lang.org) to see more about Rust.
You can install Rust components with the following one line:
```bash
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
The above installation is done by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), and Rustup enables to easily switch a compiler in use.
As ConsAlifold's dependencies, you need to install [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) and [LocARNA-P (if you wish to use instead of ConsProb)](https://github.com/s-will/LocARNA).
You can install ConsAlifold as follows: 
```bash
$ RUSTFLAGS='--emit asm -C target-feature=+avx -C target-feature=+ssse3 -C target-feature=+mmx' cargo install consalifold # AVX, SSE, and MMX enabled for rustc (another example: RUSTFLAGS='--emit asm -C target-feature=+avx2 -C target-feature=+ssse3 -C target-feature=+mmx -C target-feature=+fma')
```
Check if you have installed ConsAlifold properly as follows:
```bash
$ consalifold # Its available command options will be displayed.
```
By the following Python script, you can reproduce the figures shown in the paper describing ConsAlifold's principle:
```bash
$ cd scripts
$ ./run_all.py # Please install python packages required to this reproduction. Saved figures will appear at the "../assets/images" directory.
```

# Docker Playground <img src="https://www.docker.com/sites/default/files/d8/styles/role_icon/public/2019-07/Moby-logo.png?itok=sYH_JEaJ" width="40">
Replaying computational experiments in academic papers is the first but troublesome step to understand developed computational methods.
I provide a Ubuntu-based computational environment implemented on [Docker](https://www.docker.com/) as a playground to try out ConsProb:
```bash
$ git clone https://github.com/heartsh/consprob && cd consprob
$ docker build -t heartsh/consprob .
```
You can dive into the Docker image "heartsh/consprob" built by the above commands, using Zsh:
```bash
$ docker run -it heartsh/consprob zsh
```

# Method Digest
[RNAalifold](https://www.tbi.univie.ac.at/RNA/) folds each RNA sequence alignment, minimizing the average free energy of a predicted consensus secondary structure.
Based on posterior column base-pairing probabilities on RNA consensus secondary structures, [PETfold](https://rth.dk/resources/petfold/) and [CentroidAlifold](https://github.com/satoken/centroid-rna-package) fold each RNA sequence alignment.
PETfold and CentroidAlifold correct potential errors in each input sequence alignment utilizing posterior nucleotide base-pairing probabilities on RNA secondary structures.
To achieve better alignment error correction than PETfold and CentroidAlifold, I developed ConsAlifold implemented in this repository.
ConsAlifold folds each RNA sequence alignment correcting its potential errors with average probabilistic consistency.

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).
