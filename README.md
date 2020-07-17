# ConsAlifold, which Predicts Global Consensus RNA Secondary Structures to Consider Sparse Global Pairwise RNA Structural Alignments
# Installation
This project has been written in mainly Rust, a systems programming language.
You need to install the Rust components, which are rustc (the Rust compiler), cargo (the Rust package manager), and the Rust standard library.
Visit [the Rust website](https://www.rust-lang.org) to see more about the language.
You can install the components with one line as follows:
```bash
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
The above installation is done by [Rustup](https://github.com/rust-lang-nursery/rustup.rs), so you can easily switch a compiler to use. 
Also you need to install the [Centroid RNA package](https://github.com/satoken/centroid-rna-package) to predict posterior pairing probabilities on global consensus RNA secondary structure.
You can install ConsAlifold as follows: 
```bash
$ cargo install consalifold
```
Check if the program has been installed properly as follows:
```bash
$ consalifold # Its available command options will be displayed.
```
The figures shown in the paper of the program can be reproduced:
```bash
$ cd src
$ ./run_all.py # Install python packages required to the reproduction. Saved figures will appear at the "../assets/images" directory.
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).
