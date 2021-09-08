# Consensus Secondary Structure Predictor Engaging Structural Alignment-based Error Correction
# Installation
This project is written in Rust, a systems programming language.
You need to install Rust components, i.e., rustc (the Rust compiler), cargo (the Rust package manager), and the Rust standard library.
Visit [the Rust website](https://www.rust-lang.org) to see more about Rust.
You can install Rust components with the following one line:
```bash
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
[Rustup](https://github.com/rust-lang-nursery/rustup.rs) arranges the above installation and enables to switch a compiler in use easily.
As ConsAlifold's dependencies, you need to install [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) and [LocARNA-P (if you wish to use instead of ConsProb)](https://github.com/s-will/LocARNA).
You can install ConsAlifold as follows: 
```bash
$ # AVX, SSE, and MMX enabled for rustc (another example: RUSTFLAGS='--emit asm -C target-feature=+avx2 -C target-feature=+ssse3 -C target-feature=+mmx -C target-feature=+fma')
$ RUSTFLAGS='--emit asm -C target-feature=+avx -C target-feature=+ssse3 -C target-feature=+mmx' cargo install consalifold
```
Check if you have installed ConsAlifold properly as follows:
```bash
$ consalifold # Its available command options will be displayed.
```
By the following Python script, you can reproduce the figures shown in the paper describing ConsAlifold's principle:
```bash
$ cd scripts
$ # Please install python packages required to this reproduction.
$ # Saved figures will appear at the "../assets/images" directory.
$ ./run_all.py
```

# Author
[Heartsh](https://github.com/heartsh)

# License
Copyright (c) 2018 Heartsh  
Licensed under [the MIT license](http://opensource.org/licenses/MIT).
