# Stochastic Completion Field
<p align="left">
	<img src="https://img.shields.io/badge/platform-ubuntu-blueviolet?style=for-the-badge"
			 alt="platform">
	<img src="https://img.shields.io/badge/language-haskell | stack-blueviolet?style=for-the-badge"
			 alt="language">
</p>
 Implementation of paper "Stochastic Completion Field: A Neural Model of Illusory Contour Shapes and Salience". Williams LR, Jacobs DW. 1997. <a href="https://www.cs.unm.edu/~williams/williams95stochastic.pdf">Link to paper</a>.

## Further Description
To be written

## Installation and build guide
This installation guide is for Linux/Ubuntu operating system. There are no plans for Windows and Apple though we don't believe there should be much problems with self-modifying the steps.
1. Install <a href="https://docs.haskellstack.org/en/stable/README/#how-to-install>">Stack</a>. (Stack will include Haskell). 
   - Stack is a tool to build and manage packages for Haskell programs.
   * For Ubuntu, this is as simple as running on the command line: 
     ```bash
     curl -sSL https://get.haskellstack.org/ | sh
     ```
     or
     ```bash
     wget -qO- https://get.haskellstack.org/ | sh
     ```
   * You may want to run the following command to install/update system dependencies on Ubuntu:
     ```bash
     sudo apt-get install g++ gcc libc6-dev libffi-dev libgmp-dev make xz-utils zlib1g-dev git gnupg netbase
     ```
   * <a href="https://docs.haskellstack.org/en/stable/install_and_upgrade/">Stack full installation guide</a>
   
2. Install <a href="http://www.fftw.org/">fftw3</a>. 
   - FFTW "is a C subroutine library for computing the discrete Fourier transform (DFT)".
   * Download and/or extract the latest stable release.
   * Open a new command line terminal and navigate to the folder with the downloaded files.
   * In succession, input the following commands:
     ```bash
     ./configure
     ```
     ```bash
     make
     ```
     ```bash
     sudo make install
     ```
   * <a href="http://www.fftw.org/fftw3_doc/Installation-and-Customization.html">FFTW3 full installation guide<a/>. Linux/Ubuntu is an Unix system.
	
3. Install <a href="http://www.netlib.org/blas/">BLAS</a> and <a href="http://www.netlib.org/lapack/">LAPACK<a/>. 
   - BLAS "are routines that provide standard building blocks for performing basic vector and matrix operations".
   - LAPACK "provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems".
   * For Ubuntu:
     ```bash
     sudo apt-get install libblas-dev liblapack-dev
     ```
	
4. Download/clone this repository. <a href="https://help.github.com/en/articles/cloning-a-repository">GitHub guide to cloning directory</a>.

5. On the command line, navigate to your local version of this repo.

6. Build program with 
   ```bash
   stack build
   ```
   
## Authors and Acknowledgement
Author: Xinhua Zhang. Department of Computer Science, University of New Mexico, US. <a href="https://github.com/XinhuaZhang">GitHub profile</a>.<br>
Supervisor: Lance R Williams, Prof. Department of Computer Science, University of New Mexico, US.

with documentation by Viet Than, Department of Electrical Engineering and Computer Science, Vanderbilt University, US. <a href="https://github.com/VietThan">GitHub</a>.

## License
<a href="LICENSE">MIT License</a>

## Project status
Work in progress.
