# Image-Denoising-with-MPI
Parallel processing by using C/C++ with MPI library.

```bash
mpic++ main.cpp -o mpi_project -std=c++11
mpiexec  -n 5 mpi_project input_output/lena200_noisy.txt input_output/out.txt 0.8 0.15
```

## Introduction
In this project, I have experienced parallel programming with C/C++ using MPI library. I have implemented a parallel algorithm for image denoising with the Ising model using Metropolis Hastings algorithm. 
## The Ising Model
The Ising model, named after the physicist Ernst Ising, is a mathematical model of ferromagnetism in statistical mechanics. The model consists of discrete variables that represent magnetic dipole moments of atomic spins that can be in one of two states (+1 or -1). The spins are arranged in a graph, usually a lattice, allowing each spin to interact with its neighbors.
- The Ising model models things that have two states (+1 or -1)
- Things interact with its neighbors
I assume that a black and white image is generated using The Ising model, then it means that if we take random black pixel from the image, it is more likely that this pixel is surrounded by black pixels (same for the white pixels).
## Parallel Image Denoising
In this project, my aim is to implement the parallel image denoising using MPI. I assume that the image is a square of size n×n. There will be 1 master and p slave processors and we will assume that n is divisible by the number of slave processors p.
In the second approach, the grid is divided into (n/p × n/p) blocks and each process is responsible from a block. Now the updates are more trickier since each process has more than 1 adjacent process.
## Implementation
1. The program reads the input from a text file and print the result in another text file. The input text file will be a 2D array representation of a black and white noisy image. The input file must only be read by master processor and distributed to slave (rest) processors by the master processor. The whole array should not be stored in each processor locally.
2. The program start by distributing the input among the processesors and let each processor work on its pixels without any communication. 
3. Any functioning of the program regarding the whole program such as printing the output are done by the master processor.

4. The names of the input and output files and the values of β and π priors will be given on the command line, as well as the number of processes. An example execution for a Windows user that runs on 4 processors and uses input.txt and output.txt files with β = 0.6 and π = 0.1 would be:
mpiexec -n 4 project.exe input.txt output.txt 0.6 0.1


- image_to_text.py: Converts your image into noise-free project ready text input file.
```bash
python_image_to text.py input image output file
```
- make_noise.py: Converts your noise-free project ready text input file into noisy project ready text input file. The code takes pi as a parameter that determines the noise rate. Try with your own images.
```bash
python make_noise.py input_file pi_value output_file
```

- text_to_image.py: Converts your project ready text file into image. So you can see your noisy image or the output of your code.
```bash
python text_to_image.py input_file output_image
```

A random algorithm you may not end up with the same outputs every time. Also more importantly, as my assumption of the real image coming from the Ising Model which is not really true, you will not end up with the exact same image with the original one, but a very close one.
