# Density Decomposition of Multilayer Graphs

This project focuses on the density decomposition of multilayer graphs. 

## Dataset

Download datasets and put them into the `/dataset` folder. An example file `homo` is provided.

## Usage

#### First Step

1. Create a build directory and navigate into it:
    ```bash
    mkdir build
    cd build
    ```

2. Compile the project using CMake:
    ```bash
    cmake ..
    make
    ```

### Choice of Algorithm

You can select different algorithms by modifying the `CMakefile`:

- To run the `vsht` algorithm, add the following line to the `CMakefile`:
    ```cmakefile
    project(den_decom_vsht VERSION 1.0)
    ```

- To run the `hpdd` algorithm, add the following line to the `CMakefile`:
    ```cmakefile
    project(den_decom_hpdd VERSION 1.0)
    ```

- To run the `hpdd_plus` algorithm, add the following line to the `CMakefile`:
    ```cmakefile
    project(den_decom_hpdd_plus VERSION 1.0)
    ```

After modifying the `CMakefile`, recompile the project by running `make` again.

### Running the Project

Once compiled, you can run the project using the following command:

```bash
./den_decom_vsht