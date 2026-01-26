# Density Decomposition of Multilayer Graphs

This project focuses on the density decomposition of multilayer graphs. 

## Dataset

Download datasets and put them into the `/dataset` folder. An example file `homo` is provided.

The missing datasets are available at [FirmCore](https://github.com/joint-em/FTCS/tree/main/Code/Datasets)  and [multilayer kCore](https://github.com/egalimberti/multilayer_core_decomposition).

### Data Preparation

The algorithm requires the multilayer graph to be decomposed into separate files for each layer, accompanied by a configuration file (`mlg.conf`). 

If your dataset is stored in a **single file** format, please use the provided preprocessing script to convert it.

###  Directory Structure
Before running the script, ensure your project directory is organized as follows:
```text
dataset/
├── split_graph.py              # The preprocessing script
└── single_multilayer_file/     # Place your raw .txt files here
    ├── DBLP2.txt
    └── ...
```

### Input Format
The input single file should follow this structure:
- **Header (1st line):** `#layers #nodes #edges`
- **Body (Subsequent lines):** `layer_id u v` (where `u` and `v` are node IDs)


**Example (`raw_graph.txt`):**
```text
3 100 200
0 1 2
0 2 3
1 1 5
2 4 6
...
```
###  Process the DBLP2 dataset
```
python split_graph.py DBLP2
```

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
./den_decom_vsht  ../dataset/homo