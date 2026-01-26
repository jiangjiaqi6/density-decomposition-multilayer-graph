
import os
import argparse
import sys

def split_multilayer_graph(dataset_name):
    """
    Reads a raw dataset from './single_multilayer_file/{dataset_name}.txt'
    and splits it into './{dataset_name}/layer_x.txt' + mlg.conf
    """
    
    # 1. Define Paths
    # Get the directory where this script is located
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Input path: ./single_multilayer_file/DATASET.txt
    input_dir = os.path.join(base_dir, 'single_multilayer_file')
    input_file_path = os.path.join(input_dir, f"{dataset_name}.txt")
    
    # Output path: ./DATASET/
    output_dir = os.path.join(base_dir, dataset_name)

    # 2. Basic Validation
    if not os.path.exists(input_file_path):
        print(f"[Error] Input file not found: {input_file_path}")
        print(f"Please ensure '{dataset_name}.txt' exists inside 'single_multilayer_file' folder.")
        sys.exit(1)

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"[Info] Created output directory: {output_dir}")
    else:
        print(f"[Info] Output directory exists: {output_dir}")

    print(f"[Info] Processing dataset: {dataset_name}...")

    # 3. Processing Logic
    # We keep file handles open to avoid overhead of opening/closing files repeatedly
    layer_files = {} 
    
    try:
        with open(input_file_path, 'r') as f_in:
            # Read Header
            header = f_in.readline().strip()
            if not header:
                print("[Error] File is empty.")
                return
            
            # Optional: Parse header for info (not strictly needed for splitting, but good for check)
            try:
                meta = list(map(int, header.split()))
                print(f"[Info] Metadata -> Layers: {meta[0]}, Nodes: {meta[1]}, Edges: {meta[2]}")
            except ValueError:
                print("[Warning] Header format is not standard (integers expected). Proceeding anyway...")

            # Read Body
            line_count = 0
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) < 3:
                    continue

                # Parse: layer_id u v
                layer_id = parts[0]
                u = parts[1]
                v = parts[2]

                # If we haven't seen this layer yet, open a new file handle
                if layer_id not in layer_files:
                    filename = f"layer_{layer_id}.txt"
                    file_path = os.path.join(output_dir, filename)
                    layer_files[layer_id] = open(file_path, 'w')
                
                # Write to the specific layer file
                layer_files[layer_id].write(f"{u} {v}\n")
                line_count += 1
            
            print(f"[Success] Processed {line_count} edges.")

    except Exception as e:
        print(f"[Error] An error occurred during processing: {e}")
        return
    finally:
        # Close all open file handles
        for f in layer_files.values():
            f.close()

    # 4. Generate mlg.conf
    # Sort layer IDs numerically to ensure correct order in config file (0, 1, 2, 10...)
    try:
        sorted_layers = sorted(layer_files.keys(), key=int)
    except ValueError:
        sorted_layers = sorted(layer_files.keys()) # Fallback for non-numeric IDs

    conf_path = os.path.join(output_dir, "mlg.conf")
    with open(conf_path, 'w') as f_conf:
        for lid in sorted_layers:
            f_conf.write(f"layer_{lid}.txt\n")

    print(f"[Success] Configuration file created at: {conf_path}")
    print(f"[Done] Dataset '{dataset_name}' is ready.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a raw multilayer graph file into single-layer files.")
    
    parser.add_argument('dataset_name', type=str, 
                        help='The name of the dataset file (without .txt extensions), e.g., DBLP2')

    args = parser.parse_args()
    
    split_multilayer_graph(args.dataset_name)
