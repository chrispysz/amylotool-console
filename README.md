# AmyloTool: Console Application

This tool is intended for binary classification of potentionally amyloidogenic protein sequences from `FASTA` files. 

## Features
- Read protein sequences from a `.fa` file.
- Predict protein sequences with an optional prediction threshold.
- Generate CSV output with sequence ID, start index, end index, and prediction values.

## Build and Run using Docker

1. **Build the Docker image**:
   ```bash
   docker build -t image-name .
   ```
2. **Run the application**:
   You need to mount your local directory to access data from inside the Docker container. The following command mounts the `path/to/local/dir` directory and processes a file named `proteins.fa` located inside of it:
   ```bash
   docker run --name container-name -v path/to/local/dir:/data image-name /data/proteins.fa
   ```
   Example (Windows):
   ```bash
   docker run --name amylotool-console -v C:/Users/chris:/data amylotool-console /data/proteins.fa --threshold 0.5 --output_csv results.csv
   ```

### Arguments
`sequence_file` (mandatory): This is the path to the .fa file containing protein sequences you wish to predict.

`--threshold` (optional, *default: 0.0*): A float value ranging between 0 to 1. Values below it will not be saved in the final output.

`--output_csv` (optional, *default: output.csv*): Specify the name of the output CSV file where predictions will be written.

