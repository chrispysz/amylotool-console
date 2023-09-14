import numpy as np
import tensorflow as tf
import argparse
import csv
import time

model = tf.keras.models.load_model('./models/ProteinBERT')

ALL_AAS = 'ACDEFGHIKLMNPQRSTUVWXY'
ADDITIONAL_TOKENS = ['<OTHER>', '<START>', '<END>', '<PAD>']

ADDED_TOKENS_PER_SEQ = 2

last_displayed_percent = 0

n_aas = len(ALL_AAS)
aa_to_token_index = {aa: i for i, aa in enumerate(ALL_AAS)}
additional_token_to_index = {token: i + n_aas for i, token in enumerate(ADDITIONAL_TOKENS)}
token_to_index = {**aa_to_token_index, **additional_token_to_index}
index_to_token = {index: token for token, index in token_to_index.items()}
n_tokens = len(token_to_index)

def tokenize_seq(seq):
    other_token_index = additional_token_to_index['<OTHER>']
    return [additional_token_to_index['<START>']] + [aa_to_token_index.get(aa, other_token_index) for aa in parse_seq(seq)] + [additional_token_to_index['<END>']]
            
def parse_seq(seq):
    if isinstance(seq, str):
        return seq
    elif isinstance(seq, bytes):
        return seq.decode('utf8')
    else:
        raise TypeError('Unexpected sequence type: %s' % type(seq))

def tokenize_seqs(seqs, seq_len):
    return np.array([seq_tokens + (seq_len - len(seq_tokens)) * [additional_token_to_index['<PAD>']] for seq_tokens in map(tokenize_seq, seqs)], dtype=np.int32)

def predict_window(seq, seq_cutoff=39, threshold=0.0):
    seqqs = [seq[i:seq_cutoff+i+1] for i in range(len(seq)-seq_cutoff)]
    preds = model.predict([tokenize_seqs(seqqs, 42), np.zeros((len(seqqs), 8943), dtype=np.int8)], verbose=0)
    seq_dicts = [{"startIndex": i, "endIndex": seq_cutoff+i, "prediction": float(pred[0])} for i, pred in enumerate(preds) if float(pred[0]) >= threshold]
    return seq_dicts

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        sequences = []
        sequence = ''
        seq_id = None
        for line in f:
            if line.startswith('>'):
                if sequence:
                    sequences.append((seq_id, sequence))
                    sequence = ''
                seq_id = line.strip()
            else:
                sequence += line.strip()
        if sequence:
            sequences.append((seq_id, sequence))
    return sequences


def seconds_to_hms(seconds):
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return hours, minutes, seconds

def display_progress(current, total, start_time):
    global last_displayed_percent
    percent_done = (current / total) * 100
    elapsed_time = time.time() - start_time

    estimated_total_time = elapsed_time * (total / current) if current != 0 else 0
    estimated_time_left = estimated_total_time - elapsed_time

    hours, minutes, seconds = seconds_to_hms(estimated_time_left)

    if percent_done - last_displayed_percent >= 2.5:
        print('\rProcessed: %d/%d sequences (%.2f%%) - Estimated time left: %d hours, %d minutes, %.2f seconds' % 
              (current, total, percent_done, hours, minutes, seconds))
        last_displayed_percent = percent_done

def main():
    global last_displayed_percent
    last_displayed_percent = 0

    parser = argparse.ArgumentParser(description="Predict protein sequences")
    parser.add_argument("sequence_file", type=str, help="The .fa file containing protein sequences to predict")
    parser.add_argument("--threshold", type=float, default=0.0, help="Prediction threshold from 0 to 1.")
    parser.add_argument("--output_csv", type=str, default="output.csv", help="Name of the output CSV file.")
    args = parser.parse_args()

    print("Reading sequences from %s" % args.sequence_file)

    sequences = read_fasta(args.sequence_file)

    with open("/data/" + args.output_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Sequence ID", "Start index", "End index", "Prediction"])

        total_sequences = len(sequences)
        processed_sequences = 0

        print("Predicting sequences...")

        start_time = time.time()  # Start the timer

        for seq_id, sequence in sequences:
            if len(sequence) < 40:
                continue

            results = predict_window(sequence, threshold=args.threshold)

            for result in results:
                writer.writerow([seq_id, result['startIndex'], result['endIndex'], result['prediction']])

            processed_sequences += 1
            display_progress(processed_sequences, total_sequences, start_time)
        
        print("Done!")

if __name__ == '__main__':
    main()
