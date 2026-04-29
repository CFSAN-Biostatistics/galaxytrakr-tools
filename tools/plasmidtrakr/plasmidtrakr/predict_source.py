import argparse
import sys
import os
import pandas as pd
import joblib
import warnings

# Suppress scikit-learn warnings about feature names if they pop up
warnings.filterwarnings("ignore", category=UserWarning)

def main():
    parser = argparse.ArgumentParser(description="Predict the source of an isolate using a trained Random Forest model and Mash distances.")
    parser.add_argument("-i", "--input", required=True, help="Input Mash screen/dist file for one or more isolates.")
    
    # --- KEY FIX 1: Replaced -m and -f with a single -b (bundle) argument ---
    parser.add_argument("-b", "--bundle", required=True, help="Path to the bundled model and features (.joblib file)")
    
    parser.add_argument("-t", "--threshold", type=float, default=0.95, help="Mash identity threshold (default: 0.95)")
    parser.add_argument("-o", "--output", default="predictions.tsv", help="Output file for predictions (default: predictions.tsv)")
    args = parser.parse_args()

    print(f"Loading model bundle: {args.bundle}")
    
    try:
        # --- KEY FIX 2: Load the dictionary and extract both pieces ---
        bundle = joblib.load(args.bundle)
        rf_model = bundle['model']
        training_features = bundle['features']
        print(f"Successfully loaded model and {len(training_features)} features.")
    except Exception as e:
        print(f"FATAL: Error loading model bundle: {e}")
        sys.exit(1)

    print(f"Loading and processing input data: {args.input}")
    
    try:
        df = pd.read_csv(args.input, sep='\s+', header=None, engine='python')

        # Your format is from 'mash screen', where the columns are:
        # Identity, Shared-hashes, Median-multiplicity, P-value, Query-ID
        if len(df.columns) >= 5:
            print("--> Standard headerless Mash output detected.")
            # Keep only the first 5 columns to be safe
            df = df.iloc[:, :5]
            df.columns = ['Identity', 'Shared_Hashes', 'Median_Multiplicity', 'P_value', 'Plasmid_ID']
            
            # The 'Identity' is already the first column, just convert it to numeric
            df['Identity'] = pd.to_numeric(df['Identity'], errors='coerce')
            
            # We need to manually add the 'Run' column. For screen output, the Query-ID (isolate name)
            # is not present in the file itself. We must get it from the filename.
            run_id = os.path.splitext(os.path.basename(args.input))[0]
            df['Run'] = run_id
        else:
            print(f"FATAL: Input file format not recognized. Expected at least 5 columns for Mash output, but got {len(df.columns)}.")
            sys.exit(1)
            
        df.dropna(subset=['Identity'], inplace=True)
        df['Run'] = df['Run'].astype(str).str.strip()
        
    except Exception as e: 
        print(f"FATAL: Error reading input file '{args.input}'. Error: {e}")
        sys.exit(1)

    print(f"Filtering features (Identity >= {args.threshold})...")
    filtered_df = df[df['Identity'] >= args.threshold].copy()
    
    if filtered_df.empty:
        print("Warning: No plasmid hits met the identity threshold. Cannot make a prediction.")
        sys.exit(0)

    new_data_matrix = filtered_df.pivot_table(index='Run', columns='Plasmid_ID', values='Identity', fill_value=0)
    
    print("Aligning input features with the trained model...")
    aligned_matrix = pd.DataFrame(0, index=new_data_matrix.index, columns=training_features)
    common_plasmids = new_data_matrix.columns.intersection(training_features)
    aligned_matrix[common_plasmids] = new_data_matrix[common_plasmids]

    print(f"Making predictions for {len(aligned_matrix)} isolate(s)...")
    predictions = rf_model.predict(aligned_matrix)
    probabilities = rf_model.predict_proba(aligned_matrix)
    max_probs = probabilities.max(axis=1)

    results_df = pd.DataFrame({
        'Run': aligned_matrix.index,
        'Predicted_Source': predictions,
        'Confidence_Score': max_probs
    })

    results_df.to_csv(args.output, sep='\t', index=False)
    print(f"\n✅ Predictions complete! Saved to {args.output}")
    print("--- PREDICTION RESULTS ---")
    print(results_df.to_string(index=False))

if __name__ == "__main__":
    main()
