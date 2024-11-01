import pandas as pd
from tqdm import tqdm
import mygene
import time

def get_string_to_hugo_mapping(string_ids):
    """
    Get mapping from STRING protein IDs to HUGO symbols using mygene
    
    Args:
        string_ids (set): Set of STRING protein IDs to map
        
    Returns:
        dict: Mapping from STRING protein IDs to HUGO symbols
    """
    print("Getting STRING ID to HUGO mapping using mygene...")
    
    # Initialize mygene client
    mg = mygene.MyGeneInfo()
    
    # Convert STRING IDs to Ensembl protein IDs
    ensembl_ids = set()
    string_to_ensembl = {}
    for string_id in string_ids:
        if string_id.startswith('9606.ENSP'):
            ensembl_id = string_id.replace('9606.ENSP', 'ENSP')
            ensembl_ids.add(ensembl_id)
            string_to_ensembl[string_id] = ensembl_id
    
    # Query mygene in batches to avoid timeout
    batch_size = 1000
    ensembl_list = list(ensembl_ids)
    mapping = {}
    
    print(f"Querying mygene for {len(ensembl_list)} proteins...")
    for i in tqdm(range(0, len(ensembl_list), batch_size)):
        batch = ensembl_list[i:i + batch_size]
        
        try:
            # Query mygene
            results = mg.querymany(batch, 
                                 scopes='ensembl.protein',
                                 fields='symbol',
                                 species='human',
                                 returnall=True)
            
            # Process results
            for hit in results['out']:
                if 'symbol' in hit:
                    ensembl_id = hit['query']
                    hugo = hit['symbol']
                    # Convert back to STRING ID format
                    string_id = [sid for sid, eid in string_to_ensembl.items() if eid == ensembl_id]
                    if string_id:
                        mapping[string_id[0]] = hugo
            
            # Small delay to avoid overwhelming the server
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Warning: Error in batch {i}-{i+batch_size}: {str(e)}")
            continue
    
    print(f"Successfully mapped {len(mapping)} proteins to HUGO symbols")
    return mapping

def process_string_network(input_file, output_file, min_score=400):
    """
    Process STRING network file:
    1. Filter by minimum score
    2. Convert protein IDs to HUGO symbols
    
    Args:
        input_file (str): Path to input STRING network file
        output_file (str): Path to output file
        min_score (int): Minimum combined score to keep (default: 400)
    """
    # First pass: collect all unique protein IDs
    print("Collecting unique protein IDs...")
    unique_proteins = set()
    total_lines = 0
    
    with open(input_file, 'r') as f:
        # Skip header line
        header = f.readline()
        
        for line in tqdm(f):
            try:
                protein1, protein2, score = line.strip().split()
                if int(score) >= min_score:
                    unique_proteins.add(protein1)
                    unique_proteins.add(protein2)
                total_lines += 1
            except ValueError:
                print(f"Warning: Skipping malformed line: {line.strip()}")
                continue
    
    print(f"\nFound {len(unique_proteins)} unique proteins in network")
    
    # Get mapping
    string_to_hugo = get_string_to_hugo_mapping(unique_proteins)
    
    print(f"\nProcessing network file (minimum score: {min_score})...")
    
    # Initialize counters
    total_processed = 0
    total_kept = 0
    
    # Open output file
    with open(output_file, 'w') as out_f:
        # Write header
        out_f.write("protein1_hugo\tprotein2_hugo\tcombined_score\n")
        
        # Process input file line by line
        with open(input_file, 'r') as in_f:
            # Skip header line
            next(in_f)
            
            for line in tqdm(in_f, total=total_lines):
                try:
                    total_processed += 1
                    
                    # Parse line
                    protein1, protein2, score = line.strip().split()
                    score = int(score)
                    
                    # Check score
                    if score >= min_score:
                        # Get HUGO symbols
                        hugo1 = string_to_hugo.get(protein1)
                        hugo2 = string_to_hugo.get(protein2)
                        
                        # Write if both HUGO symbols exist
                        if hugo1 and hugo2:
                            out_f.write(f"{hugo1}\t{hugo2}\t{score}\n")
                            total_kept += 1
                    
                    # Print progress every million lines
                    if total_processed % 1000000 == 0:
                        print(f"\nProcessed: {total_processed:,} pairs")
                        print(f"Kept: {total_kept:,} pairs")
                
                except ValueError:
                    print(f"Warning: Skipping malformed line: {line.strip()}")
                    continue
    
    print(f"\nFinal counts:")
    print(f"Total processed: {total_processed:,} pairs")
    print(f"Total kept: {total_kept:,} pairs")

def analyze_network_statistics(output_file):
    """
    Print basic statistics about the processed network
    
    Args:
        output_file (str): Path to processed network file
    """
    print("\nAnalyzing network statistics...")
    
    # Read the processed network
    df = pd.read_csv(output_file, sep='\t')
    
    # Calculate statistics
    total_interactions = len(df)
    unique_genes = set(df['protein1_hugo'].unique()) | set(df['protein2_hugo'].unique())
    avg_score = df['combined_score'].mean()
    
    print(f"\nNetwork Statistics:")
    print(f"Total interactions: {total_interactions:,}")
    print(f"Unique genes: {len(unique_genes):,}")
    print(f"Average combined score: {avg_score:.2f}")

def main():
    # File paths
    input_file = "9606.protein.links.v12.0.txt"  #
    output_file = "network.tsv"
    min_score = 400
    
    try:
        # Process network
        process_string_network(input_file, output_file, min_score)
        
        # Analyze results
        analyze_network_statistics(output_file)
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()
