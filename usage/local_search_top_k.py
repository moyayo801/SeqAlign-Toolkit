import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.engine import AlignerEngine
from src.scoring import get_blosum62

def main():

    seq1 = "VSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"
    seq2 = "VSPADKTNVKAAWGKV" 
    
    engine = AlignerEngine(mode='local', gap_penalty=-4)
    
    score_dict = get_blosum62()
    H = engine.compute_matrix(seq1, seq2, score_dict)
    
    T = 20
    hits = engine.find_threshold_hits(H, T)
    
    print(f"--- Assignment 4: Local Alignment with Threshold T={T} ---")
    print(f"Max Score in Matrix: {np.max(H)}")
    print(f"Number of cells reaching threshold {T}: {len(hits)}")
    print("\nCoordinates of top hits (Row, Col):")
    for hit in hits[:10]: 
        print(f"Score {H[hit[0], hit[1]]} at position {hit}")

if __name__ == "__main__":
    main()