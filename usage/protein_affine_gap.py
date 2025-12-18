import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.engine import AlignerEngine
from src.scoring import get_blosum62

def main():

    seq1 = "MKALWALLL"
    seq2 = "MKALWLLLV"
    
    engine = AlignerEngine(mode='global')
    
    score_dict = get_blosum62()
    
    M, I, D = engine.compute_affine_matrix(seq1, seq2, score_dict, gop=-11, gep=-1)
    
    final_score = max(M[-1, -1], I[-1, -1], D[-1, -1])
    
    print("--- Assignment 3: Protein Alignment with Affine Gaps ---")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Gap Opening Penalty (GOP): -11")
    print(f"Gap Extension Penalty (GEP): -1")
    print(f"\nFinal Affine Alignment Score: {final_score}")

if __name__ == "__main__":
    main()