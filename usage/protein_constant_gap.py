import sys
import os

# Ensure src is in the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.engine import AlignerEngine
from src.scoring import get_blosum62

def main():

    seq1 = "MKALWALLL"
    seq2 = "MKALWLLLV"
    
    engine = AlignerEngine(mode='global', gap_penalty=-4)
    
    score_dict = get_blosum62()
    
    matrix = engine.compute_matrix(seq1, seq2, score_dict)
    
    print("--- Assignment 2: Protein Global Alignment (BLOSUM62) ---")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nFinal alignment score using BLOSUM62: {matrix[-1, -1]}")

if __name__ == "__main__":
    main()