import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.engine import AlignerEngine
from src.scoring import get_nucleotide_matrix

def main():

    seq1 = "ATCGTACGTA"
    seq2 = "ATCGGACGDA"
    
    engine = AlignerEngine(mode='global', gap_penalty=-2)
    
    score_dict = get_nucleotide_matrix(match=1, mismatch=-1)
    
    matrix = engine.compute_matrix(seq1, seq2, score_dict)
    
    print("--- Assignment 1: DNA Global Alignment ---")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print("\nResulting Scoring Matrix:")
    print(matrix)
    print(f"\nFinal Alignment Score: {matrix[-1, -1]}")

if __name__ == "__main__":
    main()