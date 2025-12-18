import numpy as np

class AlignerEngine:
    def __init__(self, mode='global', gap_penalty=-1):
        """
        Initialize the engine.
        mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman)
        gap_penalty: constant cost for a gap
        """
        self.mode = mode
        self.gap_penalty = gap_penalty

    def _initialize_matrix(self, n, m):
        """Creates and initializes the DP matrix based on the mode."""
        H = np.zeros((n + 1, m + 1))
        
        if self.mode == 'global':
            for i in range(1, n + 1):
                H[i, 0] = i * self.gap_penalty
            for j in range(1, m + 1):
                H[0, j] = j * self.gap_penalty
        
        return H

    def compute_matrix(self, seq1, seq2, score_dict):
        """
        Fills the scoring matrix.
        score_dict: the dictionary from scoring.py (BLOSUM62 or Nucleotide)
        """
        n, m = len(seq1), len(seq2)
        H = self._initialize_matrix(n, m)
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                char1, char2 = seq1[i-1], seq2[j-1]
                
                match = H[i-1, j-1] + score_dict.get((char1, char2), -4)
                
                deletion = H[i-1, j] + self.gap_penalty
                
                insertion = H[i, j-1] + self.gap_penalty
                
                if self.mode == 'global':
                    H[i, j] = max(match, deletion, insertion)
                else:
                    H[i, j] = max(0, match, deletion, insertion)
        
        return H

    def find_threshold_hits(self, H, threshold):
        """
        Returns coordinates of all cells with score >= threshold.
        """
        return np.argwhere(H >= threshold)