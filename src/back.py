import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda):
        """
        Calculates the normalized Array Factor (AF).
        
        Args:
            N (int): Number of antennas [cite: 11]
            d_lambda (float): Separation in wavelengths (d/lambda) [cite: 12]
        
        Returns:
            theta (array): Angles in radians
            af_norm (array): Normalized array factor
        """
        theta = np.linspace(0, 2 * np.pi, 360)
        
        k = 2 * np.pi
        beta = 0 # Broadside array [cite: 13]
        
        AF = np.zeros_like(theta, dtype=complex)
        
        for n in range(N):
            # Assuming uniform current I=1
            psi = k * d_lambda * np.cos(theta) + beta
            AF += 1 * np.exp(1j * n * psi)
            
        af_mag = np.abs(AF)
        
        # Avoid division by zero
        if np.max(af_mag) > 0:
            af_norm = af_mag / np.max(af_mag)
        else:
            af_norm = af_mag

        return theta, af_norm
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        """
        Converts the linear pattern to dB and applies a dynamic range.
        
        Args:
            af_linear (array): Normalized linear values (0 to 1).
            dynamic_range (float): Range in dB to display (e.g., 40).
        """
        # Avoid log(0) by adding a small epsilon
        af_safe = af_linear + 1e-12
        
        # Convert to dB
        af_db = 20 * np.log10(af_safe)
        
        # Apply the floor based on the selected dynamic range
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        
        return af_db_clipped