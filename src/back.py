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