import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda, beta_deg, element_type):
        """
        Calculates the total pattern using the Multiplication Theorem.
        
        Args:
            N (int): Number of antennas
            d_lambda (float): Separation d/lambda
            beta_deg (float): Phase difference in degrees
            element_type (str): "Isotropic" or "Dipole (lambda/2)"
        """
        # Theta from 0 to 2pi
        theta = np.linspace(0, 2 * np.pi, 360)
        
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg) # Convert user input to radians
        
        # --- 1. ARRAY FACTOR (AF) ---
        # Use centered indices for symmetry: -(N-1)/2 to (N-1)/2
        # Example N=4: indices = [-1.5, -0.5, 0.5, 1.5]
        indices = np.arange(N) - (N - 1) / 2
        
        AF = np.zeros_like(theta, dtype=complex)
        
        # Vectorization of the calculation (faster than 'for' loop)
        # psi = k * d * cos(theta) + beta
        # Exponent = j * n * psi
        
        # NOTE: To match the convention where theta=0 is the Z-axis
        # and the array is along the Z-axis:
        psi = k * d_lambda * np.cos(theta) + beta
        
        for n in indices:
            AF += 1 * np.exp(1j * n * psi)
            
        # --- 2. ELEMENT FACTOR (EF) ---
        EF = self._get_element_factor(theta, element_type)
        
        # --- 3. MULTIPLICATION (Total Field) ---
        # Pattern = |AF| * |EF|
        total_field = np.abs(AF) * EF
        
        # Normalization
        if np.max(total_field) > 0:
            total_norm = total_field / np.max(total_field)
        else:
            total_norm = total_field

        return theta, total_norm

    def _get_element_factor(self, theta, el_type):
        """Calculates the normalized pattern of the individual element."""
        if el_type == "Isotropic":
            return np.ones_like(theta)
        
        elif "Dipole" in el_type: # Half-wave dipole
            # Assuming dipole is aligned with the Z-axis
            
            numerator = np.cos((np.pi / 2) * np.cos(theta))
            denominator = np.sin(theta)
            
            # Handle division by zero (at theta = 0 and theta = pi)
            # Using 'ignore' so it doesn't print warnings, then fix NaNs
            with np.errstate(divide='ignore', invalid='ignore'):
                ef = np.abs(numerator / denominator)
            
            # Replace NaNs or Infs with 0 (dipole nulls)
            ef = np.nan_to_num(ef, nan=0.0, posinf=0.0, neginf=0.0)
            return ef
            
        return np.ones_like(theta)
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        # Avoid log(0) by adding a small epsilon
        af_safe = af_linear + 1e-12
        # Convert to dB
        af_db = 20 * np.log10(af_safe)
        # Apply the floor based on the selected dynamic range
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        return af_db_clipped