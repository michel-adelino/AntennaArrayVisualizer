import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda, beta_deg, element_type, currents):
        """
        Calculates the total pattern using the Multiplication Theorem.
        
        Args:
            N (int): Number of antennas
            d_lambda (float): Separation d/lambda
            beta_deg (float): Phase difference in degrees
            element_type (str): "Isotropic", "Dipole (lambda/2)", or "Monopole (lambda/4)"
            currents (array): Amplitude/Intensity for each element
        """
        # Theta from 0 to 2pi
        theta = np.linspace(0, 2 * np.pi, 1000)
        
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        
        # --- 1. ARRAY FACTOR (AF) ---
        # Use centered indices for symmetry: -(N-1)/2 to (N-1)/2
        indices = np.arange(N) - (N - 1) / 2
        
        AF = np.zeros_like(theta, dtype=complex)
        
        # psi = k * d * cos(theta) + beta
        # NOTE: theta=0 is the Z-axis (Up)
        psi = k * d_lambda * np.cos(theta) + beta
        
        # Summation with individual currents
        for i, n in enumerate(indices):
            # currents[i] is the magnitude for the n-th element
            AF += currents[i] * np.exp(1j * n * psi)
            
        # --- 2. ELEMENT FACTOR (EF) ---
        EF = self._get_element_factor(theta, element_type)
        
        # --- 3. TOTAL FIELD ---
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
        
        elif "Dipole" in el_type or "Monopole" in el_type:
            
            # 1. Calculate denominator
            denominator = np.sin(theta)
            
            # 2. FIX SINGULARITY (Avoid division by zero)
            tol = 1e-5 
            ef = np.zeros_like(theta)
            valid_mask = np.abs(denominator) > tol
            
            numerator = np.cos((np.pi / 2) * np.cos(theta[valid_mask]))
            ef[valid_mask] = np.abs(numerator / denominator[valid_mask])
            
            # --- MONOPOLE SPECIFIC LOGIC ---
            if "Monopole" in el_type:
                # Mask out the bottom hemisphere (90 to 270 deg)
                mask_ground = (theta <= np.pi/2) | (theta >= 3*np.pi/2)
                ef = ef * mask_ground
                
            return ef
            
        return np.ones_like(theta)
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        # Avoid log(0)
        af_safe = af_linear + 1e-12
        # Convert to dB
        af_db = 20 * np.log10(af_safe)
        # Apply floor
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        return af_db_clipped