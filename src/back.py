import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda, beta_deg, element_type, currents, view="Vertical (XZ)"):
        """
        Calculates the total pattern using the Multiplication Theorem.
        
        Args:
            N (int): Number of antennas
            d_lambda (float): Separation d/lambda
            beta_deg (float): Phase difference in degrees
            element_type (str): "Isotropic", "Dipole (lambda/2)", or "Monopole (lambda/4)"
            currents (array): Amplitude/Intensity for each element
            view (str): The plane to calculate, "Vertical (XZ)" or "Horizontal (XY)"
        """
        n_points = 1000
        is_horizontal = "Horizontal" in view

        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        AF = np.zeros(n_points, dtype=complex)

        if is_horizontal:
            # Horizontal Plane (XY): Assume array is along the X-axis
            # The plot variable is phi (azimuth), theta is fixed at 90 deg.
            phi_plot = np.linspace(0, 2 * np.pi, n_points)
            theta_for_ef_calc = np.full(n_points, np.pi / 2)
            
            # Phase term depends on phi
            psi = k * d_lambda * np.sin(phi_plot) + beta
            theta_plot = phi_plot

        else: # Vertical Plane (XZ): Assume array is along the Z-axis
            # The plot variable is theta (elevation), phi is fixed at 0 deg.
            theta_plot = np.linspace(0, 2 * np.pi, n_points)
            theta_for_ef_calc = theta_plot
            
            # Phase term depends on theta
            psi = k * d_lambda * np.cos(theta_plot) + beta
        
        # --- 1. ARRAY FACTOR (AF) ---
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        # --- 2. ELEMENT FACTOR (EF) ---
        EF = self._get_element_factor(theta_for_ef_calc, element_type)
        
        # --- 3. TOTAL FIELD ---
        total_field = np.abs(AF) * EF
        
        # Normalization
        if np.max(total_field) > 0:
            total_norm = total_field / np.max(total_field)
        else:
            total_norm = total_field

        return theta_plot, total_norm

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

    def calculate_directivity(self, pattern, theta):
        integral = np.trapezoid(pattern**2, theta)
        D = 2 * np.pi / integral
        return D

    def calculate_hpbw(self, pattern, theta):
        max_val = np.max(pattern)
        half_power = max_val / np.sqrt(2)
        peak_idx = np.argmax(pattern)
        peak_theta = theta[peak_idx]
        
        left_idx = peak_idx
        while left_idx > 0 and pattern[left_idx] > half_power:
            left_idx -= 1
        if left_idx < peak_idx:
            if pattern[left_idx] < half_power < pattern[left_idx+1]:
                frac = (half_power - pattern[left_idx]) / (pattern[left_idx+1] - pattern[left_idx])
                left_theta = theta[left_idx] + frac * (theta[left_idx+1] - theta[left_idx])
            else:
                left_theta = theta[left_idx]
        else:
            left_theta = theta[0]
        
        right_idx = peak_idx
        while right_idx < len(pattern)-1 and pattern[right_idx] > half_power:
            right_idx += 1
        if right_idx > peak_idx:
            if pattern[right_idx-1] > half_power > pattern[right_idx]:
                frac = (pattern[right_idx-1] - half_power) / (pattern[right_idx-1] - pattern[right_idx])
                right_theta = theta[right_idx-1] + frac * (theta[right_idx] - theta[right_idx-1])
            else:
                right_theta = theta[right_idx]
        else:
            right_theta = theta[-1]
        
        hpbw_rad = right_theta - left_theta
        hpbw_deg = np.rad2deg(hpbw_rad)
        return hpbw_deg