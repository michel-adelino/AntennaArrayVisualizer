import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda, beta_deg, element_type, currents, view="Vertical (XZ)"):
        """
        Calculates the normalized pattern for visualization.
        """
        n_points = 1000
        is_horizontal = "Horizontal" in view

        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        AF = np.zeros(n_points, dtype=complex)

        if is_horizontal:
            # Horizontal Plane (XY): Array along X logic for visualization
            phi_plot = np.linspace(0, 2 * np.pi, n_points)
            theta_for_ef_calc = np.full(n_points, np.pi / 2) # Theta is 90 deg on horizon
            psi = k * d_lambda * np.cos(phi_plot) + beta
            plot_angle = phi_plot
        else: 
            # Vertical Plane (XZ): Array along Z
            theta_plot = np.linspace(0, 2 * np.pi, n_points)
            theta_for_ef_calc = theta_plot
            psi = k * d_lambda * np.cos(theta_plot) + beta
            plot_angle = theta_plot
        
        # 1. Array Factor Summation
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        # 2. Element Factor
        EF = self._get_element_factor(theta_for_ef_calc, element_type)
        
        # 3. Total Field
        total_field = np.abs(AF) * EF
        
        # Normalize
        max_val = np.max(total_field)
        if max_val > 0:
            total_norm = total_field / max_val
        else:
            total_norm = total_field

        return plot_angle, total_norm

    def calculate_metrics(self, N, d_lambda, beta_deg, element_type, currents):
        """
        Calculates Directivity (D0) and HPBW based on the 3D physics of the array (Z-axis).
        """
        # High resolution theta for accurate integration
        theta = np.linspace(0, np.pi, 2000) 
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        # Calculate Field on the Principal Vertical Cut
        AF = np.zeros_like(theta, dtype=complex)
        psi = k * d_lambda * np.cos(theta) + beta
        
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        EF = self._get_element_factor(theta, element_type)
        E_total = np.abs(AF) * EF
        power_pattern = E_total**2
        max_power = np.max(power_pattern)
        
        if max_power == 0:
            return 0.0, 0.0

        # --- Directivity Calculation ---
        # Formula: D = 4*pi*Umax / Prad
        # For body of revolution (Z-axis array): Prad = 2*pi * Integral( |E|^2 * sin(theta) dtheta )
        # Therefore: D = 2 * max(|E|^2) / Integral( |E|^2 * sin(theta) )
        
        integrand = power_pattern * np.sin(theta)
        integral_val = np.trapz(integrand, theta) # Numerical integration
        
        if integral_val > 0:
            D_linear = 2 * max_power / integral_val
        else:
            D_linear = 0

        # --- HPBW Calculation ---
        # Find angular width where power >= 0.5 * max
        half_power = 0.5 * max_power
        idx_max = np.argmax(power_pattern)
        
        # Search left from peak
        idx_left = idx_max
        while idx_left > 0 and power_pattern[idx_left] > half_power:
            idx_left -= 1
            
        # Search right from peak
        idx_right = idx_max
        while idx_right < len(theta) - 1 and power_pattern[idx_right] > half_power:
            idx_right += 1
            
        hpbw_rad = theta[idx_right] - theta[idx_left]
        hpbw_deg = np.rad2deg(hpbw_rad)

        return D_linear, hpbw_deg

    def _get_element_factor(self, theta, el_type):
        if el_type == "Isotropic":
            return np.ones_like(theta)
        
        elif "Dipole" in el_type or "Monopole" in el_type:
            denominator = np.sin(theta)
            
            # Singularity handling at poles (theta = 0, pi)
            tol = 1e-5 
            ef = np.zeros_like(theta)
            valid_mask = np.abs(denominator) > tol
            
            numerator = np.cos((np.pi / 2) * np.cos(theta[valid_mask]))
            ef[valid_mask] = np.abs(numerator / denominator[valid_mask])
            
            if "Monopole" in el_type:
                # Zero out lower hemisphere
                mask_ground = (theta <= np.pi/2) | (theta >= 3*np.pi/2)
                ef = ef * mask_ground
                
            return ef
            
        return np.ones_like(theta)
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        af_safe = af_linear + 1e-12
        af_db = 20 * np.log10(af_safe)
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        return af_db_clipped