import numpy as np

class AntennaCalculator:
    def __init__(self):
        pass

    def calculate_pattern(self, N, d_lambda, beta_deg, element_type, currents, view="Vertical (XZ)", array_axis="Z"):
        """
        Calculates pattern based on Array Alignment (Z or X) and View (Vertical/Horizontal).
        """
        n_points = 1000
        # View identification
        is_view_horizontal = "Horizontal" in view  # Are we looking at Azimuth?
        
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        # --- LOGIC SELECTION ---
        if "Z" in array_axis: # (Vertical Array)
            if is_view_horizontal: # Azimuth Cut (XY)
                # Theta is fixed at 90 deg. Phi varies 0-360.
                plot_angle = np.linspace(0, 2*np.pi, n_points) # Phi
                theta_for_ef = np.full(n_points, np.pi/2)
                # Array on Z -> Depends on cos(theta). Theta is 90 -> cos(90)=0.
                # Phase is constant (Omnidirectional Array Factor)
                psi = k * d_lambda * 0 + beta 
            else: # Elevation Cut (XZ)
                # Theta varies 0-360. Phi fixed at 0.
                plot_angle = np.linspace(0, 2*np.pi, n_points) # Theta
                theta_for_ef = plot_angle
                # Array on Z -> Depends on cos(theta)
                psi = k * d_lambda * np.cos(plot_angle) + beta

        else: # array_axis == "X" (Horizontal Array)
            if is_view_horizontal: # Azimuth Cut (XY)
                # We scan Phi 0-360. Theta fixed at 90.
                plot_angle = np.linspace(0, 2*np.pi, n_points) # Phi
                theta_for_ef = np.full(n_points, np.pi/2)
                # Array on X -> Depends on sin(theta)*cos(phi).
                # sin(90)=1 -> Depends on cos(phi)
                psi = k * d_lambda * np.cos(plot_angle) + beta
            else: # Elevation Cut (XZ)
                # We scan Theta 0-360. Phi fixed at 0 (X-axis direction)
                plot_angle = np.linspace(0, 2*np.pi, n_points) # Theta
                theta_for_ef = plot_angle
                # Array on X -> Depends on sin(theta)*cos(0) = sin(theta)
                psi = k * d_lambda * np.sin(plot_angle) + beta

        # 1. Array Factor
        AF = np.zeros(n_points, dtype=complex)
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        # 2. Element Factor
        EF = self._get_element_factor(theta_for_ef, element_type)
        
        # 3. Total
        total = np.abs(AF) * EF
        max_val = np.max(total)
        norm = total / max_val if max_val > 0 else total
        
        return plot_angle, norm

    def calculate_metrics(self, N, d_lambda, beta_deg, element_type, currents, array_axis="Z"):
        """
        Metrics calculation adapted for array orientation.
        """
        theta = np.linspace(0, np.pi, 2000) # Integration variable
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        # Determine the "Principal Cut" variable for integration based on axis
        if "X" in array_axis:
            # For X-array, the "Principal Cut" is the Azimuth (XY) plane (Phi variation).
            # It behaves mathematically identical to the Z-array Elevation cut.
            var_angle = theta # Treating as 'phi' effectively for calc
            psi = k * d_lambda * np.cos(var_angle) + beta
            # Element factor in Azimuth (theta=90) for Vertical Dipole is 1.0 (Omni).
            EF = np.ones_like(var_angle) 
        else:
            # Z-array (Standard)
            var_angle = theta
            psi = k * d_lambda * np.cos(var_angle) + beta
            EF = self._get_element_factor(var_angle, element_type)

        AF = np.zeros_like(var_angle, dtype=complex)
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        total = np.abs(AF) * EF
        power = total**2
        max_p = np.max(power)
        
        if max_p == 0: return 0, 0
        
        # HPBW
        half = 0.5 * max_p
        idx_max = np.argmax(power)
        l, r = idx_max, idx_max
        while l > 0 and power[l] > half: l -= 1
        while r < len(power)-1 and power[r] > half: r += 1
        hpbw = np.rad2deg(var_angle[r] - var_angle[l])
        
        # Directivity Estimation
        if "X" in array_axis:
            # For horizontal array (X-axis) of vertical dipoles, the radiation pattern is not a body of revolution.
            # The directivity estimation here is based on integrating the principal cut (azimuth plane) only.
            # This is a rough engineering approximation and does not account for the full 3D pattern.
            # For more accurate directivity, a full 3D integration over all solid angles is required.
            # See e.g. Balanis, "Antenna Theory", 4th Ed., Section 2.5, for directivity definitions.
            denom = np.trapz(power, var_angle) # Simple 1D integral over principal cut
            d_lin = 2 * max_p * np.pi / denom if denom > 0 else 0 # Approximate directivity (linear array, principal cut only)
        else:
            # Standard Body of Revolution integral
            integrand = power * np.sin(var_angle)
            denom = 2 * np.pi * np.trapz(integrand, var_angle)
            d_lin = 4 * np.pi * max_p / denom if denom > 0 else 0

        d_dbi = 10*np.log10(d_lin) if d_lin > 0 else 0
        return d_lin, hpbw

    def _get_element_factor(self, theta, el_type):
        if el_type == "Isotropic": return np.ones_like(theta)
        elif "Dipole" in el_type or "Monopole" in el_type:
            denominator = np.sin(theta)
            tol = 1e-5
            ef = np.zeros_like(theta)
            mask = np.abs(denominator) > tol
            num = np.cos((np.pi/2)*np.cos(theta[mask]))
            ef[mask] = np.abs(num/denominator[mask])
            if "Monopole" in el_type: ef = ef * (theta <= np.pi/2)
            return ef
        return np.ones_like(theta)
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        af_safe = af_linear + 1e-12
        af_db = 20 * np.log10(af_safe)
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        return af_db_clipped