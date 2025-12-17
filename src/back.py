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
        is_view_horizontal = "Azimuth" in view  # Are we looking at Azimuth?
        
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
                # Phase is constant (Omnidirectional Array Factor) when theta=90° for vertical arrays
                psi = beta 
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
        Calculates the exact directivity using numerical integration over the complete sphere.
        Valid for any orientation (Z or X) and element type.
        """
        # 1. Define a spherical grid (Theta: Elevation, Phi: Azimuth)
        # A resolution of 180x360 is sufficient for engineering precision (<0.1 dB error)
        n_theta = 180
        n_phi = 360
        theta = np.linspace(0, np.pi, n_theta)
        phi = np.linspace(0, 2 * np.pi, n_phi)
        
        # Create 2D meshes (Meshgrid)
        THETA, PHI = np.meshgrid(theta, phi, indexing='ij') # Shape: (180, 360)
        
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        # 2. Calculate PSI (Spatial phase) in 3D
        if "X" in array_axis:
            # Array on X: Phase depends on projection onto X: sin(theta)*cos(phi)
            psi_spatial = k * d_lambda * np.sin(THETA) * np.cos(PHI)
        else:
            # Array on Z: Phase depends on projection onto Z: cos(theta)
            psi_spatial = k * d_lambda * np.cos(THETA)
            
        psi = psi_spatial + beta

        # 3. Calculate Array Factor (AF) vectorized
        AF = np.zeros_like(THETA, dtype=complex)
        
        # Vectorized summation (Broadcast indices over the mesh)
        # Expand dims to multiply: currents[i] * exp(...)
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        # 4. Calculate Element Factor (EF)
        # We assume dipoles are always oriented in Z (Vertical),
        # so EF only depends on Theta, independent of array orientation.
        EF_1D = self._get_element_factor(theta, element_type)
        # Expand EF to 2D to match the mesh (copy along phi)
        EF = EF_1D[:, np.newaxis] * np.ones_like(PHI) # Explicit broadcasting

        # 5. Total Power U(theta, phi)
        total_field = np.abs(AF) * EF
        power = total_field**2
        max_p = np.max(power)
        
        if max_p == 0: return 0, 0

        # 6. Numerical Integration (Double Trapezoidal Integral)
        # Integral = integral( integral(U * sin(theta) dTheta) dPhi )
        integrand = power * np.sin(THETA)
        
        # Integrate first over Theta (axis=0)
        int_theta = np.trapz(integrand, theta, axis=0)
        # Integrate the result over Phi (axis=0 of the resulting vector)
        total_radiated_power = np.trapz(int_theta, phi)
        
        # 7. Directivity
        if total_radiated_power > 0:
            d_lin = 4 * np.pi * max_p / total_radiated_power
        else:
            d_lin = 0

        # 8. HPBW (Half Power Beam Width) for both planes
        hpbw_elevation = self._calculate_hpbw_from_cut(N, d_lambda, beta_deg, element_type, currents, array_axis, cut="elevation")
        hpbw_azimuth = self._calculate_hpbw_from_cut(N, d_lambda, beta_deg, element_type, currents, array_axis, cut="azimuth")

        return d_lin, hpbw_elevation, hpbw_azimuth

    def _calculate_hpbw_from_cut(self, N, d_lambda, beta_deg, element_type, currents, array_axis, cut="elevation"):
        """Auxiliary method to extract HPBW from the specified cut plane"""
        k = 2 * np.pi
        beta = np.deg2rad(beta_deg)
        indices = np.arange(N) - (N - 1) / 2
        
        # --- FIX 1: Siempre escanear 0 a 2pi (360°) para evitar cortes en los polos ---
        n_points = 2000
        var_angle = np.linspace(0, 2 * np.pi, n_points) 
        
        if cut == "elevation":
            if "X" in array_axis:
                psi = k * d_lambda * np.sin(var_angle) + beta
            else:
                psi = k * d_lambda * np.cos(var_angle) + beta
            EF = self._get_element_factor(var_angle, element_type)
            
        elif cut == "azimuth":
            if "X" in array_axis:
                psi = k * d_lambda * np.cos(var_angle) + beta
            else:
                psi = beta
            EF = np.ones_like(var_angle)
        else:
            raise ValueError("Invalid cut type")
        
        AF = np.zeros_like(var_angle, dtype=complex)
        for i, n in enumerate(indices):
            AF += currents[i] * np.exp(1j * n * psi)
            
        power = (np.abs(AF) * EF)**2
        max_p = np.max(power)
        if max_p == 0: return 0.0
        
        # --- FIX 2: Manejo de "Wrap-around" (Rotar el array para centrar el pico) ---
        # Esto soluciona problemas cuando el haz está en 0°, 180° o 360°
        
        idx_max = np.argmax(power)
        center_idx = n_points // 2
        shift = center_idx - idx_max
        
        # Rotamos el array de potencia para que el máximo quede en el centro
        power_shifted = np.roll(power, shift)
        
        half = 0.5 * max_p
        
        # Buscar hacia la izquierda desde el centro
        l = center_idx
        while l > 0 and power_shifted[l] > half:
            l -= 1
        
        # Buscar hacia la derecha desde el centro
        r = center_idx
        while r < n_points - 1 and power_shifted[r] > half:
            r += 1
            
        # Calcular ancho en índices y convertir a grados
        # Nota: Usamos 360 (no 2pi) porque el resultado final se espera en grados
        # El rango es 0 a 360, por lo que la resolución es 360 / (n_points - 1)
        angle_res = 360.0 / (n_points - 1)
        width_indices = r - l
        hpbw_deg = width_indices * angle_res
        
        if hpbw_deg >= 359: 
            return 360.0
        
        return hpbw_deg

    def _get_element_factor(self, theta, el_type):
        """Return element factor vs theta for given element type.

        NOTE: This implementation assumes element orientation is vertical (Z-axis).
        Therefore the element factor depends only on theta (elevation). If you
        need to support arbitrary element orientations, this method should be
        extended and documented accordingly.
        """
        if el_type == "Isotropic": return np.ones_like(theta)
        elif "Dipole" in el_type or "Monopole" in el_type:
            denominator = np.sin(theta)
            tol = 1e-5
            ef = np.zeros_like(theta)
            mask = np.abs(denominator) > tol
            num = np.cos((np.pi/2)*np.cos(theta[mask]))
            ef[mask] = np.abs(num/denominator[mask])
            if "Monopole" in el_type: ef = ef * ((theta <= np.pi/2) | (theta >= 3*np.pi/2))
            return ef
        return np.ones_like(theta)
    
    def convert_to_db(self, af_linear, dynamic_range=40):
        af_safe = af_linear + 1e-12
        af_db = 20 * np.log10(af_safe)
        floor_val = -abs(dynamic_range)
        af_db_clipped = np.clip(af_db, floor_val, 0)
        return af_db_clipped