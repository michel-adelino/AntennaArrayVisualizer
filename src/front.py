import customtkinter as ctk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import tkinter as tk
from src.back import AntennaCalculator
import os
import sys
from mpl_toolkits.mplot3d import Axes3D

# --- CUSTOM TOOLTIP CLASS ---
class CTkTooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.widget.bind("<Enter>", self.show_tooltip)
        self.widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event=None):
        if self.tooltip_window or not self.text: return
        x = self.widget.winfo_rootx() + 25
        y = self.widget.winfo_rooty() + 25
        self.tooltip_window = tk.Toplevel(self.widget)
        self.tooltip_window.wm_overrideredirect(True)
        self.tooltip_window.wm_geometry(f"+{x}+{y}")
        label = tk.Label(self.tooltip_window, text=self.text, justify='left',
                         background="#2b2b2b", fg="#ffffff", relief='solid', borderwidth=1,
                         font=("Arial", 10, "normal"), padx=5, pady=2)
        label.pack(ipadx=1)

    def hide_tooltip(self, event=None):
        if self.tooltip_window:
            self.tooltip_window.destroy()
            self.tooltip_window = None

# --- MAIN APPLICATION ---
class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Antenna Array Visualizer - ITBA 22.21")
        self.geometry("1100x750")
        ctk.set_appearance_mode("Dark")
        ctk.set_default_color_theme("blue")

        # Set window icon
        icon_path = os.path.join(sys._MEIPASS, 'docs', 'icon.ico') if hasattr(sys, '_MEIPASS') else 'docs/icon.ico'
        try:
            self.iconbitmap(icon_path)
        except Exception:
            pass  # Fallback if icon not found

        self.calculator = AntennaCalculator()
        self.update_timer = None

        # State variables for cursor
        self.D_linear = 1.0
        self.theta_plot = None
        self.af_db_plot = None
        self.theta_deg_sorted = None
        self.af_db_sorted = None
        self.fixed_cursor = False
        self.fixed_x = None
        self.fixed_db = None
        self.last_x = None
        self.last_db = None
        self.cursor_point = None
        self.cursor_line = None
        self.cursor_text = ""

        # --- LAYOUT ---
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.frame_controls = ctk.CTkFrame(self, width=280)
        self.frame_controls.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.frame_controls.grid_columnconfigure(0, weight=1) 
        self.frame_controls.grid_columnconfigure(1, weight=2) 

        # --- PARAMETERS ---
        r = 0 
        self.lbl_title = ctk.CTkLabel(self.frame_controls, text="Parameters", font=("Arial", 20, "bold"))
        self.lbl_title.grid(row=r, column=0, columnspan=2, pady=(15, 10))
        r += 1
        
        self.lbl_hint = ctk.CTkLabel(self.frame_controls, text="(Hover labels for help)", font=("Arial", 11), text_color="gray")
        self.lbl_hint.grid(row=r, column=0, columnspan=2, pady=(0, 10))
        r += 1

        # Inputs
        self._add_input(r, "N Antennas:", "entry_n", "4", "Total number of radiating elements\nin the array.")
        r += 1
        self._add_input(r, "Separation (d/λ):", "entry_d", "0.5", "Distance between elements relative\nto wavelength (e.g., 0.5 = λ/2).")
        r += 1
        self._add_input(r, "Phase (β°):", "entry_beta", "0", "Progressive phase shift between\nconsecutive elements (degrees).\nUsed for beam steering.")
        r += 1
        self._add_input(r, "Intensities:", "entry_currents", "1", "Current intensity for each element.\nExample for N=4: '1, 0.5, 0.5, 0.33'\nSingle value applies to all.")
        r += 1

        # Combos
        self.lbl_type = ctk.CTkLabel(self.frame_controls, text="Type:", cursor="hand2")
        self.lbl_type.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_type, "Radiation pattern of the individual element.\n(Multiplied by Array Factor).")
        self.combo_type = ctk.CTkComboBox(self.frame_controls, values=["Isotropic", "Dipole (λ/2)", "Monopole (λ/4)"], command=lambda v: self.update_plot())
        self.combo_type.set("Isotropic")
        self.combo_type.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        self.lbl_axis = ctk.CTkLabel(self.frame_controls, text="Array Axis:", cursor="hand2")
        self.lbl_axis.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_axis, "Axis along which the array is oriented.\nZ: Vertical (elevation beamforming)\nX: Horizontal (azimuth beamforming)")
        self.combo_axis = ctk.CTkComboBox(self.frame_controls, values=["Z (Vertical)", "X (Horizontal)"], command=lambda v: self.update_plot())
        self.combo_axis.set("Z (Vertical)")
        self.combo_axis.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        self.lbl_view = ctk.CTkLabel(self.frame_controls, text="View:", cursor="hand2")
        self.lbl_view.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_view, "Cut Plane:\nVertical (Elevation/Theta)\nHorizontal (Azimuth/Phi)\nBoth: Show both views")
        self.combo_view = ctk.CTkComboBox(self.frame_controls, values=["Vertical (XZ)", "Horizontal (XY)", "Both"], command=lambda v: self.update_plot())
        self.combo_view.set("Vertical (XZ)")
        self.combo_view.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        self.lbl_plot = ctk.CTkLabel(self.frame_controls, text="Plot:", cursor="hand2")
        self.lbl_plot.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_plot, "Coordinate system:\nPolar (Directional)\nCartesian (Rectangular analysis)")
        self.seg_plot_type = ctk.CTkSegmentedButton(self.frame_controls, values=["Polar", "Cartesian"], command=lambda v: self.update_plot())
        self.seg_plot_type.set("Polar")
        self.seg_plot_type.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        self.chk_3d = ctk.CTkCheckBox(self.frame_controls, text="Show 3D Orientation", command=self.update_plot)
        self.chk_3d.select()  # Default to checked
        self.chk_3d.grid(row=r, column=0, columnspan=2, pady=(5, 5))
        r += 1

        # Separator
        ctk.CTkFrame(self.frame_controls, height=2, fg_color="gray30").grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        r += 1

        # Sliders
        self._add_slider(r, "Dynamic Range:", "slider_range", 10, 80, 40, "Minimum dB floor for the plot.\nValues below this are clipped.")
        r += 2
        self._add_slider(r, "Tick Step:", "slider_tick", 1, 20, 10, "Spacing between magnitude grid lines (dB).")
        r += 2
        self._add_slider(r, "Angle Step:", "slider_angle", 5, 90, 30, "Spacing between angle grid lines (degrees).")
        r += 2

        # Calculate Button
        self.btn_calc = ctk.CTkButton(self.frame_controls, text="CALCULATE PATTERN", 
                                      command=self.update_plot,
                                      height=40,
                                      fg_color="#1f6aa5")
        self.btn_calc.grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=20)
        r += 1
        
        self.lbl_status = ctk.CTkLabel(self.frame_controls, text="Ready.", text_color="gray", height=30)
        self.lbl_status.grid(row=r, column=0, columnspan=2, pady=5)

        # PLOT FRAME
        self.frame_plot = ctk.CTkFrame(self)
        self.frame_plot.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame_plot)
        self.toolbar.update()

        self.update_plot()

    # --- HELPER FUNCTIONS ---
    def _add_input(self, r, label, var_name, default, tooltip):
        lbl = ctk.CTkLabel(self.frame_controls, text=label, cursor="hand2")
        lbl.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(lbl, tooltip)
        entry = ctk.CTkEntry(self.frame_controls)
        entry.insert(0, default)
        entry.bind("<FocusOut>", lambda e: self.update_plot())
        entry.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        setattr(self, var_name, entry)

    def _add_slider(self, r, label, var_name, vmin, vmax, vdef, tooltip):
        lbl = ctk.CTkLabel(self.frame_controls, text=label, cursor="hand2")
        lbl.grid(row=r, column=0, sticky="w", padx=10)
        CTkTooltip(lbl, tooltip)
        val_lbl = ctk.CTkLabel(self.frame_controls, text=f"{vdef}{'°' if 'Angle' in label else ' dB'}", width=50, anchor="e")
        val_lbl.grid(row=r, column=1, sticky="e", padx=10)
        
        is_deg = "Angle" in label
        def handler(val):
            val_lbl.configure(text=f"{int(val)}{'°' if is_deg else ' dB'}")
            if self.update_timer: self.after_cancel(self.update_timer)
            self.update_timer = self.after(400, self.update_plot)
        
        slider = ctk.CTkSlider(self.frame_controls, from_=vmin, to=vmax, number_of_steps=(vmax-vmin))
        slider.set(vdef)
        slider.configure(command=handler)
        slider.grid(row=r+1, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        setattr(self, var_name, slider)

    # --- CURSOR LOGIC ---
    def get_cursor_data(self, event):
        """Returns (angle_deg, db_val, plot_x) from event."""
        if self.combo_view.get() == "Both": return None
        if self.seg_plot_type.get() == "Polar":
            if event.xdata is None: return None
            if self.theta_plot is None or self.af_db_plot is None:
                return None
            # Polar xdata is radians. Normalize to 0-360
            angle_rad = (event.xdata % (2 * np.pi) + 2 * np.pi) % (2 * np.pi)
            db = np.interp(angle_rad, self.theta_plot, self.af_db_plot)
            return np.rad2deg(angle_rad), db, angle_rad
        else:
            if event.xdata is None: return None
            if self.theta_deg_sorted is None or self.af_db_sorted is None:
                return None
            angle_deg = event.xdata
            db = np.interp(angle_deg, self.theta_deg_sorted, self.af_db_sorted)
            return angle_deg, db, angle_deg

    def update_cursor_visuals(self, angle_deg, db, x_val, is_fixed=False):
        # Calculate Directivity at this angle: D_max(dBi) + NormalizedPattern(dB)
        d_max_dbi = 10 * np.log10(self.D_linear) if self.D_linear > 0 else 0
        local_dir = d_max_dbi + db
        
        # Text format
        text = f"Cursor: ({angle_deg:.1f}°, {db:.1f}dB), D(°)={local_dir:.2f}dBi"
        self.cursor_text = text
        
        # Update point/line
        if is_fixed:
            if self.cursor_point: self.cursor_point.set_visible(True)
            if self.cursor_line: self.cursor_line.set_visible(True)
            # We create distinct objects for fixed if needed, but for simplicity reusing logic
            # Here we just update title
        else:
            if self.cursor_point:
                self.cursor_point.set_offsets([x_val, db])
                self.cursor_point.set_visible(True)
            else:
                self.cursor_point = self.ax.scatter(x_val, db, color='red', s=50, zorder=10)
                
            if self.cursor_line:
                self.cursor_line.set_xdata([x_val, x_val])
                self.cursor_line.set_visible(True)
            else:
                self.cursor_line = self.ax.axvline(x_val, color='red', linestyle='--')

        self.update_title()

    def on_mouse_move(self, event):
        if event.inaxes != self.ax or self.fixed_cursor or self.theta_plot is None: return
        data = self.get_cursor_data(event)
        if data:
            self.last_x, self.last_db = data[2], data[1]
            self.update_cursor_visuals(data[0], data[1], data[2])
            self.canvas.draw_idle() # Optimized redraw

    def on_click(self, event):
        if event.inaxes != self.ax: return
        
        if self.fixed_cursor:
            # Unfreeze
            self.fixed_cursor = False
            self.cursor_text = ""
            self.update_title()
        else:
            # Freeze
            data = self.get_cursor_data(event)
            if data:
                self.fixed_cursor = True
                self.fixed_x, self.fixed_db = data[2], data[1]
                self.update_cursor_visuals(data[0], data[1], data[2], is_fixed=True)
                self.canvas.draw()

    def update_title(self):
        if hasattr(self, 'ax') and self.ax:
            # Preserve original title lines (first 2)
            curr = self.ax.get_title()
            base = '\n'.join(curr.split('\n')[:2])
            self.ax.set_title(f"{base}\n{self.cursor_text}", va='bottom', fontsize=10)

    # --- 3D INSET LOGIC ---
    def draw_3d_inset(self, is_horizontal, array_axis='Z'):
        """Draws a mini 3D plot to show array orientation and cut plane."""
        # Add inset axes at bottom left (0.0, 0.0) with width/height 0.25
        ax_geo = self.fig.add_axes([0.0, 0.0, 0.25, 0.25], projection='3d')
        ax_geo.set_axis_off() # floating look
        
        # 1. Draw Axis Arrows
        len_arrow = 1.5
        ax_geo.quiver(0,0,0, len_arrow,0,0, color='r', arrow_length_ratio=0.2, linewidth=1) # X
        ax_geo.text(len_arrow, 0, 0, "X", color='r', fontsize=8)
        
        ax_geo.quiver(0,0,0, 0,len_arrow,0, color='g', arrow_length_ratio=0.2, linewidth=1) # Y
        ax_geo.text(0, len_arrow, 0, "Y", color='g', fontsize=8)
        
        # Array Axis - Thicker
        if array_axis == 'X':
            ax_geo.quiver(0,0,0, len_arrow,0,0, color='k', arrow_length_ratio=0.2, linewidth=2) # X
            ax_geo.text(len_arrow, 0, 0, "X", color='k', fontsize=9, fontweight='bold')
            # Draw Array Elements along X
            x_pos = np.linspace(-0.8, 0.8, 4)
            ax_geo.scatter(x_pos, np.zeros(4), np.zeros(4), color='black', s=15, alpha=1.0)
        else:
            ax_geo.quiver(0,0,0, 0,0,len_arrow, color='k', arrow_length_ratio=0.2, linewidth=2) # Z
            ax_geo.text(0, 0, len_arrow, "Z", color='k', fontsize=9, fontweight='bold')
            # Draw Array Elements along Z
            z_pos = np.linspace(-0.8, 0.8, 4)
            ax_geo.scatter(np.zeros(4), np.zeros(4), z_pos, color='black', s=15, alpha=1.0)

        # 3. Draw Cut Plane Indicator
        t = np.linspace(0, 2*np.pi, 60)
        if is_horizontal:
            # Horizontal (XY) cut - Blue circle on XY plane
            x_c = np.cos(t)
            y_c = np.sin(t)
            z_c = np.zeros_like(t)
            # Filled disk
            r = np.linspace(0, 1, 2)
            R, T = np.meshgrid(r, t)
            X_fill = R * np.cos(T)
            Y_fill = R * np.sin(T)
            Z_fill = np.zeros_like(X_fill)
            ax_geo.plot_surface(X_fill, Y_fill, Z_fill, color='blue', alpha=0.2)
            # Outline
            ax_geo.plot(x_c, y_c, z_c, color='blue', linestyle='--', linewidth=1, alpha=0.8)
        else:
            # Vertical (XZ) cut - Blue circle on XZ plane
            x_c = np.sin(t)
            y_c = np.zeros_like(t)
            z_c = np.cos(t)
            # Filled disk
            r = np.linspace(0, 1, 2)
            R, T = np.meshgrid(r, t)
            X_fill = R * np.sin(T)
            Y_fill = np.zeros_like(X_fill)
            Z_fill = R * np.cos(T)
            ax_geo.plot_surface(X_fill, Y_fill, Z_fill, color='blue', alpha=0.2)
            # Outline
            ax_geo.plot(x_c, y_c, z_c, color='blue', linestyle='--', linewidth=1, alpha=0.8)

        # Set limits to ensure aspect ratio
        limit = 1.2
        ax_geo.set_xlim(-limit, limit)
        ax_geo.set_ylim(-limit, limit)
        ax_geo.set_zlim(-limit, limit)
        ax_geo.view_init(elev=20, azim=45)

    def update_plot(self):
        if hasattr(self, 'update_timer') and self.update_timer:
            self.after_cancel(self.update_timer)
            self.update_timer = None
        try:
            # Disable cursor data during update to prevent race conditions
            self.theta_plot = None
            self.af_db_plot = None
            self.theta_deg_sorted = None
            self.af_db_sorted = None
            
            # Inputs
            N = int(self.entry_n.get())
            d = float(self.entry_d.get())
            beta = float(self.entry_beta.get())
            el_type = self.combo_type.get()
            view = self.combo_view.get()
            dyn_range = int(self.slider_range.get())
            tick_step = int(self.slider_tick.get())
            angle_step = int(self.slider_angle.get())
            plot_type = self.seg_plot_type.get()
            array_axis = self.combo_axis.get()

            # Currents Parsing
            raw_curr = self.entry_currents.get()
            if not raw_curr.strip(): curr_list = [1.0]
            else: curr_list = [float(x.strip()) for x in raw_curr.split(',')]
            
            if len(curr_list) == 1: currents = np.full(N, curr_list[0])
            elif len(curr_list) == N: currents = np.array(curr_list)
            else: raise ValueError(f"Currents count != N")

            # --- CALCULATE METRICS (PHYSICS) ---
            # Calculate D and HPBW using the rigorous 3D model (Vertical cut integration)
            d_lin, hpbw = self.calculator.calculate_metrics(N, d, beta, el_type, currents, array_axis=array_axis)
            d_dbi = 10 * np.log10(d_lin) if d_lin > 0 else 0
            
            self.D_linear = d_lin # Store for cursor usage

            # --- CALCULATE PLOT DATA (VISUALIZATION) ---
            theta, total_linear = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view=view, array_axis=array_axis)
            af_db = self.calculator.convert_to_db(total_linear, dynamic_range=dyn_range)
            
            # Store for cursor interpolation
            self.theta_plot = theta
            self.af_db_plot = af_db

            # --- PLOTTING ---
            self.fig.clear()
            self.cursor_point = None
            self.cursor_line = None
            self.cursor_text = ""

            if plot_type == "Polar":
                if view == "Both":
                    # Plot both views side by side
                    self.ax1 = self.fig.add_subplot(121, projection='polar')
                    self.ax2 = self.fig.add_subplot(122, projection='polar')
                    
                    # Vertical view
                    theta_v, total_linear_v = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view="Vertical (XZ)", array_axis=array_axis)
                    af_db_v = self.calculator.convert_to_db(total_linear_v, dynamic_range=dyn_range)
                    
                    self.ax1.set_theta_zero_location('N')
                    self.ax1.set_theta_direction(-1)
                    angles_deg = np.arange(0, 360, angle_step)
                    labels = [f'{deg}°' if deg <= 180 else f'{deg - 360}°' for deg in angles_deg]
                    if 180 in angles_deg: labels[list(angles_deg).index(180)] = '±180°'
                    self.ax1.set_thetagrids(angles_deg, labels)
                    self.ax1.set_ylim(-dyn_range, 0)
                    self.ax1.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax1.plot(theta_v, af_db_v, color='#1f77b4', linewidth=2)
                    self.ax1.grid(True, alpha=0.5)
                    
                    directions_v = {'+X': 0, '+Z': np.pi/2, '-X': np.pi, '-Z': 3*np.pi/2}
                    for label, theta_rad in directions_v.items():
                        self.ax1.text(theta_rad, -2, label, ha='center', va='center', fontsize=8, color='red', fontweight='bold')
                    self.ax1.set_title("Vertical (XZ)", fontsize=10)
                    
                    # Horizontal view
                    theta_h, total_linear_h = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view="Horizontal (XY)", array_axis=array_axis)
                    af_db_h = self.calculator.convert_to_db(total_linear_h, dynamic_range=dyn_range)
                    
                    self.ax2.set_theta_zero_location('N')
                    self.ax2.set_theta_direction(-1)
                    self.ax2.set_thetagrids(angles_deg, labels)
                    self.ax2.set_ylim(-dyn_range, 0)
                    self.ax2.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax2.plot(theta_h, af_db_h, color='#1f77b4', linewidth=2)
                    self.ax2.grid(True, alpha=0.5)
                    
                    directions_h = {'+X': 0, '+Y': np.pi/2, '-X': np.pi, '-Y': 3*np.pi/2}
                    for label, theta_rad in directions_h.items():
                        self.ax2.text(theta_rad, -2, label, ha='center', va='center', fontsize=8, color='red', fontweight='bold')
                    self.ax2.set_title("Horizontal (XY)", fontsize=10)
                    
                    self.ax = self.ax1  # For cursor, use first one, but disable cursor for Both
                    
                else:
                    self.ax = self.fig.add_subplot(111, projection='polar')
                    self.ax.set_theta_zero_location('N')
                    self.ax.set_theta_direction(-1)
                    
                    angles_deg = np.arange(0, 360, angle_step)
                    labels = [f'{deg}°' if deg <= 180 else f'{deg - 360}°' for deg in angles_deg]
                    if 180 in angles_deg: labels[list(angles_deg).index(180)] = '±180°'
                    self.ax.set_thetagrids(angles_deg, labels)
                    
                    self.ax.set_ylim(-dyn_range, 0)
                    self.ax.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax.plot(theta, af_db, color='#1f77b4', linewidth=2)
                    self.ax.grid(True, alpha=0.5)
                    
                    # Add axis direction labels
                    if "Horizontal" in view:
                        directions = {'+X': 0, '+Y': np.pi/2, '-X': np.pi, '-Y': 3*np.pi/2}
                    else:
                        directions = {'+X': 0, '+Z': np.pi/2, '-X': np.pi, '-Z': 3*np.pi/2}
                    for label, theta_rad in directions.items():
                        self.ax.text(theta_rad, -2, label, ha='center', va='center', fontsize=8, color='red', fontweight='bold')
            
            else: # Cartesian
                if view == "Both":
                    # Plot both views side by side
                    self.ax1 = self.fig.add_subplot(121)
                    
                    # Vertical view
                    theta_v, total_linear_v = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view="Vertical (XZ)", array_axis=array_axis)
                    af_db_v = self.calculator.convert_to_db(total_linear_v, dynamic_range=dyn_range)
                    theta_deg_v = np.rad2deg(theta_v)
                    theta_deg_shifted_v = np.copy(theta_deg_v)
                    theta_deg_shifted_v[theta_deg_shifted_v >= 180] -= 360
                    sort_idx_v = np.argsort(theta_deg_shifted_v)
                    self.theta_deg_sorted_v = theta_deg_shifted_v[sort_idx_v]
                    self.af_db_sorted_v = af_db_v[sort_idx_v]
                    
                    self.ax1.plot(self.theta_deg_sorted_v, self.af_db_sorted_v, color='#1f77b4', linewidth=2)
                    self.ax1.set_xlabel("Angle (°)")
                    self.ax1.set_ylabel("Normalized Power (dB)")
                    self.ax1.set_xlim(-180, 180)
                    self.ax1.margins(x=0)
                    self.ax1.set_ylim(-dyn_range, 0)
                    self.ax1.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax1.set_xticks(np.arange(-180, 181, angle_step))
                    self.ax1.grid(True, alpha=0.5)
                    self.ax1.set_title("Vertical (XZ)", fontsize=10)
                    
                    self.ax2 = self.fig.add_subplot(122)
                    
                    # Horizontal view
                    theta_h, total_linear_h = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view="Horizontal (XY)", array_axis=array_axis)
                    af_db_h = self.calculator.convert_to_db(total_linear_h, dynamic_range=dyn_range)
                    theta_deg_h = np.rad2deg(theta_h)
                    theta_deg_shifted_h = np.copy(theta_deg_h)
                    theta_deg_shifted_h[theta_deg_shifted_h >= 180] -= 360
                    sort_idx_h = np.argsort(theta_deg_shifted_h)
                    self.theta_deg_sorted_h = theta_deg_shifted_h[sort_idx_h]
                    self.af_db_sorted_h = af_db_h[sort_idx_h]
                    
                    self.ax2.plot(self.theta_deg_sorted_h, self.af_db_sorted_h, color='#1f77b4', linewidth=2)
                    self.ax2.set_xlabel("Angle (°)")
                    self.ax2.set_ylabel("Normalized Power (dB)")
                    self.ax2.set_xlim(-180, 180)
                    self.ax2.margins(x=0)
                    self.ax2.set_ylim(-dyn_range, 0)
                    self.ax2.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax2.set_xticks(np.arange(-180, 181, angle_step))
                    self.ax2.grid(True, alpha=0.5)
                    self.ax2.set_title("Horizontal (XY)", fontsize=10)
                    
                    self.ax = self.ax1  # For compatibility
                    
                else:
                    self.ax = self.fig.add_subplot(111)
                    theta_deg = np.rad2deg(theta)
                    # Sort for interpolation/plotting
                    theta_deg_shifted = np.copy(theta_deg)
                    theta_deg_shifted[theta_deg_shifted >= 180] -= 360
                    sort_idx = np.argsort(theta_deg_shifted)
                    
                    self.theta_deg_sorted = theta_deg_shifted[sort_idx]
                    self.af_db_sorted = af_db[sort_idx]
                    
                    self.ax.plot(self.theta_deg_sorted, self.af_db_sorted, color='#1f77b4', linewidth=2)
                    self.ax.set_xlabel("Angle (°)")
                    self.ax.set_ylabel("Normalized Power (dB)")
                    self.ax.set_xlim(-180, 180)
                    self.ax.margins(x=0)
                    self.ax.set_ylim(-dyn_range, 0)
                    self.ax.set_yticks(np.arange(-dyn_range, 1, tick_step))
                    self.ax.set_xticks(np.arange(-180, 181, angle_step))
                    self.ax.grid(True, alpha=0.5)

            if view == "Both":
                title_text = f"Pattern (Both Views, Array on {array_axis}): {el_type}\nN={N}, d={d}λ, β={beta}°, Dmax={d_dbi:.2f}dBi, HPBW={hpbw:.1f}°"
                self.fig.suptitle(title_text, fontsize=10)
            else:
                title_text = f"Pattern ({view}, Array on {array_axis}): {el_type}\nN={N}, d={d}λ, β={beta}°, Dmax={d_dbi:.2f}dBi, HPBW={hpbw:.1f}°"
                self.ax.set_title(title_text, va='bottom', fontsize=10)
            self.fig.tight_layout()
            
            # --- DRAW 3D ORIENTATION INSET ---
            if self.chk_3d.get():
                is_horiz = "Horizontal" in view if view != "Both" else False  # Default to Vertical for Both
                axis_letter = 'X' if 'X' in array_axis else 'Z'
                self.draw_3d_inset(is_horiz, axis_letter)
            
            # Update fixed cursor db and create visuals
            if self.fixed_cursor and self.fixed_x is not None:
                if plot_type == "Polar":
                    db = np.interp(self.fixed_x, self.theta_plot, self.af_db_plot)
                else:
                    db = np.interp(self.fixed_x, self.theta_deg_sorted, self.af_db_sorted)
                self.fixed_db = db
                self.cursor_point = self.ax.scatter(self.fixed_x, self.fixed_db, color='red', s=50, zorder=10)
                self.cursor_line = self.ax.axvline(self.fixed_x, color='red', linestyle='--')
                # Update cursor text
                angle_deg = np.rad2deg(self.fixed_x) if plot_type == "Polar" else self.fixed_x
                dir_dbi = 10 * np.log10(self.D_linear) + db if self.D_linear > 0 else 0
                self.cursor_text = f"Cursor: ({angle_deg:.1f}°, {db:.1f}dB), D(°)={dir_dbi:.1f}dBi"
            
            self.update_title()
            self.canvas.draw()
            
            self.lbl_status.configure(text="Calculation successful.", text_color="green", font=("Arial", 12, "normal"))

        except ValueError as ve:
            self.lbl_status.configure(text=f"Error: {str(ve)}", text_color="red", font=("Arial", 13, "bold"))
        except Exception as e:
            self.lbl_status.configure(text=f"Error: {str(e)}", text_color="red", font=("Arial", 13, "bold"))