import customtkinter as ctk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import tkinter as tk
from src.back import AntennaCalculator

# --- CUSTOM TOOLTIP CLASS ---
class CTkTooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.widget.bind("<Enter>", self.show_tooltip)
        self.widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event=None):
        if self.tooltip_window or not self.text:
            return
        
        # Coordinates
        x = self.widget.winfo_rootx() + 25
        y = self.widget.winfo_rooty() + 25
        
        
        # Create floating window
        self.tooltip_window = tk.Toplevel(self.widget)
        self.tooltip_window.wm_overrideredirect(True) # No borders
        self.tooltip_window.wm_geometry(f"+{x}+{y}")
        
        # Internal label
        label = tk.Label(self.tooltip_window, text=self.text, justify='left',
                         background="#2b2b2b", fg="#ffffff", # Dark Theme Colors
                         relief='solid', borderwidth=1,
                         font=("Arial", 10, "normal"),
                         padx=5, pady=2)
        label.pack(ipadx=1)

    def hide_tooltip(self, event=None):
        if self.tooltip_window:
            self.tooltip_window.destroy()
            self.tooltip_window = None

# --- MAIN APPLICATION ---
class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Main window configuration
        self.title("Antenna Array Visualizer - ITBA 22.21")
        self.geometry("1100x750")
        ctk.set_appearance_mode("Dark")
        ctk.set_default_color_theme("blue")

        self.calculator = AntennaCalculator()
        self.update_timer = None

        # --- MAIN LAYOUT ---
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # CONTROLS FRAME
        self.frame_controls = ctk.CTkFrame(self, width=280)
        self.frame_controls.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        self.frame_controls.grid_columnconfigure(0, weight=1) 
        self.frame_controls.grid_columnconfigure(1, weight=2) 

        # --- PARAMETERS (GRID LAYOUT) ---
        r = 0 

        # Title
        self.lbl_title = ctk.CTkLabel(self.frame_controls, text="Parameters", font=("Arial", 20, "bold"))
        self.lbl_title.grid(row=r, column=0, columnspan=2, pady=(15, 15))
        r += 1
        
        # Hint
        self.lbl_hint = ctk.CTkLabel(self.frame_controls, text="(Hover labels for help)", font=("Arial", 12), text_color="gray")
        self.lbl_hint.grid(row=r, column=0, columnspan=2, pady=(0,5))
        r += 1

        # 1. N Antennas
        self.lbl_n = ctk.CTkLabel(self.frame_controls, text="N Antennas:", cursor="hand2") 
        self.lbl_n.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_n, "Total number of radiating elements\nin the array.")
        
        self.entry_n = ctk.CTkEntry(self.frame_controls)
        self.entry_n.insert(0, "4")
        self.entry_n.bind("<FocusOut>", lambda e: self.update_plot())
        self.entry_n.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # 2. Separation
        self.lbl_d = ctk.CTkLabel(self.frame_controls, text="Separation (d/λ):", cursor="hand2")
        self.lbl_d.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_d, "Distance between elements relative\nto wavelength (e.g., 0.5 = λ/2).")

        self.entry_d = ctk.CTkEntry(self.frame_controls)
        self.entry_d.insert(0, "0.5")
        self.entry_d.bind("<FocusOut>", lambda e: self.update_plot())
        self.entry_d.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # 3. Phase Diff
        self.lbl_beta = ctk.CTkLabel(self.frame_controls, text="Phase (β°):", cursor="hand2")
        self.lbl_beta.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_beta, "Progressive phase shift between\nconsecutive elements (degrees).\nUsed for beam steering.")

        self.entry_beta = ctk.CTkEntry(self.frame_controls)
        self.entry_beta.insert(0, "0")
        self.entry_beta.bind("<FocusOut>", lambda e: self.update_plot())
        self.entry_beta.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1
        
        # 4. Currents
        self.lbl_currents = ctk.CTkLabel(self.frame_controls, text="Intensities:", cursor="hand2")
        self.lbl_currents.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_currents, "Current intensity for each element.\nExample for N=4: '1, 0.5, 0.5, 0.33'\nSingle value applies to all.")
        
        self.entry_currents = ctk.CTkEntry(self.frame_controls)
        self.entry_currents.insert(0, "1") 
        self.entry_currents.bind("<FocusOut>", lambda e: self.update_plot())
        self.entry_currents.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # 5. Element Type
        self.lbl_type = ctk.CTkLabel(self.frame_controls, text="Type:", cursor="hand2")
        self.lbl_type.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_type, "Radiation pattern of the individual element.\n(Multiplied by Array Factor).")

        self.combo_type = ctk.CTkComboBox(self.frame_controls, 
                                          values=["Isotropic", "Dipole (λ/2)", "Monopole (λ/4)"],
                                          command=lambda value: self.update_plot())
        self.combo_type.set("Isotropic")
        self.combo_type.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # View (Plane)
        self.lbl_view = ctk.CTkLabel(self.frame_controls, text="View:", cursor="hand2")
        self.lbl_view.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_view, "Cut Plane:\nVertical (Elevation/Theta)\nHorizontal (Azimuth/Phi)")

        self.combo_view = ctk.CTkComboBox(self.frame_controls,
                                          values=["Vertical (XZ)", "Horizontal (XY)"],
                                          command=lambda value: self.update_plot())
        self.combo_view.set("Vertical (XZ)")
        self.combo_view.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # Plot Type
        self.lbl_plot_type = ctk.CTkLabel(self.frame_controls, text="Plot:", cursor="hand2")
        self.lbl_plot_type.grid(row=r, column=0, sticky="e", padx=5, pady=5)
        CTkTooltip(self.lbl_plot_type, "Coordinate system:\nPolar (Directional)\nCartesian (Rectangular analysis)")

        self.seg_plot_type = ctk.CTkSegmentedButton(self.frame_controls,
                                                     values=["Polar", "Cartesian"],
                                                     command=lambda value: self.update_plot())
        self.seg_plot_type.set("Polar")
        self.seg_plot_type.grid(row=r, column=1, sticky="ew", padx=5, pady=5)
        r += 1

        # --- SLIDERS ---
        ctk.CTkFrame(self.frame_controls, height=2, fg_color="gray30").grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        r += 1

        # Dynamic Range
        self.lbl_range = ctk.CTkLabel(self.frame_controls, text="Dynamic Range:", cursor="hand2")
        self.lbl_range.grid(row=r, column=0, sticky="w", padx=10)
        CTkTooltip(self.lbl_range, "Minimum dB floor for the plot.\nValues below this are clipped.")
        
        self.lbl_range_val = ctk.CTkLabel(self.frame_controls, text="40 dB", width=50, anchor="e")
        self.lbl_range_val.grid(row=r, column=1, sticky="e", padx=10)
        r += 1
        
        self.slider_range = ctk.CTkSlider(self.frame_controls, from_=10, to=80, number_of_steps=70)
        self.slider_range.set(40)
        self.slider_range.configure(command=self._slider_event_handler)
        self.slider_range.grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        r += 1

        # Tick Step
        self.lbl_tick = ctk.CTkLabel(self.frame_controls, text="Tick Step:", cursor="hand2")
        self.lbl_tick.grid(row=r, column=0, sticky="w", padx=10)
        CTkTooltip(self.lbl_tick, "Spacing between magnitude grid lines (dB).")

        self.lbl_tick_val = ctk.CTkLabel(self.frame_controls, text="10 dB", width=50, anchor="e")
        self.lbl_tick_val.grid(row=r, column=1, sticky="e", padx=10)
        r += 1

        self.slider_tick = ctk.CTkSlider(self.frame_controls, from_=1, to=20, number_of_steps=19)
        self.slider_tick.set(10)
        self.slider_tick.configure(command=self._tick_slider_event_handler)
        self.slider_tick.grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        r += 1

        # Angle Step
        self.lbl_angle = ctk.CTkLabel(self.frame_controls, text="Angle Step:", cursor="hand2")
        self.lbl_angle.grid(row=r, column=0, sticky="w", padx=10)
        CTkTooltip(self.lbl_angle, "Spacing between angle grid lines (degrees).")

        self.lbl_angle_val = ctk.CTkLabel(self.frame_controls, text="30°", width=50, anchor="e")
        self.lbl_angle_val.grid(row=r, column=1, sticky="e", padx=10)
        r += 1

        self.slider_angle = ctk.CTkSlider(self.frame_controls, from_=5, to=90, number_of_steps=17)
        self.slider_angle.set(30)
        self.slider_angle.configure(command=self._angle_slider_event_handler)
        self.slider_angle.grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=(0, 10))
        r += 1

        # Calculate Button
        self.btn_calc = ctk.CTkButton(self.frame_controls, text="CALCULATE PATTERN", 
                                      command=self.update_plot,
                                      height=40,
                                      fg_color="#1f6aa5")
        self.btn_calc.grid(row=r, column=0, columnspan=2, sticky="ew", padx=10, pady=20)
        r += 1
        
        # Status Label
        self.lbl_status = ctk.CTkLabel(self.frame_controls, text="Ready.", text_color="gray")
        self.lbl_status.grid(row=r, column=0, columnspan=2, pady=5)


        # PLOT FRAME
        self.frame_plot = ctk.CTkFrame(self)
        self.frame_plot.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        # --- GRAPH CONFIGURATION ---
        self.fig = Figure(figsize=(5, 5), dpi=100)
        
        # 1. Canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # 2. Toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame_plot)
        self.toolbar.update()

        self.cursor_line = None
        self.fixed_cursor = False
        self.cursor_text = ""
        self.D = 0
        self.fixed_x = None
        self.fixed_db = None
        self.cursor_point = None
        self.theta_plot = None
        self.af_db_plot = None
        self.theta_deg_sorted = None
        self.af_db_sorted = None
        self.last_x = None
        self.last_db = None

        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        self.update_plot()

    # --- HANDLERS ---
    def _slider_event_handler(self, value):
        self.lbl_range_val.configure(text=f"{int(value)} dB")
        if self.update_timer: self.after_cancel(self.update_timer)
        self.update_timer = self.after(400, self.update_plot)

    def _tick_slider_event_handler(self, value):
        self.lbl_tick_val.configure(text=f"{int(value)} dB")
        if self.update_timer: self.after_cancel(self.update_timer)
        self.update_timer = self.after(400, self.update_plot)

    def _angle_slider_event_handler(self, value):
        self.lbl_angle_val.configure(text=f"{int(value)}°")
        if self.update_timer: self.after_cancel(self.update_timer)
        self.update_timer = self.after(400, self.update_plot)

    def on_mouse_move(self, event):
        if event.inaxes == self.ax and not self.fixed_cursor and self.theta_plot is not None:
            if self.seg_plot_type.get() == "Polar":
                if event.xdata is not None:
                    angle_rad = (event.xdata % (2 * np.pi) + 2 * np.pi) % (2 * np.pi)
                    db = np.interp(angle_rad, self.theta_plot, self.af_db_plot)
                    angle_deg = np.rad2deg(angle_rad)
                    x_val = angle_rad
                else:
                    return
            else:
                if event.xdata is not None:
                    angle_deg = event.xdata
                    db = np.interp(angle_deg, self.theta_deg_sorted, self.af_db_sorted)
                    x_val = angle_deg
                else:
                    return
            dir_dbi = 10 * np.log10(self.D) + db if self.D > 0 else 0
            self.cursor_text = f"Cursor: ({angle_deg:.1f}°, {db:.2}dBi), D(°)={dir_dbi:.2f}dBi"
            self.last_x = x_val
            self.last_db = db
            if self.cursor_point:
                self.cursor_point.set_offsets([x_val, db])
            else:
                self.cursor_point = self.ax.scatter(x_val, db, color='red', s=50, zorder=10)
            if self.cursor_line:
                self.cursor_line.set_xdata([x_val, x_val])
            else:
                self.cursor_line = self.ax.axvline(x_val, color='red', linestyle='--')
            self.canvas.draw()
            self.update_title()

    def on_click(self, event):
        if event.inaxes == self.ax:
            if self.seg_plot_type.get() == "Polar":
                if event.xdata is not None:
                    angle_rad = (event.xdata % (2 * np.pi) + 2 * np.pi) % (2 * np.pi)
                    db = np.interp(angle_rad, self.theta_plot, self.af_db_plot)
                    angle_deg = np.rad2deg(angle_rad)
                else:
                    return
            else:
                if event.xdata is not None:
                    angle_deg = event.xdata
                    db = np.interp(angle_deg, self.theta_deg_sorted, self.af_db_sorted)
                else:
                    return
            if self.fixed_cursor:
                # Save current fixed position as last mobile
                if self.fixed_x is not None and self.fixed_db is not None:
                    self.last_x = self.fixed_x
                    self.last_db = self.fixed_db
                self.fixed_cursor = False
                self.cursor_text = ""
                self.fixed_x = None
                self.fixed_db = None
                if self.cursor_line:
                    self.cursor_line.set_visible(False)
                if self.cursor_point:
                    self.cursor_point.set_visible(False)
                # Restore last mobile cursor
                if self.last_x is not None and self.last_db is not None:
                    self.cursor_point = self.ax.scatter(self.last_x, self.last_db, color='red', s=50, zorder=10)
                    self.cursor_line = self.ax.axvline(self.last_x, color='red', linestyle='--')
                    angle_deg_for_text = np.rad2deg(self.last_x) if self.seg_plot_type.get() == "Polar" else self.last_x
                    dir_dbi = 10 * np.log10(self.D) + self.last_db if self.D > 0 else 0
                    self.cursor_text = f"Cursor: {angle_deg_for_text:.1f}°, {self.last_db:.1f} dB, Dir: {dir_dbi:.1f} dBi"
            else:
                self.fixed_cursor = True
                if self.seg_plot_type.get() == "Polar":
                    x_val = angle_rad
                else:
                    x_val = angle_deg
                dir_dbi = 10 * np.log10(self.D) + db if self.D > 0 else 0
                self.cursor_text = f"Cursor: {angle_deg:.1f}°, {db:.1f} dB, Dir: {dir_dbi:.1f} dBi"
                self.fixed_x = x_val
                self.fixed_db = db
                if self.cursor_point:
                    self.cursor_point.set_offsets([x_val, db])
                else:
                    self.cursor_point = self.ax.scatter(x_val, db, color='red', s=50, zorder=10)
                if self.cursor_line:
                    self.cursor_line.set_xdata([x_val, x_val])
                else:
                    self.cursor_line = self.ax.axvline(x_val, color='red', linestyle='--')
            self.canvas.draw()
            self.update_title()

    def update_title(self):
        if hasattr(self, 'ax') and self.ax:
            current_title = self.ax.get_title()
            lines = current_title.split('\n')
            base_title = '\n'.join(lines[:2])  # First two lines
            self.ax.set_title(f"{base_title}\n{self.cursor_text}", va='bottom', fontsize=10)
            self.canvas.draw()

    def update_plot(self):
        if hasattr(self, 'update_timer') and self.update_timer:
            self.after_cancel(self.update_timer)
            self.update_timer = None
        try:
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

            # Currents
            raw_currents = self.entry_currents.get()
            if not raw_currents.strip(): 
                current_list = [1.0]
            else:
                current_list = [float(x) for x in raw_currents.split(',')]
            
            if len(current_list) == 1:
                currents = np.full(N, current_list[0])
            elif len(current_list) == N:
                currents = np.array(current_list)
            else:
                raise ValueError(f"Currents count ({len(current_list)}) must match N ({N})")

            # Calculate
            theta, total_linear = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view=view)

            # Convert to dB
            af_db = self.calculator.convert_to_db(total_linear, dynamic_range=dyn_range)

            D = self.calculator.calculate_directivity(total_linear, theta)
            hpbw = self.calculator.calculate_hpbw(total_linear, theta)

            self.D = D
            self.theta_plot = theta
            self.af_db_plot = af_db

            # Plotting
            self.fig.clear()

            if plot_type == "Polar":
                self.ax = self.fig.add_subplot(111, projection='polar')
                self.ax.set_theta_zero_location('N')
                self.ax.set_theta_direction(-1)
                
                angles_deg = np.arange(0, 360, angle_step)
                labels = [f'{d}°' if d <= 180 else f'{d - 360}°' for d in angles_deg]
                if 180 in angles_deg: labels[list(angles_deg).index(180)] = '±180°'
                self.ax.set_thetagrids(angles_deg, labels)
                
                self.ax.set_ylim(-dyn_range, 0)
                self.ax.set_yticks(np.arange(-dyn_range, 1, tick_step))
                self.ax.plot(theta, af_db, color='#1f77b4', linewidth=2)
            
            else: # Cartesian
                self.ax = self.fig.add_subplot(111)
                theta_deg = np.rad2deg(theta)
                theta_deg_shifted = np.copy(theta_deg)
                theta_deg_shifted[theta_deg_shifted >= 180] -= 360
                
                sort_indices = np.argsort(theta_deg_shifted)
                theta_deg_sorted = theta_deg_shifted[sort_indices]
                af_db_sorted = af_db[sort_indices]
                
                self.theta_deg_sorted = theta_deg_sorted
                self.af_db_sorted = af_db_sorted
                
                self.ax.plot(theta_deg_sorted, af_db_sorted, color='#1f77b4', linewidth=2)
                self.ax.set_xlabel("Angle (°)")
                self.ax.set_ylabel("Normalized Power (dB)")
                self.ax.set_xlim(-180, 180)
                self.ax.margins(x=0)
                self.ax.set_ylim(-dyn_range, 0)
                self.ax.set_yticks(np.arange(-dyn_range, 1, tick_step))
                self.ax.set_xticks(np.arange(-180, 181, angle_step))
                self.ax.tick_params(axis='x', rotation=45)

            # Title & Layout
            title_text = f"Pattern ({'Horizontal' if 'Horizontal' in view else 'Vertical'}): {el_type}\nN={N}, d={d}λ, β={beta}°, Dmax={D:.2f}dBi, HPBW={hpbw:.1f}°\n{self.cursor_text}"
            self.ax.set_title(title_text, va='bottom', fontsize=10)
            self.ax.grid(True, alpha=0.5)
            self.fig.tight_layout()
            
            if self.fixed_cursor and self.fixed_x is not None and self.fixed_db is not None:
                self.cursor_point = self.ax.scatter(self.fixed_x, self.fixed_db, color='red', s=50, zorder=10)
                self.cursor_line = self.ax.axvline(self.fixed_x, color='red', linestyle='--')
            elif not self.fixed_cursor and self.last_x is not None and self.last_db is not None:
                self.cursor_point = self.ax.scatter(self.last_x, self.last_db, color='red', s=50, zorder=10)
                self.cursor_line = self.ax.axvline(self.last_x, color='red', linestyle='--')
            
            self.canvas.draw()
            self.lbl_status.configure(text="Calculation successful.", text_color="green", font=("Arial", 12, "normal"))

            self.cursor_text = ""
            if self.cursor_line and not self.fixed_cursor and self.last_x is None:
                self.cursor_line.set_visible(False)
                self.cursor_line = None
            if self.cursor_point and not self.fixed_cursor and self.last_x is None:
                self.cursor_point.set_visible(False)
                self.cursor_point = None

        except ValueError as ve:
            self.lbl_status.configure(text=f"Error: {str(ve)}", text_color="red", font=("Arial", 13, "bold"))
        except Exception as e:
            self.lbl_status.configure(text=f"Error: {str(e)}", text_color="red", font=("Arial", 13, "bold"))