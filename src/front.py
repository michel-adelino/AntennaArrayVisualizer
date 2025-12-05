import customtkinter as ctk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from src.back import AntennaCalculator

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

        self.lbl_title = ctk.CTkLabel(self.frame_controls, text="Parameters", font=("Arial", 20, "bold"))
        self.lbl_title.pack(pady=20)

        # 1. Number of antennas (N)
        self.lbl_n = ctk.CTkLabel(self.frame_controls, text="Number of Antennas (N):")
        self.lbl_n.pack(pady=(5, 0))
        self.entry_n = ctk.CTkEntry(self.frame_controls)
        self.entry_n.insert(0, "4")
        self.entry_n.pack(pady=5)
        self.entry_n.bind("<FocusOut>", lambda e: self.update_plot())

        # 2. Separation (d/lambda)
        self.lbl_d = ctk.CTkLabel(self.frame_controls, text="Separation (d/λ):")
        self.lbl_d.pack(pady=(5, 0))
        self.entry_d = ctk.CTkEntry(self.frame_controls)
        self.entry_d.insert(0, "0.5")
        self.entry_d.pack(pady=5)
        self.entry_d.bind("<FocusOut>", lambda e: self.update_plot())

        # 3. Phase Difference (Beta)
        self.lbl_beta = ctk.CTkLabel(self.frame_controls, text="Phase Diff (β in degrees):")
        self.lbl_beta.pack(pady=(5, 0))
        self.entry_beta = ctk.CTkEntry(self.frame_controls)
        self.entry_beta.insert(0, "0")
        self.entry_beta.pack(pady=5)
        self.entry_beta.bind("<FocusOut>", lambda e: self.update_plot())
        
        # 4. Currents
        self.lbl_currents = ctk.CTkLabel(self.frame_controls, text="Currents (CSV, e.g., '1, 0.5, 1'):")
        self.lbl_currents.pack(pady=(5, 0))
        self.entry_currents = ctk.CTkEntry(self.frame_controls)
        self.entry_currents.insert(0, "1") # Default uniform distribution
        self.entry_currents.pack(pady=5)
        self.entry_currents.bind("<FocusOut>", lambda e: self.update_plot())

        # 5. Element Type
        self.lbl_type = ctk.CTkLabel(self.frame_controls, text="Element Type:")
        self.lbl_type.pack(pady=(15, 0))
        self.combo_type = ctk.CTkComboBox(self.frame_controls, 
                                          values=["Isotropic", 
                                                  "Dipole (λ/2)", 
                                                  "Monopole (λ/4)"],
                                          command=lambda value: self.update_plot())
        self.combo_type.set("Isotropic")
        self.combo_type.pack(pady=5)

        # View (Plane)
        self.lbl_view = ctk.CTkLabel(self.frame_controls, text="View (Plane):")
        self.lbl_view.pack(pady=(10, 0))
        self.combo_view = ctk.CTkComboBox(self.frame_controls,
                                          values=["Vertical (XZ)", "Horizontal (XY)"],
                                          command=lambda value: self.update_plot())
        self.combo_view.set("Vertical (XZ)")
        self.combo_view.pack(pady=5)

        # Plot Type
        self.lbl_plot_type = ctk.CTkLabel(self.frame_controls, text="Plot Type:")
        self.lbl_plot_type.pack(pady=(10, 0))
        self.seg_plot_type = ctk.CTkSegmentedButton(self.frame_controls,
                                                     values=["Polar", "Cartesian"],
                                                     command=lambda value: self.update_plot())
        self.seg_plot_type.set("Polar")
        self.seg_plot_type.pack(pady=5)

        # Dynamic Range
        self.lbl_range = ctk.CTkLabel(self.frame_controls, text="Dynamic Range (dB):")
        self.lbl_range.pack(pady=(20, 0))
        self.slider_range = ctk.CTkSlider(self.frame_controls, from_=10, to=80, number_of_steps=70)
        self.slider_range.set(40)
        self.slider_range.pack(pady=5)
        
        self.lbl_range_val = ctk.CTkLabel(self.frame_controls, text="40 dB")
        self.lbl_range_val.pack(pady=(0, 5))
        self.slider_range.configure(command=self._slider_event_handler)

        # Tick Step
        self.lbl_tick = ctk.CTkLabel(self.frame_controls, text="Tick Step (dB):")
        self.lbl_tick.pack(pady=(10, 0))
        self.slider_tick = ctk.CTkSlider(self.frame_controls, from_=1, to=20, number_of_steps=19)
        self.slider_tick.set(10)
        self.slider_tick.pack(pady=5)
        
        self.lbl_tick_val = ctk.CTkLabel(self.frame_controls, text="10 dB")
        self.lbl_tick_val.pack(pady=(0, 5))
        self.slider_tick.configure(command=self._tick_slider_event_handler)

        # Calculate Button
        self.btn_calc = ctk.CTkButton(self.frame_controls, text="CALCULATE PATTERN", 
                                      command=self.update_plot,
                                      height=40,
                                      fg_color="#1f6aa5")
        self.btn_calc.pack(pady=30)
        
        self.lbl_status = ctk.CTkLabel(self.frame_controls, text="Ready.", text_color="gray")
        self.lbl_status.pack(side="bottom", pady=10)

        # PLOT FRAME
        self.frame_plot = ctk.CTkFrame(self)
        self.frame_plot.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.update_plot()

    def _slider_event_handler(self, value):
        """Handles slider value changes with a debounce timer to update the plot."""
        self.lbl_range_val.configure(text=f"{int(value)} dB")
        
        if self.update_timer is not None:
            self.after_cancel(self.update_timer)
            
        self.update_timer = self.after(400, self.update_plot)

    def _tick_slider_event_handler(self, value):
        """Handles tick slider value changes with a debounce timer to update the plot."""
        self.lbl_tick_val.configure(text=f"{int(value)} dB")
        
        if self.update_timer is not None:
            self.after_cancel(self.update_timer)
            
        self.update_timer = self.after(400, self.update_plot)

    def update_plot(self):
        try:
            # Get simple numeric inputs
            N = int(self.entry_n.get())
            d = float(self.entry_d.get())
            beta = float(self.entry_beta.get())
            el_type = self.combo_type.get()
            view = self.combo_view.get()
            dyn_range = int(self.slider_range.get())
            tick_step = int(self.slider_tick.get())
            plot_type = self.seg_plot_type.get()

            # Parse Currents Input (CSV string)
            raw_currents = self.entry_currents.get()
            
            # Convert string "1, 0.5" -> list [1.0, 0.5]
            if not raw_currents.strip(): 
                current_list = [1.0] # Fallback if empty
            else:
                current_list = [float(x) for x in raw_currents.split(',')]
            
            # Logic: If user enters 1 value, broadcast to N. Else, require length N.
            if len(current_list) == 1:
                currents = np.full(N, current_list[0])
            elif len(current_list) == N:
                currents = np.array(current_list)
            else:
                raise ValueError(f"Currents count ({len(current_list)}) must match N ({N}) or be a single value.")

            # Calculate pattern
            theta, total_linear = self.calculator.calculate_pattern(N, d, beta, el_type, currents, view=view)

            # Convert to dB
            af_db = self.calculator.convert_to_db(total_linear, dynamic_range=dyn_range)

            # Plot Config
            self.fig.clear()

            if plot_type == "Polar":
                self.ax = self.fig.add_subplot(111, projection='polar')
                self.ax.set_theta_zero_location('N')
                self.ax.set_theta_direction(-1)
                
                angles_deg = np.arange(0, 360, 30)
                labels = [f'{d}°' if d <= 180 else f'{d - 360}°' for d in angles_deg]
                if 180 in angles_deg:
                    labels[list(angles_deg).index(180)] = '±180°'
                self.ax.set_thetagrids(angles_deg, labels)
                
                self.ax.set_ylim(-dyn_range, 0)
                step = tick_step
                self.ax.set_yticks(np.arange(-dyn_range, 1, step))
                self.ax.plot(theta, af_db, color='#1f77b4', linewidth=2)
            
            else: # Cartesian
                self.ax = self.fig.add_subplot(111)
                theta_deg = np.rad2deg(theta)
                theta_deg_shifted = np.copy(theta_deg)
                theta_deg_shifted[theta_deg_shifted >= 180] -= 360
                
                sort_indices = np.argsort(theta_deg_shifted)
                theta_deg_sorted = theta_deg_shifted[sort_indices]
                af_db_sorted = af_db[sort_indices]
                
                self.ax.plot(theta_deg_sorted, af_db_sorted, color='#1f77b4', linewidth=2)
                self.ax.set_xlabel("Angle (°)")
                self.ax.set_ylabel("Normalized Power (dB)")
                self.ax.set_xlim(-180, 180)
                self.ax.set_ylim(-dyn_range, 0)
                self.ax.set_yticks(np.arange(-dyn_range, 1, tick_step))

            # Update title
            title_text = f"Pattern ({'Horizontal' if 'Horizontal' in view else 'Vertical'}): {el_type}\nN={N}, d={d}λ, β={beta}°"
            self.ax.set_title(title_text, va='bottom', fontsize=10)
            self.ax.grid(True, alpha=0.5)
            
            self.canvas.draw()
            self.lbl_status.configure(text="Calculation successful.", text_color="green")

        except ValueError as ve:
            self.lbl_status.configure(text=f"Error: {str(ve)}", text_color="red")
        except Exception as e:
            self.lbl_status.configure(text=f"Error: {str(e)}", text_color="red")