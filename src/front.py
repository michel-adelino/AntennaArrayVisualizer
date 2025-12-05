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
        self.geometry("1000x650")
        ctk.set_appearance_mode("Dark")
        ctk.set_default_color_theme("blue")

        # Backend instance
        self.calculator = AntennaCalculator()

        # --- MAIN LAYOUT ---
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # CONTROLS FRAME
        self.frame_controls = ctk.CTkFrame(self, width=250)
        self.frame_controls.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.lbl_title = ctk.CTkLabel(self.frame_controls, text="Parameters", font=("Arial", 20, "bold"))
        self.lbl_title.pack(pady=20)

        # Input: Number of antennas (N)
        self.lbl_n = ctk.CTkLabel(self.frame_controls, text="Number of Antennas (N):")
        self.lbl_n.pack(pady=(10, 0))
        self.entry_n = ctk.CTkEntry(self.frame_controls)
        self.entry_n.insert(0, "4")
        self.entry_n.pack(pady=5)

        # Input: Separation (d/lambda)
        self.lbl_d = ctk.CTkLabel(self.frame_controls, text="Separation (d/λ):")
        self.lbl_d.pack(pady=(10, 0))
        self.entry_d = ctk.CTkEntry(self.frame_controls)
        self.entry_d.insert(0, "0.5")
        self.entry_d.pack(pady=5)

        # --- Dynamic Range Control ---
        self.lbl_range = ctk.CTkLabel(self.frame_controls, text="Dynamic Range (dB):")
        self.lbl_range.pack(pady=(20, 0))
        
        self.slider_range = ctk.CTkSlider(self.frame_controls, from_=10, to=80, number_of_steps=70)
        self.slider_range.set(40) # Default value
        self.slider_range.pack(pady=5)
        
        # Label to show the slider's current value
        self.lbl_range_val = ctk.CTkLabel(self.frame_controls, text="40 dB")
        self.lbl_range_val.pack(pady=(0, 5))
        
        # Callback to update the label text
        self.slider_range.configure(command=lambda val: self.lbl_range_val.configure(text=f"{int(val)} dB"))

        # Calculate Button
        self.btn_calc = ctk.CTkButton(self.frame_controls, text="CALCULATE PATTERN", command=self.update_plot)
        self.btn_calc.pack(pady=30)
        
        # Status label
        self.lbl_status = ctk.CTkLabel(self.frame_controls, text="Ready.", text_color="gray")
        self.lbl_status.pack(side="bottom", pady=10)

        # PLOT FRAME
        self.frame_plot = ctk.CTkFrame(self)
        self.frame_plot.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        self.initialize_plot()

    def initialize_plot(self):
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='polar')
        self.ax.set_title("Radiation Pattern", pad=20)
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def update_plot(self):
        try:
            # Get values from UI
            N = int(self.entry_n.get())
            d = float(self.entry_d.get())
            dyn_range = int(self.slider_range.get())

            # Calculate linear pattern
            theta, af_linear = self.calculator.calculate_pattern(N, d)

            # Convert to dB using the backend function
            af_db = self.calculator.convert_to_db(af_linear, dynamic_range=dyn_range)

            # Configure plot
            self.ax.clear()
            self.ax.set_theta_zero_location('N') # 0 degrees on top
            self.ax.set_theta_direction(-1)      # Clockwise
            
            # --- Adjustments for dB view ---
            self.ax.set_ylim(-dyn_range, 0)
            
            # Adjust grid circles (every 10dB or 5dB)
            step = 10 if dyn_range > 30 else 5
            self.ax.set_yticks(np.arange(-dyn_range, 1, step))
            
            # Plot
            self.ax.plot(theta, af_db, color='#1f77b4', linewidth=2)
            self.ax.set_title(f"Pattern (dB) for N={N}, d={d}λ", va='bottom')
            self.ax.grid(True, alpha=0.5)
            
            self.canvas.draw()
            self.lbl_status.configure(text="Calculation successful.", text_color="green")

        except ValueError:
            self.lbl_status.configure(text="Error: Check numeric values.", text_color="red")
        except Exception as e:
            self.lbl_status.configure(text=f"Error: {str(e)}", text_color="red")