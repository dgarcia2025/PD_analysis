import subprocess
import sys
from pathlib import Path
import os
import tkinter as tk

class BotonSelector:
    def __init__(self):
        self.selected = None

        self.window = tk.Tk()
        self.window.title("Probelm type")
        self.window.configure(bg='#2e2e2e')
        self.window.geometry('300x150')

        self.button1 = tk.Button(
            self.window,
            text="Boundaries not grounded",
            command=lambda: self.on_button_click("Neumann"),
            bg="#444444", fg="white",
            font=("Arial", 14)
        )
        self.button1.pack(fill='x', padx=20, pady=(30, 10))

        self.button2 = tk.Button(
            self.window,
            text="Boundaries grounded",
            command=lambda: self.on_button_click("Dirichlet"),
            bg="#292626", fg="white",
            font=("Arial", 14)
        )
        self.button2.pack(fill='x', padx=20, pady=(0, 10))

        self.window.mainloop()

    def on_button_click(self, value):
        self.selected = value
        print(f"Selected: {value}")
        self.window.after(1000, self.window.destroy)

# Get the base directory of the current script
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths of the scripts and excecutables
data_prepare_file = os.path.join(BASE_DIR, "data_prepare.py")
visualization_script = os.path.join(BASE_DIR, "solution_display.py")

solver = BotonSelector()

# Executing the data preparation script
try:
    subprocess.run([sys.executable, data_prepare_file, solver.selected], check=True)
except subprocess.CalledProcessError as e:
    print(f"Error executing {data_prepare_file}: {e}")
    
# Executing the visualization script
try:
    subprocess.run([sys.executable, str(visualization_script)], check=True)
except subprocess.CalledProcessError as e:
    print(f"Error executing {visualization_script.name}: {e}")