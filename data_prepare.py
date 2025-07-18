import trimesh
import numpy as np
import os
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.font as tkFont
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from tkinter import ttk, filedialog, messagebox
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.figure import Figure
import subprocess
import sys  
from math import ceil

problem_type = sys.argv[1]

def get_file_path():
    # root = tk.Tk()
    # root.withdraw()  # Hide root window
    file_path = filedialog.askopenfilename(title="Select a file")
    return file_path

root = tk.Tk()
root.withdraw()  # Ocultar ventana raíz

# Usage
file = get_file_path()

def show_plot_with_voltage_input(root, solid, solids, i):
    voltage_result = {'value': None}

    def on_submit():
        try:
            v = float(entry.get())
            voltage_result['value'] = v
            window.destroy()
        except ValueError:
            messagebox.showerror("Entrada inválida", "Por favor, ingrese un número válido.")

    window = tk.Toplevel(root)
    window.title(f"Sólido {i+1} - Visualización y Voltaje")
    window.configure(bg="#2b2b2b")

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(solid.vertices[:, 0], solid.vertices[:, 1], solid.vertices[:, 2], color='r', s=1)

    for j, other_solid in enumerate(solids):
        if j != i:
            ax.scatter(other_solid.vertices[:, 0], other_solid.vertices[:, 1], other_solid.vertices[:, 2], color='b', s=1)

    all_vertices = np.vstack([s.vertices for s in solids])
    x_min, y_min, z_min = all_vertices.min(axis=0)
    x_max, y_max, z_max = all_vertices.max(axis=0)

    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min) / 2

    ax.set_xlim(x_mid - max_range, x_mid + max_range)
    ax.set_ylim(y_mid - max_range, y_mid + max_range)
    ax.set_zlim(z_mid - max_range, z_mid + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Sólido {i+1} en rojo, demás en azul')

    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, rowspan=4, padx=10, pady=10)

    label = tk.Label(window, text=f"Ingrese voltaje para el sólido {i+1}:", font=("Segoe UI", 10), bg="#2e2c2c", fg="#ffffff")
    label.grid(row=0, column=1, padx=10, pady=(10, 5), sticky='w')

    entry = tk.Entry(window, font=("Segoe UI", 10), bg="#3c3f41", fg="#ffffff", insertbackground="white")
    entry.grid(row=1, column=1, padx=10, pady=5, sticky='w')

    submit_btn = tk.Button(window, text="Aceptar", command=on_submit,
                           bg="#313354", fg="white", font=("Segoe UI", 10, "bold"),
                           activebackground="#3d3f81", relief="flat")
    submit_btn.grid(row=2, column=1, padx=10, pady=5, sticky='w')

    window.grab_set()
    window.wait_window()
    return voltage_result['value']

# Carpeta donde está el script actual
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Cargar el archivo STL
# Load original STL file
# file = r"C:\Users\dgara\Documents\AutoCAD\Demo\Drawing3.stl"
mesh = trimesh.load(file)  # Reemplaza con la ruta de tu archivo

# Separar en sólidos individuales
solids = mesh.split(only_watertight=False)

# Asignar voltajes constantes
voltages = []
for i, solid in enumerate(solids):
    print(f"Procesando sólido {i+1}/{len(solids)}...")

    # show_plot_with_voltage_input(solid, solids, i, voltages)
    v = show_plot_with_voltage_input(root, solid, solids, i)
    voltages.append(float(v))

root.destroy()  # Cerrar ventana raíz después de editar el bounding box
root = tk.Tk()
# root.withdraw()  # Ocultar ventana raíz

def ask_grid_resolution(parent=None):
    """
    Creates a professional dialog asking whether grid should be automatic or manual
    Args: parent - parent window (optional, uses root if None)
    Returns: 'automatic', 'manual', or None if cancelled
    """
    result = None
    
    # Create the dialog window
    dialog = tk.Toplevel(parent)
    dialog.title("Grid Resolution Settings")
    dialog.geometry("400x250")
    dialog.configure(bg="#313354")
    dialog.resizable(False, False)
    
    # Center the dialog relative to parent or screen
    dialog.update_idletasks()
    if parent:
        x = parent.winfo_x() + (parent.winfo_width() // 2) - (dialog.winfo_width() // 2)
        y = parent.winfo_y() + (parent.winfo_height() // 2) - (dialog.winfo_height() // 2)
    else:
        x = (dialog.winfo_screenwidth() // 2) - (dialog.winfo_width() // 2)
        y = (dialog.winfo_screenheight() // 2) - (dialog.winfo_height() // 2)
    dialog.geometry(f"+{x}+{y}")
    
    # Configure custom styles
    style = ttk.Style()
    style.theme_use('clam')
    
    # Configure button styles
    style.configure('Custom.TButton',
                   background='#4a4d6b',
                   foreground='white',
                   borderwidth=0,
                   focuscolor='none',
                   font=('Arial', 10, 'bold'))
    
    style.map('Custom.TButton',
              background=[('active', '#5a5d7b'),
                         ('pressed', '#3a3d5b')])
    
    # Configure automatic button style (highlighted)
    style.configure('Auto.TButton',
                   background='#6c7ae0',
                   foreground='white',
                   borderwidth=0,
                   focuscolor='none',
                   font=('Arial', 10, 'bold'))
    
    style.map('Auto.TButton',
              background=[('active', '#7c8af0'),
                         ('pressed', '#5c6ad0')])
    
    # Main frame
    main_frame = tk.Frame(dialog, bg="#313354", padx=30, pady=30)
    main_frame.pack(fill=tk.BOTH, expand=True)
    
    # Title
    title_font = tkFont.Font(family="Arial", size=16, weight="bold")
    title_label = tk.Label(main_frame, 
                          text="Grid Resolution", 
                          font=title_font,
                          bg="#313354", 
                          fg="white")
    title_label.pack(pady=(0, 15))
    
    # Question text
    question_font = tkFont.Font(family="Arial", size=11)
    question_label = tk.Label(main_frame, 
                             text="How would you like the grid to be configured?",
                             font=question_font,
                             bg="#313354", 
                             fg="#c0c0c0",
                             wraplength=320)
    question_label.pack(pady=(0, 25))
    
    # Button frame
    button_frame = tk.Frame(main_frame, bg="#313354")
    button_frame.pack(pady=(0, 10))
    
    def on_automatic():
        nonlocal result
        result = 'automatic'
        dialog.destroy()
    
    def on_manual():
        nonlocal result
        result = 'manual'
        dialog.destroy()
    
    def on_cancel():
        nonlocal result
        result = None
        dialog.destroy()
    
    # Automatic button (highlighted as recommended)
    auto_btn = ttk.Button(button_frame, 
                         text="Automatic",
                         style='Auto.TButton',
                         command=on_automatic,
                         width=12)
    auto_btn.pack(side=tk.LEFT, padx=(0, 10))
    
    # Manual button
    manual_btn = ttk.Button(button_frame, 
                           text="Manual",
                           style='Custom.TButton',
                           command=on_manual,
                           width=12)
    manual_btn.pack(side=tk.LEFT, padx=(10, 0))
    
    # Description text
    desc_font = tkFont.Font(family="Arial", size=9)
    desc_label = tk.Label(main_frame, 
                         text="Automatic: Grid adjusts to fit solids automatically\nManual: You control grid size and positioning",
                         font=desc_font,
                         bg="#313354", 
                         fg="#9090a0",
                         justify=tk.CENTER)
    desc_label.pack(pady=(15, 0))
    
    # Cancel button (small, at bottom)
    cancel_btn = ttk.Button(main_frame, 
                           text="Cancel",
                           style='Custom.TButton',
                           command=on_cancel,
                           width=8)
    cancel_btn.pack(pady=(15, 0))
    
    # Handle window close button
    dialog.protocol("WM_DELETE_WINDOW", on_cancel)
    
    # Set focus and grab - make dialog modal
    dialog.transient(parent)
    dialog.grab_set()
    dialog.focus_set()
    
    # Wait for dialog to close
    dialog.wait_window(dialog)
    
    return result

grid_type = ask_grid_resolution(root)

root.destroy()  # Cerrar ventana raíz después de editar el bounding box

print('Generating grid...')

if grid_type == 'automatic':
    # Calcular el bounding box
    bounds = mesh.bounds
    min_corner = bounds[0]
    max_corner = bounds[1]

elif grid_type == 'manual':
    import ast
    box_py_path = os.path.join(BASE_DIR, 'manual box.py')
    # subprocess.run([sys.executable, box_py_path], check=True)
    # Run the script and capture output
    result = subprocess.run(
        [sys.executable, box_py_path, file],
        check=True,
        capture_output=True,  # Capture stdout and stderr
        text=True  # Return output as string instead of bytes
    )

    print(result.stdout)  # Print the output of the script
    
    # Check stdout for the cube parameters
    output = result.stdout
    cube_params = None
    
    # Look for the line starting with "CUBE_PARAMS:"
    for line in output.splitlines():
        if line.startswith("CUBE_PARAMS:"):
            # Extract the dictionary string and parse it
            params_str = line[len("CUBE_PARAMS:"):]  # Remove prefix
            try:
                cube_params = ast.literal_eval(params_str)  # Safely parse string to dict
                break
            except (ValueError, SyntaxError) as e:
                print(f"Error parsing cube parameters: {e}")
                cube_params = None
    
    if cube_params:
        # print("Cube Parameters:", cube_params)
        # Access individual parameters
        x = cube_params.get('x')
        y = cube_params.get('y')
        z = cube_params.get('z')
        width = cube_params.get('width')
        height = cube_params.get('height')
        depth = cube_params.get('depth')
        step = cube_params.get('step')  # Default step size if not provided
        # print(f"X: {x}, Y: {y}, Z: {z}, Width: {width}, Height: {height}, Depth: {depth}")
        min_corner = (x - width/2, y - height/2, z - depth/2)
        max_corner = (x + width/2, y + height/2, z + depth/2)
        print(f"Min Corner: {min_corner}, Max Corner: {max_corner}")
    else:
        print("No cube parameters found in output")

# Calculate cube extents
dx = max_corner[0] - min_corner[0]  # Extent along x
dy = max_corner[1] - min_corner[1]  # Extent along y
dz = max_corner[2] - min_corner[2]  # Extent along z

# Define target resolution (step size)
step_size = step # ref_step  # Maintain this resolution across all dimensions

# Calculate number of points for each dimension
Nx = ceil(dx / step_size) # + 1  # Add 1 to include endpoint
Ny = ceil(dy / step_size) # + 1
Nz = ceil(dz / step_size) # + 1

# Create grid with different number of points per dimension
x = np.linspace(min_corner[0], max_corner[0], Nx)
y = np.linspace(min_corner[1], max_corner[1], Ny)
z = np.linspace(min_corner[2], max_corner[2], Nz)
grid = np.meshgrid(x, y, z, indexing='ij')
grid_points = np.vstack([grid[0].ravel(), grid[1].ravel(), grid[2].ravel()]).T

# Create arrays for types and Dirichlet values
type_grid = np.zeros((Nx, Ny, Nz), dtype=int)  # 0: interior, 1: Dirichlet, 2: Neumann
dirichlet_grid = np.zeros((Nx, Ny, Nz))

# Mark points inside solids (Dirichlet)
for idx, p in enumerate(grid_points):
    i, j, k = np.unravel_index(idx, (Nx, Ny, Nz))
    for s, solid in enumerate(solids):
        if solid.contains([p]):  # Verify if point is inside the solid
            type_grid[i, j, k] = 1
            dirichlet_grid[i, j, k] = voltages[s]
            break

# Mark boundary points as Neumann if not Dirichlet
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            if (i == 0 or i == Nx-1 or j == 0 or j == Ny-1 or k == 0 or k == Nz-1) and type_grid[i, j, k] == 0:
                type_grid[i, j, k] = 2

print('Saving data...')

os.makedirs(os.path.join(BASE_DIR, 'Grid_data'), exist_ok=True)
grid_folder = os.path.join(BASE_DIR, 'Grid_data')

# Save data
# Save grid dimensions (Nx, Ny, Nz)
np.savetxt(os.path.join(grid_folder, 'grid_dims.txt'), [Nx, Ny, Nz], fmt='%d', header='Nx Ny Nz')

# Save step sizes for each dimension
step_x = (max_corner[0] - min_corner[0]) / (Nx - 1) if Nx > 1 else 0
step_y = (max_corner[1] - min_corner[1]) / (Ny - 1) if Ny > 1 else 0
step_z = (max_corner[2] - min_corner[2]) / (Nz - 1) if Nz > 1 else 0
np.savetxt(os.path.join(grid_folder, 'grid_h.txt'), [step_x, step_y, step_z], 
           fmt='%.10e')

# Save type_grid and dirichlet_grid with dimension header
np.savetxt(os.path.join(grid_folder, 'type_grid.txt'), type_grid.flatten(), 
           fmt='%d', header=f'Nx={Nx} Ny={Ny} Nz={Nz}')
np.savetxt(os.path.join(grid_folder, 'dirichlet_grid.txt'), dirichlet_grid.flatten(), 
           fmt='%.10e', header=f'Nx={Nx} Ny={Ny} Nz={Nz}')

print('Guardado exitosamente!')

# Solve for dirichlet
if problem_type == 'Dirichlet':
    p = subprocess.Popen([os.path.join(BASE_DIR, './solver_dirichelet')], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

# Solve for neumann
elif problem_type == 'Neumann':
    p = subprocess.Popen([os.path.join(BASE_DIR, './solver_neumann')], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

else:
    raise ValueError(f"Unknown problem type: {problem_type}. Expected 'dirichlet' or 'neumann'.")

hx = step_x
hy = step_y
hz = step_z

# First send dims
p.stdin.write(f"{Nx} {Ny} {Nz}\n".encode())
# Then send h_x h_y h_z
p.stdin.write(f"{hx:.10e} {hy:.10e} {hz:.10e}\n".encode())
# Then send type_grid (flattened)
for t in type_grid.flatten():
    p.stdin.write(f"{t} ".encode())
p.stdin.write(b"\n")
# Then send dirichlet_grid (flattened)
for v in dirichlet_grid.flatten():
    p.stdin.write(f"{v:.10e} ".encode())
p.stdin.write(b"\n")
p.stdin.flush()

# Read back the solution (flat)
phi_data = p.stdout.read().split()
# print(phi_data) 
phi = np.array(phi_data, dtype=float).reshape((Nx,Ny,Nz))

os.makedirs(os.path.join(BASE_DIR, 'solution'), exist_ok=True)
solution_folder = os.path.join(BASE_DIR, 'solution')
np.save(os.path.join(solution_folder, 'phi.npy'), phi)