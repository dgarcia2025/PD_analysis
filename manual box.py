import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import trimesh
import tkinter as tk
from tkinter import filedialog, messagebox
import sys
import os

pio.renderers.default = "browser"

def load_and_process_stl(file_path):
    try:
        mesh = trimesh.load(file_path)
        if hasattr(mesh, 'split'):
            split_meshes = mesh.split(only_watertight=False)
            print(f"STL file split into {len(split_meshes)} components")
            return split_meshes
        else:
            return [mesh]
    except Exception as e:
        messagebox.showerror("STL Load Error", str(e))
        return None

def create_cube_wireframe(x_center=0, y_center=0, z_center=0, width=1, height=1, depth=1):
    vertices = np.array([
        [x_center - width/2, y_center - height/2, z_center - depth/2],
        [x_center + width/2, y_center - height/2, z_center - depth/2],
        [x_center + width/2, y_center + height/2, z_center - depth/2],
        [x_center - width/2, y_center + height/2, z_center - depth/2],
        [x_center - width/2, y_center - height/2, z_center + depth/2],
        [x_center + width/2, y_center - height/2, z_center + depth/2],
        [x_center + width/2, y_center + height/2, z_center + depth/2],
        [x_center - width/2, y_center + height/2, z_center + depth/2]
    ])
    edges = np.array([
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7]
    ])
    return vertices, edges

def create_plotly_wireframe(vertices, edges, name="Wireframe", color="red", line_width=3):
    x_lines, y_lines, z_lines = [], [], []
    for edge in edges:
        start, end = vertices[edge[0]], vertices[edge[1]]
        x_lines += [start[0], end[0], None]
        y_lines += [start[1], end[1], None]
        z_lines += [start[2], end[2], None]
    return go.Scatter3d(x=x_lines, y=y_lines, z=z_lines,
                        mode='lines', line=dict(color=color, width=line_width),
                        name=name, showlegend=True)

def create_plotly_mesh(vertices, faces, name="Mesh", color="lightblue", opacity=0.7):
    return go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        name=name, color=color, opacity=opacity, showscale=False
    )

class STLViewerApp:
    def __init__(self, root, stl_path):
        self.root = root
        self.root.title("STL Viewer")
        self.root.configure(bg='#313354')
        self.root.geometry("400x450")
        self.meshes = None
        self.cube_params = None
        self.step_size = 1.0
        self.stl_path = stl_path

        self.style = {
            'bg': '#313354',
            'fg': '#FFFFFF',
            'button_bg': '#4A4A7A',
            'button_fg': '#FFFFFF',
            'entry_bg': '#2A2A4A',
            'entry_fg': '#FFFFFF',
            'font': ('Arial', 10)
        }

        self._create_gui()
        self.load_file(path=self.stl_path)  # Load STL on start

    def _create_gui(self):
        main_frame = tk.Frame(self.root, bg=self.style['bg'], padx=20, pady=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        title = tk.Label(main_frame, text="STL Viewer", font=('Arial', 14, 'bold'),
                         bg=self.style['bg'], fg=self.style['fg'])
        title.pack(pady=10)

        self.fields = {}
        input_frame = tk.Frame(main_frame, bg=self.style['bg'])
        input_frame.pack(fill=tk.X, pady=10)

        for label in ["X Center", "Y Center", "Z Center",
                      "Width (x_axis)", "Height (y_axis)", "Depth (z_axis)",
                      "Step Size"]:
            row = tk.Frame(input_frame, bg=self.style['bg'])
            row.pack(fill=tk.X, pady=5)
            lab = tk.Label(row, width=15, text=label, anchor='w',
                           bg=self.style['bg'], fg=self.style['fg'],
                           font=self.style['font'])
            ent = tk.Entry(row, bg=self.style['entry_bg'], fg=self.style['entry_fg'],
                           insertbackground=self.style['fg'], font=self.style['font'])
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            self.fields[label] = ent

        button_frame = tk.Frame(main_frame, bg=self.style['bg'])
        button_frame.pack(fill=tk.X, pady=15)

        update_btn = tk.Button(button_frame,
                  text="Update Cube",
                  command=self.update_plot,
                  bg=self.style['button_bg'], fg=self.style['button_fg'],
                  font=self.style['font'], relief=tk.FLAT,
                  padx=10, pady=5)

        select_btn = tk.Button(button_frame,
                  text="Select Cube",
                  command=self.select_cube,
                  bg=self.style['button_bg'], fg=self.style['button_fg'],
                  font=self.style['font'], relief=tk.FLAT,
                  padx=10, pady=5)

        update_btn.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=10)
        select_btn.pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=10)

    def load_file(self, path):
        if not path or not os.path.isfile(path):
            messagebox.showerror("STL Error", "Invalid STL path")
            return

        self.meshes = load_and_process_stl(path)
        if not self.meshes:
            return

        all_verts = np.vstack([m.vertices for m in self.meshes])
        xc, yc, zc = all_verts[:, 0].mean(), all_verts[:, 1].mean(), all_verts[:, 2].mean()

        self.fields["X Center"].insert(0, f"{xc:.2f}")
        self.fields["Y Center"].insert(0, f"{yc:.2f}")
        self.fields["Z Center"].insert(0, f"{zc:.2f}")
        for dim in ["Width (x_axis)", "Height (y_axis)", "Depth (z_axis)"]:
            self.fields[dim].insert(0, "1.0")
        self.fields["Step Size"].insert(0, "1.0")

        self.update_plot()

    def update_plot(self):
        if not self.meshes:
            return

        try:
            x = float(self.fields["X Center"].get())
            y = float(self.fields["Y Center"].get())
            z = float(self.fields["Z Center"].get())
            w = float(self.fields["Width (x_axis)"].get())
            h = float(self.fields["Height (y_axis)"].get())
            d = float(self.fields["Depth (z_axis)"].get())
            step = float(self.fields["Step Size"].get())
            self.cube_params = {'x': x, 'y': y, 'z': z, 'width': w, 'height': h, 'depth': d}
            self.step_size = step
        except ValueError:
            messagebox.showerror("Input Error", "Invalid float input")
            return

        fig = go.Figure()
        colors = px.colors.qualitative.Set3
        for i, m in enumerate(self.meshes):
            fig.add_trace(create_plotly_mesh(m.vertices, m.faces, name=f"Part {i+1}",
                                             color=colors[i % len(colors)], opacity=0.6))

        cube_vertices, cube_edges = create_cube_wireframe(x, y, z, w, h, d)
        fig.add_trace(create_plotly_wireframe(cube_vertices, cube_edges,
                                              name="Cube", color="red", line_width=4))

        fig.update_layout(title="STL Viewer", scene=dict(aspectmode='data'), width=800, height=600)
        fig.show()

    def select_cube(self):
        try:
            x = float(self.fields["X Center"].get())
            y = float(self.fields["Y Center"].get())
            z = float(self.fields["Z Center"].get())
            w = float(self.fields["Width (x_axis)"].get())
            h = float(self.fields["Height (y_axis)"].get())
            d = float(self.fields["Depth (z_axis)"].get())
            step = float(self.fields["Step Size"].get())
            self.cube_params = {'x': x, 'y': y, 'z': z, 'width': w, 'height': h, 'depth': d, 'step': step}
            self.step_size = step
            print(f"CUBE_PARAMS:{self.cube_params}")
            # print(f"CUBE_PARAMS:{self.cube_params}\nSTEP:{self.step_size}")
            self.root.destroy()
            sys.exit(0)
        except ValueError:
            messagebox.showerror("Input Error", "Invalid float input")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <path_to_stl>")
        sys.exit(1)

    stl_path = sys.argv[1]
    root = tk.Tk()
    app = STLViewerApp(root, stl_path=stl_path)
    root.mainloop()