import numpy as np
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Load voltage data
V = np.load(os.path.join(BASE_DIR, 'solution', "phi.npy"))
print("Loaded voltage data with shape:", V.shape)

# Transpose for different slicing axes
V_x = np.transpose(V, (2, 1, 0))  # Shape: (x, y, z)
V_y = np.transpose(V, (1, 2, 0))  # Shape: (y, x, z)
V_z = V  # Original shape: (z, y, x)

class ElectricFieldVisualizer:
    def __init__(self, voltage_array, axis, axis_names):
        """
        Initialize visualizer with 3D voltage array
        
        Parameters:
        voltage_array: 3D numpy array of voltage values
        axis: Slicing axis ('x', 'y', or 'z')
        axis_names: String representing spatial order (e.g., 'xyz')
        """
        self.voltage = voltage_array
        self.n0, self.n1, self.n2 = voltage_array.shape
        self.axis = axis
        self.axis_names = axis_names
        
        # Create coordinates
        self.coords0 = np.arange(self.n0)
        self.coords1 = np.arange(self.n1)
        self.coords2 = np.arange(self.n2)
        
        # Calculate electric field
        self._calculate_electric_field()

    def _calculate_electric_field(self):
        """Calculate electric field E = -∇V"""
        # Compute gradient components
        grad0, grad1, grad2 = np.gradient(self.voltage)
        
        # Map gradients to spatial components
        self.Ex = np.zeros_like(self.voltage)
        self.Ey = np.zeros_like(self.voltage)
        self.Ez = np.zeros_like(self.voltage)
        
        for i, ax in enumerate(self.axis_names):
            if ax == 'x':
                self.Ex = -grad0 if i == 0 else (-grad1 if i == 1 else -grad2)
            elif ax == 'y':
                self.Ey = -grad0 if i == 0 else (-grad1 if i == 1 else -grad2)
            elif ax == 'z':
                self.Ez = -grad0 if i == 0 else (-grad1 if i == 1 else -grad2)
        
        # Calculate field magnitude
        self.E_magnitude = np.sqrt(self.Ex**2 + self.Ey**2 + self.Ez**2)

    def create_2d_slice_visualization(self):
        """Create 2D visualization with slider for planar slices"""
        # Create figure with subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                f'Voltage (V) in {self.axis_names[1]}{self.axis_names[2]} plane',
                f'Electric Field Magnitude (|E|)',
                f'Electric Field Vectors',
                f'Equipotential Lines'
            ),
            specs=[
                [{"secondary_y": False}, {"secondary_y": False}],
                [{"secondary_y": False}, {"secondary_y": False}]
            ],
            horizontal_spacing=0.1,
            vertical_spacing=0.1
        )
        
        # Initial slice index (middle of the axis)
        slice_idx = self.n0 // 2
        
        # Create initial slice data
        voltage_slice = self.voltage[slice_idx, :, :]
        e_mag_slice = self.E_magnitude[slice_idx, :, :]
        ex_slice = self.Ex[slice_idx, :, :]
        ey_slice = self.Ey[slice_idx, :, :]
        
        # 1. Voltage heatmap
        fig.add_trace(
            go.Heatmap(
                z=voltage_slice,
                x=self.coords1,
                y=self.coords2,
                colorscale='RdBu',
                name='Voltage',
                colorbar=dict(title="V", len=0.45, x=0.45)
            ),
            row=1, col=1
        )
        
        # 2. Electric field magnitude heatmap
        fig.add_trace(
            go.Heatmap(
                z=e_mag_slice,
                x=self.coords1,
                y=self.coords2,
                colorscale='Plasma',
                name='|E|',
                colorbar=dict(title="|E|", len=0.45, x=1.0)
            ),
            row=1, col=2
        )
        
        # 3. Electric field vectors (downsampled)
        sample = max(1, self.n1//20, self.n2//20)  # Adaptive sampling
        x_vec = self.coords1[::sample]
        y_vec = self.coords2[::sample]
        
        # Create vector field
        for i in range(0, len(x_vec)):
            for j in range(0, len(y_vec)):
                xi = np.searchsorted(self.coords1, x_vec[i])
                yj = np.searchsorted(self.coords2, y_vec[j])
                
                if xi < ex_slice.shape[0] and yj < ex_slice.shape[1]:
                    ex_val = ex_slice[xi, yj]
                    ey_val = ey_slice[xi, yj]
                    magnitude = np.sqrt(ex_val**2 + ey_val**2)
                    
                    if magnitude > 0:
                        # Scale vectors for visibility
                        scale = 0.5 * sample
                        dx = ex_val / magnitude * scale
                        dy = ey_val / magnitude * scale
                        
                        fig.add_trace(
                            go.Scatter(
                                x=[x_vec[i], x_vec[i] + dx],
                                y=[y_vec[j], y_vec[j] + dy],
                                mode='lines',
                                line=dict(color='red', width=1.5),
                                showlegend=False
                            ),
                            row=2, col=1
                        )
        
        # 4. Equipotential lines
        fig.add_trace(
            go.Contour(
                z=voltage_slice,
                x=self.coords1,
                y=self.coords2,
                contours=dict(
                    start=np.min(voltage_slice),
                    end=np.max(voltage_slice),
                    size=(np.max(voltage_slice) - np.min(voltage_slice)) / 10
                ),
                line=dict(width=1.5),
                colorscale='RdBu',
                name='Equipotential',
                showscale=False
            ),
            row=2, col=2
        )
        
        # Create slider steps
        steps = []
        for i in range(self.n0):
            step = dict(
                method='update',
                args=[{
                    'z': [
                        self.voltage[i, :, :],  # Voltage
                        self.E_magnitude[i, :, :],  # |E|
                        self.voltage[i, :, :]   # Equipotential (same as voltage)
                    ],
                    'x': [None, None, x_vec, None],
                    'y': [None, None, y_vec, None]
                }],
                label=f'{i}'
            )
            steps.append(step)
        
        sliders = [dict(
            active=slice_idx,
            currentvalue=dict(prefix=f"{self.axis} = "),
            pad=dict(t=30),
            steps=steps
        )]
        
        # Configure layout
        fig.update_layout(
            title_text=f"Electric Field Analysis - Slicing along {self.axis}-axis",
            height=800,
            width=1200,
            margin=dict(t=100, b=50, l=50, r=50),
            sliders=sliders
        )
        
        # Set axis labels
        x_label, y_label = self.axis_names[1], self.axis_names[2]
        fig.update_xaxes(title_text=x_label, row=1, col=1)
        fig.update_yaxes(title_text=y_label, row=1, col=1)
        fig.update_xaxes(title_text=x_label, row=1, col=2)
        fig.update_yaxes(title_text=y_label, row=1, col=2)
        fig.update_xaxes(title_text=x_label, row=2, col=1)
        fig.update_yaxes(title_text=y_label, row=2, col=1)
        fig.update_xaxes(title_text=x_label, row=2, col=2)
        fig.update_yaxes(title_text=y_label, row=2, col=2)
        
        return fig

def generate_visualizations():
    """Generate all 2D slice visualizations"""
    # Z-axis slices (original orientation)
    print("Creating visualizations for Z-axis slices...")
    viz_z = ElectricFieldVisualizer(V_z, 'z', 'zyx')
    fig_z = viz_z.create_2d_slice_visualization()
    fig_z.write_html(os.path.join(BASE_DIR, 'solution', 'electric_field_z_slices.html'))
    
    # X-axis slices
    print("Creating visualizations for X-axis slices...")
    viz_x = ElectricFieldVisualizer(V_x, 'x', 'xyz')
    fig_x = viz_x.create_2d_slice_visualization()
    fig_x.write_html(os.path.join(BASE_DIR, 'solution', 'electric_field_x_slices.html'))
    
    # Y-axis slices
    print("Creating visualizations for Y-axis slices...")
    viz_y = ElectricFieldVisualizer(V_y, 'y', 'yxz')
    fig_y = viz_y.create_2d_slice_visualization()
    fig_y.write_html(os.path.join(BASE_DIR, 'solution', 'electric_field_y_slices.html'))
    
    print("All visualizations saved as HTML files")

# Generate all visualizations
generate_visualizations()