import plotly.graph_objects as go
import numpy as np

# Define cylinder parameters
r = 1  # radius
h = 5  # height
n = 100  # number of vertices around the circumference

# Create vertices
theta = np.linspace(0, 2*np.pi, n, endpoint=False)
x = r * np.cos(theta)
y = r * np.sin(theta)
z = np.zeros_like(x)
vertices_top = np.column_stack((x, y, z+h))
vertices_bottom = np.column_stack((x, y, z))

# Create mesh
mesh = go.Mesh3d(
    x=np.concatenate((vertices_top[:, 0], vertices_bottom[:, 0])),
    y=np.concatenate((vertices_top[:, 1], vertices_bottom[:, 1])),
    z=np.concatenate((vertices_top[:, 2], vertices_bottom[:, 2])),
    i=np.concatenate((np.arange(n), np.arange(n) + n, np.arange(n))),
    j=np.concatenate((np.arange(n) + n, np.arange(n) + 2*n, np.arange(n) + n)),
    k=np.concatenate((np.arange(n) + 2*n, np.arange(n), np.arange(n) + 2*n)),
    color='black',
    opacity=1
)

# Plot cylinder
fig = go.Figure(data=[mesh])
fig.update_layout(scene=dict(
    xaxis_title='X',
    yaxis_title='Y',
    zaxis_title='Z',
    aspectratio=dict(x=1, y=1, z=0.5)),
    paper_bgcolor='white',
    margin=dict(l=0, r=0, b=0, t=0))
fig.show()
