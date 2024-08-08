# MolVizr3D: Molecular Visualizer in 3D

![Downloads](https://img.shields.io/pypi/dm/MolVizr3D) ![PyPI](https://img.shields.io/pypi/v/MolVizr3D) ![Python Version](https://img.shields.io/pypi/pyversions/MolVizr3D.svg) ![License](https://img.shields.io/pypi/l/MolVizr3D)

MolVizr3D is a Python package designed to provide an intuitive and powerful solution for visualizing molecular structures in three dimensions. This package supports various molecule formats such as SDF, MOL, and more, allowing users to easily load and render complex molecular structures with high-quality 3D graphics.

**Key Features:**

* `Wide Format Support:` Load and visualize molecular structures from various file formats including SDF, MOL, and others.
* `High-Quality 3D Rendering:` Generate detailed and accurate 3D visualizations of molecular structures.
* `Interactive Visualization:` Rotate, zoom, and explore molecules interactively to better understand their geometry and properties.
* `Customization:` Customize visual aspects such as colors, styles, and labels to suit your specific needs.
* `Integration:` Seamlessly integrates with other computational chemistry tools and workflows.

## Installation

Install MolVizr3D with pip

```python
  pip install MolVizr3D
```

## Examples

Import package as:

```python
import molvizr3d as mv3d
```

To check mv3d version:

```python
# check mv3d version
print(mv3d.__version__)
```

To visualize the structure of compound using its sdf file:

```python
# sdf file
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_7979.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# visualize compound by sdf file
mv3d.td(sdf_file)
```

To visualize the structure of compound using its InChI

```python
# InChI
inchi = 'InChI=1S/C3H8/c1-3-2/h3-4H,1H2,2H3'

# visualize compound by inchi
mv3d.td_by_inchi()
```

## FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/) 

## Authors

- [@sinagilassi](https://www.github.com/sinagilassi)

