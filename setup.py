from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'MolVizr3D'
VERSION = '1.0.0'
DESCRIPTION = 'MolVizr3D: Molecular Visualizer in 3D.'
LONG_DESCRIPTION = 'MolVizr3D is a Python package designed to provide an intuitive and powerful solution for visualizing molecular structures in three dimensions. This package supports various molecule formats such as SDF, MOL, and more, allowing users to easily load and render complex molecular structures with high-quality 3D graphics.'

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author="Sina Gilassi",
    author_email="<sina.gilassi@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    license='MIT',
    install_requires=['pandas', 'pillow', 'requests',
                      'urllib3', 'matplotlib', 'PubChemQuery', 'numpy'],
    keywords=['python', 'chemistry', 'chemistry-visualization', 'MolVizr3D'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.6',
)
