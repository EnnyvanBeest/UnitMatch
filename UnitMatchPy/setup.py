from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    description = f.read()

setup(
    name="UnitMatchPy",
    version="1.2",
    description="Description of the package",
    author="Enny van Beest, Celian Bimbard and Sam Dodgson",
    author_email="e.beest@ucl.ac.uk",
    packages=find_packages(),  
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "mtscomp",
        "jobLib",
        "tk",
        "tqdm",
    ],
    long_description=description,
    long_description_content_type='text/markdown'
)