from setuptools import setup, find_packages

setup(
    name="UMPy",
    version="0.1",
    description="Description of the package",
    author="e.beest@ucl.ac.uk",
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
        "spikeinterface"
    ],
)