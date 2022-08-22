from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")


setup(
    name="ClusterOpt",
    version="0.0.0",  
    description="Atomic nanocluster optimisation using basin hopping",
    long_description=long_description,  
    long_description_content_type="text/markdown",
    url="https://github.com/sbanik2/ClusterOpt",
    author="Suvo Banik",
    author_email="sbanik2@uic.edu", 
    keywords="Nano Cluster, Basin Hopping, Optimization",
    packages=find_packages(),  
    python_requires=">=3.5, <4",
    install_requires=[
        "numpy>=1.20.1",
        "requests",
        "scipy>=1.5.0",
        "pandas",
    ],
)
