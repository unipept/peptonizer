from setuptools import setup, find_packages

setup(
    name="peptonizer",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "networkx",
        "pandas"
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A simple example package",
    url="https://github.com/yourusername/mypackage",  # Replace with your GitHub URL
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)