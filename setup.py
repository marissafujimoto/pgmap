from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="pgmap",
    version="0.0.3",
    description="""pgmap is a package to help count paired guide RNAs from double
      CRISPR KO experiments""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Marissa Fujimoto and Candace Savonen",
    packages=find_packages(include=["pgmap"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: CC-BY License",
        "Operating System :: OS Independent",
    ],
    install_requires=required,
    license="CC-BY",
    url="https://github.com/FredHutch/pgmap",
)
