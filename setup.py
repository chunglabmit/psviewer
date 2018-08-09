from setuptools import setup

version = "0.1.0"

with open("./README.md") as fd:
    long_description = fd.read()

setup(
    name="psviewer",
    version=version,
    description=
    "Pre-stitching image viewer",
    long_description=long_description,
    install_requires=[
        "matplotlib",
        "scikit-image",
    ],
    author="Kwanghun Chung Lab",
    packages=["psviewer"],
    entry_points={ 'console_scripts': [
        'psviewer=psviewer.main:main',
    ]},
    url="https://github.com/chunglabmit/psviewer",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        'Programming Language :: Python :: 3.6',
    ]
)
