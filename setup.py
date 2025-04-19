import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="demultiplexer2",
    version="1.1.5",
    author="Dominik Buchner",
    author_email="dominik.buchner@uni-due.de",
    description="A python command line interface to demultiplex illumina reads.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DominikBuchner/demultiplexer2",
    packages=setuptools.find_packages(),
    license="MIT",
    install_requires=[
        "psutil >= 5.7.3",
        "biopython >= 1.84",
        "joblib >= 0.16.0",
        "luddite >= 1.0.4",
        "pandas >= 2.2.3",
        "numpy>=1.24.0, <2.0.0",
        "tqdm>=4.66.4",
        "openpyxl>=3.1.1",
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "demultiplexer2 = demultiplexer2.__main__:main",
        ]
    },
)
