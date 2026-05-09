from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enteroscope",
    version="1.0.0",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    description="Advanced Enterobacter cloacae Complex Genomic Analysis Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/enteroscope",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'enteroscope': ['**/*'],   # include EVERYTHING under enteroscope/ folder
    },
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.85",
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "beautifulsoup4>=4.11.0",
        "lxml>=4.9.0",
    ],
    entry_points={
        "console_scripts": [
            "enteroscope=enteroscope.enteroscope:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
)