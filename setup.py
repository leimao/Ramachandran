import setuptools

with open("README.md", "r") as fhand:
  long_description = fhand.read()

setuptools.setup(
    name="ramachandran",
    version="0.0.1",
    author="Lei Mao",
    author_email="dukeleimao@gmail.com",
    description="Generating Ramachandran plot and statistics.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/leimao/Ramachandran",
    #packages=setuptools.find_packages(),
    packages=["ramachandran"],
    package_dir={"ramachandran": "ramachandran"},
    package_data={"ramachandran": ["data/*.npz"]},
    #package_data={"ramachandran": ["data/*"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)
