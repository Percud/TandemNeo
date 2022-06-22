import setuptools
import os

setuptools.setup(
    name="TandemNeo",
    version="0.0.1",
    author="Carlo De Rito, Marco Malatesta",
    author_email="carlo.derito2@gmail.com, marcomala46@gmail.com",
    description="Neofunctionalization of duplicated genes",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Percud/TandemNeo",
    packages=setuptools.find_packages(),
    classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: Unix",
            "Development Status::4 - Beta",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            #"License :: OSI Approved :: MIT License",
            ],
    python_requires='>=3.7',
    install_requires=[
                    'natsort==7.1.0', 
                    'numpy==1.22.0', 
                    'wget==3.2', 
                    'biopandas==0.2.7', 
                    'requests==2.25.1', 
                    'pandas==1.1.5', 
                    'Bio==0.3.0', 
                    'scikit_learn==0.24.1'
                    ],
    entry_points={
                'console_scripts': [
                'TandemNeo=tandemneo.TandemNeo:main',
                ],},
    )



