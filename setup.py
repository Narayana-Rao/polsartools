import os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
import sys

class CMakeBuild(build_ext):
    def run(self):
        # Ensure the build directory exists
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        
        # Run CMake to configure the build
        subprocess.check_call(['cmake', '..'], cwd=self.build_temp)
        
        # Run the build command
        subprocess.check_call(['cmake', '--build', '.'], cwd=self.build_temp)
        
        super().run()

class CustomInstall(install):
    def run(self):
        # Call the install command
        install.run(self)
        # Optionally run additional post-installation commands here

setup(
    name='polsartools',
    version='0.1.0',
    description='A package for processing geospatial raster data with PolSAR tools.',
    author='Narayana',
    author_email='your.email@example.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'xarray',
        'gdal',
        'scipy',
        'pillow',  # For image processing if needed
    ],
    extras_require={
        'dev': [
            'pytest',
            'sphinx',
            'doxygen',
        ],
    },
    entry_points={
        'console_scripts': [
            'polsartools-cli=polsartools.io.main:main',  # If you have a CLI
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    cmdclass={
        'build_ext': CMakeBuild,
        'install': CustomInstall,
    },
    ext_modules=[],
)
