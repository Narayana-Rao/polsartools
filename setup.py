import os
import subprocess
from setuptools import setup, find_packages,Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
import sys,shutil
import pybind11

class BuildExt(build_ext):
    """Custom build_ext to ensure .pyd/.so files are placed inside the package directory."""
    
    def build_extension(self, ext):
        # Get the full path where the extension is originally built
        ext_path = self.get_ext_fullpath(ext.name)
        
        # Correct the destination directory to be inside the package
        ext_name = os.path.basename(ext_path)
        dest_path = os.path.join(os.path.dirname(__file__), "polsartools", ext_name)

        # Call the original method to build the extension
        super().build_extension(ext)

        # Ensure the target directory exists
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        
        # Copy the built extension to the package directory
        shutil.copy2(ext_path, dest_path)

ext_modules = [
    Extension(
        'polsartools.refined_lee',
        ['cpp/src/refined_lee.cpp'],
        include_dirs=[pybind11.get_include(), pybind11.get_include(user=True)],
        language='c++',
        extra_compile_args=['/std:c++17' if os.name == 'nt' else '-std=c++17'],
    ),
    Extension(
        'polsartools.moving_average',
        ['cpp/src/moving_average.cpp'],
        include_dirs=[pybind11.get_include(), pybind11.get_include(user=True)],
        language='c++',
        extra_compile_args=['/std:c++17' if os.name == 'nt' else '-std=c++17'],
    ),
    Extension(
        'polsartools.sum_arrays',
        ['cpp/src/sum_arrays.cpp'],
        include_dirs=[pybind11.get_include(), pybind11.get_include(user=True)],
        language='c++',
        extra_compile_args=['/std:c++17' if os.name == 'nt' else '-std=c++17'],
    ),
    Extension(
        'polsartools.cprvicpp',
        ['cpp/src/cprvicpp.cpp'],
        include_dirs=[pybind11.get_include(), pybind11.get_include(user=True)],
        language='c++',
        extra_compile_args=['/std:c++17' if os.name == 'nt' else '-std=c++17'],
    ),
]

# class CMakeBuild(build_ext):
#     def run(self):
#         # Ensure the build directory exists
#         if not os.path.exists(self.build_temp):
#             os.makedirs(self.build_temp)
        
#         # Run CMake to configure the build
#         subprocess.check_call(['cmake', '..'], cwd=self.build_temp)
        
#         # Run the build command
#         subprocess.check_call(['cmake', '--build', '.'], cwd=self.build_temp)
        
#         super().run()

# class CustomInstall(install):
#     def run(self):
#         # Call the install command
#         install.run(self)
#         # Optionally run additional post-installation commands here

setup(
    name='polsartools',
    version='0.5.1',
    description='A python package for processing Polarimetric Synthetic Aperture Radar (PolSAR) data.',
    author='Narayanarao Bhogapurapu',
    author_email='bnarayanarao@iitb.ac.in',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'gdal',
        'scipy',
        'click',
        'simplekml',
        'tqdm',  
        'matplotlib',
        'pybind11'
        
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
            'polsartools=polsartools.cli:cli',  # For CLI
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    # cmdclass={
    #     'build_ext': CMakeBuild,
    #     'install': CustomInstall,
    # },
    package_dir={'polsartools': 'polsartools'},
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},  # Use the custom build class
    include_package_data=True,
    zip_safe=False,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
