from setuptools import setup

setup(
    name='polsartools',
    version='0.1',    
    description='A package to process Synthetic Aperture Radar data',
    url='https://github.com/Narayana-Rao/polsartools',
    author='Narayanarao Bhogapurapu',
    author_email='bnarayanarao@iitb.ac.in',
    license='GPL-3.0 license',
    packages=['polsartools'],
    install_requires=['gdal',
                      'numpy',                     
                      ],

    # classifiers=[
    #     'Development Status :: 1 - alpha',
    #     'Intended Audience :: Science/Research',
    #     'License :: OSI Approved ::  GPL-3.0 license',  
    #     'Operating System :: POSIX :: Linux',        
    #     'Programming Language :: Python :: 3',
    # ],
)