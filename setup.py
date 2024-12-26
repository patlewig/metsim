from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'metsim'
LONG_DESCRIPTION = 'Python package for generating and/or processing metabolism predictions'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="metsim", 
        version=VERSION,
        author="Grace Patlewicz",
        author_email="<patlewicz.grace@epa.gov>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'

        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Developers",
            "Programming Language :: Python :: 3",
            "Operating System :: Ubuntu :: Linux",
            
        ]
)

