from setuptools import setup, find_packages, find_namespace_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name='mirge',
        version='0.1.1',
        author='Arun Patil and Marc Halushka',
        author_email='mhalush1@jhmi.edu',
        url='https://test.pypi.org/legacy/',
        description='Comprehensive analysis of small RNA sequencing Data', 
        long_description=long_description,
        keywords=['miRge', 'small RNA analysis', 'NGS', 'bioinformatics tools'],  # arbitrary keywords
        license='MIT',
        #packages=['mirge','mirge.classes','mirge.libs', 'mirge.forgi.',],
        #packages=['mirge',],
        package_dir={'mirge': 'mirge'},
        packages=find_packages(),
        #packages=find_packages(where='mirge'),
        #packages=find_namespace_packages(include=['classes','libs']),
        #packages=find_packages('tstmirge1', exclude=['.txt']),
        package_data = {'':['models/*.pkl', 'models/*.txt', 'rScripts/*.R']},
        install_requires=['cutadapt==2.7', 'biopython', 'dnaio', 'numpy',
            'scipy', 'matplotlib', 'pandas','scikit-learn',
            'reportlab',
            ],
        entry_points={'console_scripts': ['miRge3.0 = mirge.__main__:main']},
        classifiers=[
            "Development Status :: 1 - Planning",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
            ],
        include_package_data=True,
        python_requires='>=3.8',
)
