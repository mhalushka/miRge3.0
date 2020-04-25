from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name='testmirge',
        version='3.0',
        author='Arun Patil and Marc Halushka',
        author_email='mhalush1@jhmi.edu',
        url='https://github.com/arunhpatil/Project_120919',
        description='Comprehensive analysis of small RNA sequencing Data', 
        long_description=long_description,
        keywords=['miRge', 'small RNA analysis', 'NGS', 'bioinformatics tools'],  # arbitrary keywords
        license='MIT',
        package_dir={'': 'src'},
        packages=find_packages('src', exclude=['.txt']),
        package_data = {'testmirge3':['models/*.pkl', 'models/*.txt', 'rScripts/*.R']},
        install_requires=['cutadapt>=2.7', 'biopython>=1.76', 'dnaio >= 0.4.1', 'numpy>=1.18.2',
            'scipy>=1.4.1', 'matplotlib>=3.2.1', 'pandas>=0.25.3','scikit-learn>=0.22.2',
            'reportlab>=3.3.0', 'sklearn>=0.0', 'joblib >= 0.14.1',
            ],
        entry_points={'console_scripts': ['miRge3.0 = src.mirge']},
        classifiers=[
            "Development Status :: 1 - Alpha",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Bioinformatics"
            ],
        python_requires='>=3.8',
)
