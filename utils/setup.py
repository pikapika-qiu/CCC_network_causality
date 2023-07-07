from setuptools import setup

setup(
    name='scRNA_utils',
    version='0.1',
    py_modules=['scRNA_utils'],
    # Add any additional dependencies if required
    install_requires=[
        'numpy',
        'pandas',
        'scanpy',
    ],
)
