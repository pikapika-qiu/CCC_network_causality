from setuptools import setup

setup(
    name='scRNA_utils',
    version='0.1',
    py_modules=['scRNA_utils'],# [''],
    # Add any additional dependencies if required
    install_requires=[
        'numpy',
        'pandas',
        'scanpy',
    ],
)
# command to run setup.py for generating other folders if not already exist: "python setup.py bdist_wheel"
#
# go to dist dir in terminal, type "pip install scRNA_utils-0.1-py3-none-any.whl"
# import scRNA_utils as su
# example: su.findDEG()