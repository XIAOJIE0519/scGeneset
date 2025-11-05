from setuptools import setup, find_packages

setup(
    name='scGeneset',
    version='0.0.1',
    author='Shanjie Luan',
    author_email='Luan20050519@163.com',
    description='A toolkit for single-cell gene set analysis',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/xiaojie0519/scGeneset',
    packages=find_packages(),
    install_requires=[
        'scanpy>=1.9.0',
        'pandas>=1.3.0',
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'matplotlib>=3.4.0',
        'seaborn>=0.11.0',
        'scikit-learn>=0.24.0',
        'statsmodels>=0.13.0',
    ],
    python_requires='>=3.13',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.13',
    ],
)
