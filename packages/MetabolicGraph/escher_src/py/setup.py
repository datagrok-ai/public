import sys
from os.path import join, dirname, realpath

from setuptools import setup, find_packages

directory = dirname(realpath(__file__))
sys.path.insert(0, join(directory, 'escher'))
version = __import__('version').__version__
full_version = __import__('version').__full_version__
package = __import__('version').package
port = 8789

setup(
    name='Escher',
    version=full_version,
    author=package['author'],
    url=package['homepage'],
    description=package['description'],
    keywords=', '.join(package['keywords']),
    license=package['license'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Visualization',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: OS Independent'
    ],
    packages=find_packages(),
    include_package_data=True,
    data_files=[
        (
            'share/jupyter/nbextensions/escher',
            [
                'escher/static/extension.js',
                'escher/static/escher.min.js',
                'escher/static/escher.min.js.map',
            ]
        ),
        (
            'etc/jupyter/nbconfig/notebook.d',
            ['escher.json'],
        )
    ],
    install_requires=[
        'Jinja2>=2.7.3,<3',
        'pytest>=4.0.1,<5',
        'cobra>=0.5.0',
        'jsonschema>=3.0.1,<4',
        'ipywidgets>=7.4.0,<8',
        'pandas>=0.18'
    ],
    extras_require={
        'docs': [
            'sphinx>=2.1.1,<3',
            'sphinx-rtd-theme>=0.4.3,<0.5',
            'nbsphinx>=0.4.2,<0.5'
        ],
    },
)
