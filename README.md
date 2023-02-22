# stissplice

[![Documentation Status](https://readthedocs.org/projects/stissplice/badge/?version=latest)](https://stissplice.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://app.travis-ci.com/ladsantos/stissplice.svg?branch=main)](https://app.travis-ci.com/ladsantos/stissplice) [![Coverage Status](https://coveralls.io/repos/github/ladsantos/stissplice/badge.svg?branch=main)](https://coveralls.io/github/ladsantos/stissplice?branch=main)

``stissplice`` is the splicer of Echelle Spectra from Hubble Space Telescope. This code splices Echelle spectra obtained with the Space Telescope Imaging Spectrograph (STIS) instrument. It can be adapted to work with spectra obtained with other instruments as well.

This code is currently in development and accepting suggestions of features to implement, as well as contributions.

Installation
------------

You can install `stissplice` using `pip` or by compiling it from source.

### Option 1: Using `pip` (stable version)

Simply run the following command:
```angular2html
pip install stissplice
```

### Option 2: Compile from source (development version)

First, clone the repository and then navigate to it:
```angular2html
git clone https://github.com/ladsantos/stissplice
cd stissplice
```

And then compile it from source:
```angular2html
python setup.py install
```

You can test the installation from source with ``pytest`` (you may need to install ``pytest`` first):
```angular2html
pytest tests
```

Documentation
-------------

See the page [stissplice.readthedocs.io](https://stissplice.readthedocs.io/).