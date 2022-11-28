# stissplice

``stissplice`` is the splicer of Echelle Spectra from Hubble Space Telescope. This code splices Echelle spectra obtained with the Space Telescope Imaging Spectrograph (STIS) instrument. It can be adapted to work with spectra obtained with other instruments as well.

This code is currently in development and accepting suggestions of features to implement, as well as contributions.

Installation
------------

Currently, the only way to install ``stissplice`` is to compile from source. First, clone the repository and then navigate to it:
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

To compile the documentation, you will need to install the following dependencies:
* [``sphinx``](https://www.sphinx-doc.org/)
* [``jupyter``](https://jupyter.org/install)
* [``nbsphinx``](https://nbsphinx.readthedocs.io)
* [``numpydoc``](https://numpydoc.readthedocs.io)
* ``pandoc`` (see the installation instructions of ``nbsphinx``)

Then run the following commands:
```angular2html
cd docs
make html
```

The documentation pages will be compiled in a folder called `_build`.