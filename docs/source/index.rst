`stissplice` documentation
==========================

``stissplice`` is the splicer of Echelle Spectra from Hubble Space Telescope. This code splices Echelle spectra obtained with the Space Telescope Imaging Spectrograph (STIS) instrument. It can be adapted to work with spectra obtained with other instruments as well.

This code is currently in development and accepting suggestions of features to implement, as well as contributions.

Installation
------------

Option 1: Using ``pip`` (stable version)
........................................

Simply run the following command:

.. code-block:: bash

   pip install stissplice

Option 2: Compile from source (development version)
...................................................

First, clone the repository and then navigate to it:

.. code-block:: bash

    git clone https://github.com/ladsantos/stissplice
    cd stissplice

And then compile it from source:

.. code-block:: bash

    python setup.py install

You can test the installation from source with ``pytest`` (you may need to install ``pytest`` first):

.. code-block:: bash

    pytest tests


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   usage_example
   how_it_works