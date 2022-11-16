Welcome to the `wesh` documentation!
====================================

``wesh`` is the Weaver of Echelle Spectra from Hubble (Space Telescope). This code splices Echelle spectra obtained with the Space Telescope Imaging Spectrograph (STIS) instrument. It can be adapted to work with spectra obtained with other instruments as well.

This code is currently in development and accepting suggestions of features to implement, as well as contributions.

Installation
------------

Currently, the only way to install ``wesh`` is to compile from source. First, clone the repository and then navigate to it:

.. code-block:: bash

    git clone https://grit.stsci.edu/ldsantos/wesh.git
    cd wesh

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