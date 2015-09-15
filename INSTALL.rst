.. _installing-jacquard:

Installing Jaquard
==================
Jacquard has been tested with Python 2.7 and 3.4 on Windows7, OSX, and \*nix.

Prerequisites
-------------
.. note:: Pip installs all required libraries; see [Installing] below.


* natsort (3.5.2)  
* nosetests, testfixtures (3.0.2), and numpy (>=1.7.1) are required for running
  automated tests

Installing
----------
The easiest way to install Jacquard is through PyPI. Get pip if it's 
not available in your system:

``$ pip install jacquard``

You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/Jacquard``

If you don't have root permissions, you can install locally:

``$ pip install --user jacquard``

