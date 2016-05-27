---
title: "The penKdiaph program code - README"
date: 2016-05-02
---

# The penKdiaph program code

_penKdiaph_ is a Fortran-based code that runs a Monte Carlo simulation of photons and electron transport and that outputs, among several other things, a series of factors to correct for the effects of a diaphragm on the primary determination of air kerma using a _free-air chamber_.  This type of corrections is summarized in this research article:

Burns D T and Kessler C 2009 Diaphragm correction factors for free-air chamber standards for air kerma in x-rays *Phys Med Biol* 54 2737–45

The  _penKdiaph_ code is largely based on the _penEasy_ steering programme developed by the Universitat Politecnica de Catalunya [see LICENSE.md file](LICENSE.md). _penEasy_ itself is conceived for opreation in conjunction with the PENELOPE code. Should you become interested in adopting _penKdiaph_ (there are only a few people in the world who work on these things, so I'd be very surprised!) please bear in mind that this repository has no intention to maintain and develop _penEasy_ in any way: _penEasy_ is developed by Josep Sempau from the Universitat Politecnica de la Catalunya and you should contact him for up-to-date codes.

The current version of _penKdiaph_ is based on the 2014 version of _penEasy_. Please note how you can (and you are encouraged to do so) cite the work behind the two projects:

Sempau J, Badal A and Brualla L 2011 A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries—application to far-from-axis fields *Med Phys* 38 5887, [DOI address](http://dx.doi.org/10.1118/1.3643029)

Salvat Francesc, PENELOPE-2014: A Code System for Monte Carlo Simulation of Electron and Photon Transport, PENELOPE 2014 Workshop, Barcelona, Spain, 29 June - 3 July 2014 pp 1–408

Massimo Pinto, ENEA-INMRI 2014-2016.
