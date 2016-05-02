---
title: "The penKdiaph program code - README"
date: 2016-05-02
---

# The penKdiaph program code

_penKdiaph_ is a fortran-based code that runs a Monte Carlo simulation of photons and electron transport and that outputs, among several other things, a series of factors to correct for the effects of a diaphragm on the primary determination of air kerma using a _free-air chamber_. This type of corrections is summarized here:

Burns D T and Kessler C 2009 Diaphragm correction factors for free-air chamber standards for air kerma in x-rays *Phys Med Biol* 54 2737–45

This code is largely based on the _penEasy_ steering programme developed by the Universitat Politecnica de Catalunya [see LICENSE.md file](LICENSE.md). penEasy itself is conceived for opreation in conjunction with the PENELOPE code.

The current version of penKdiaph is based on the 2014 version of penEasy.

Sempau J, Badal A and Brualla L 2011 A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries—application to far-from-axis fields *Med Phys* 38 5887, [DOI address](http://dx.doi.org/10.1118/1.3643029)

Salvat F, Fernandez-Varea J M and Sempau J 2011 PENELOPE-2011: A Code System for Monte Carlo Simulation of Electron and Photon Transport PENELOPE 2011 Workshop, Barcelona, Spain, 4-7 July 2011 pp 1–385

Massimo Pinto, ENEA-INMRI 2014-2016
