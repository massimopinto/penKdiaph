                      penEasy README FILE
                      version 2015-05-30
                       by Josep Sempau
               Universitat Politecnica de Catalunya
                  e-mail: josep.sempau@upc.es

                 Compatible with PENELOPE'2014


CONTENTS:

 0. Copyright and disclaimer
 1. Purpose of this software
 2. Documentation
 3. What you get
 4. How to install it
 5. How to use it
 6. Example
 7. Voxelized geometries
 8. Variance reduction
 9. Parallel execution
10. Restarting a simulation
11. Phase-space files in IAEA format
12. Known issues and limitations
13. Where can I find the latest update?
14. Comments and questions
15. Contributions
16. Preferred citation



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
0. Copyright and disclaimer

penEasy
Copyright Universitat Politecnica de Catalunya

Permission to use, copy, modify and re-distribute copies of the penEasy software (including all its files) or parts of it and its documentation for any purpose is hereby granted without fee, provided that this copyright notice appears in all copies. The Universitat Politecnica de Catalunya makes no representations about the suitability of this software for any purpose. It is provided "as is" without express or implied warranty.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
1. Purpose of this software

PENELOPE [2,3,1] is a Monte Carlo simulation package that describes the transport of photons, electrons and positrons in complex geometries and arbitrary materials. It is freely distributed by the OECD Nuclear Energy Agency Data Bank (http://www.nea.fr). The core of the system is a set of Fortran subroutines that deal with the intricacies of the transport process. To be operative, the system must be completed with a steering main program which, among other things, must define the initial particle states (i.e., the radiation source) and the tallies for the quantities of interest, e.g., the absorbed dose in a certain spatial region. Some examples of main programs are supplied with the distribution package.

Detailed information about PENELOPE can be found in its accompanying documentation (see [2]). Hereafter it is assumed that the reader is familiar with the basic concepts of Monte Carlo simulation of radiation transport and with the operation of PENELOPE.

penEasy is a general-purpose main program for PENELOPE. It provides users with a set of source models, tallies and variance reduction techniques that are invoked from a structured code. Users need only to input all the required information through a simple configuration file and through the usual PENELOPE data files (geometry and materials).

Unfortunately, it is impossible to devise a main program flexible enough to cope with all the imaginable situations. Therefore, users frequently find themselves in the need to adapt an existing code to meet their requirements. In this context penEasy provides a modular code that facilitates the adaptation of the existing routines and the creation of new ones, thus significantly reducing the programming effort.

penEasy is mostly written in Fortran 95, although it has recourse to some features from the Fortran 2003 standard. penEasy, like PENELOPE, is free and open software. In return, we kindly ask users to cite our preferred reference [4] in their publications. Please refer to the complete simulation system as PENELOPE/penEasy.

Thank you for using PENELOPE/penEasy.


References:

[1] J. Bar�, J. Sempau, J. M. Fern�ndez-Varea and F. Salvat, PENELOPE: An algorithm for Monte Carlo simulation of the penetration and energy loss of electrons and positrons in matter, Nucl. Instr. and Meth. B 100 (1995) 31-46.

[2] F. Salvat, J.M. Fern�ndez-Varea and J. Sempau, PENELOPE-2011: A Code System for Monte Carlo Simulation of Electron and Photon Transport, OECD-NEA, Issy-les-Moulineaux, France. Available from http://www.nea.fr.

[3] J. Sempau, E. Acosta, J. Bar�, J. M. Fern�ndez-Varea and F. Salvat, An algorithm for Monte Carlo simulation of coupled electron-photon transport, Nucl. Instr. and Meth. B 132 (1997) 377-390.

*** Preferred citation for penEasy:

[4] J. Sempau, A. Badal and L. Brualla, A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries--application to far-from-axis fields, Med. Phys. 38 (2011) 5887-5895. Available from http://dx.doi.org/10.1118/1.3643029



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
2. Documentation

penEasy documentation consists of this README file, a set of help files located under the ~/documentation directory in the penEasy distribution and instructions interspersed in the data sections contained in the config data file. Please read these instructions carefully before running the code.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
3. What you get

The penEasy package consists of an all-purpose main program, called penEasy.F, and a set of subroutine libraries. After decompressing the distribution ZIP file, your working directory (represented below as ~/) should contain the following files and subdirectories (in alphabetical order):

- ~/documentation/changeHistory.txt
  List of changes with respect to previous versions.

- ~/documentation/dependencies.txt
  List of PENELOPE subroutines that are invoked from penEasy.

- ~/documentation/sectionConfig.txt
  Instructions for the simulation config section of the configuration file.

- ~/documentation/sectionPengeom+penvox.txt
  Instructions for the geometry section of the config file.

- ~/documentation/sectionVarianceReduction.txt
  Instructions for the sections of the config file devoted to variance-reduction techniques.

- ~/documentation/sourceXXX.txt
  ~/documentation/tallyXXX.txt
  Text files containing a description of the input data for each source and tally.

- ~/fortranCode/penaux.f
  Auxiliary routines that initialize PENELOPE/penEasy.

- ~/fortranCode/penEasy.F
  The main program source code.

- ~/fortranCode/penvox.f
  Geometry package for the simulation of voxelized geometries.

- ~/fortranCode/penvr.f
  Implementation of various simple variance reduction techniques.

- ~/fortranCode/sourceXXX.f
  Each source model is coded in a different file, e.g. sourcePhaseSpaceFile.f. They are described in the accompanying ~/documentation/sourceXXX.txt files.

- ~/fortranCode/tallyXXX.f
  Each tally is coded in a different file, e.g. tallySpatialDoseDistrib.f. They are described in the accompanying ~/documentation/tallyXXX.txt files.

- ~/fortranCode/timing.f
  Time routines compatible with the Fortran 95 and subsequent standards.

- ~/gnuplotScripts/*.gpl
  Gnuplot (version 4.2 or greater) scripts to represent graphically the simulation results. Each script has a name associated to a tally. See the example provided below to learn how to execute them and refer to the gnuplot manual for details on how to modify these scripts if needed.

- ~/IAEAcode/*
  C++ code and auxiliary files intended for reading and writing phase-space files (PSFs) in the binary format defined by the International Atomic Energy Agency (IAEA) in ref [1]. See section 11 for details on PSFs and on how to compile these files if required. The latest update of these files are distributed by the IAEA from http://www-nds.iaea.org/phsp.

- ~/run/command.in
  Users may change some simulation settings by editing this file while the program is running--e.g. to stop the execution. See below for more details.

- ~/run/penEasy.exe
  The executable code for Windows systems obtained with the Intel Visual Fortran 64 Compiler XE v.12.0.2.154 with the compiler options -fpp -O2 -DIAEAPSF (see section 11).

- ~/run/penEasy.in
  Sample config file for penEasy. It also serves to run the example case described below.

- ~/run/phantom.geo
  A simple quadrics (PENGEOM model) geometry file that is used in the example case described below. It also contains, for your convenience, the general layout of a quadric geometry file (taken from the PENELOPE distribution).

- ~/run/water.mat
  This is the PENELOPE material data file for water. It is included for your convenience, so that the example in this README (see below) can be run without further preparations.

- ~/voxSample/voxels.vox
  A simple voxels geometry file that complies with the penVox syntax. It contains a detailed description of the data format.

- ~/voxSample/visualizeVoxelsDensity.gpl
  Gnuplot script to visualize the map of material densities in a voxels file.

- ~/voxSample/visualizeVoxelsMaterial.gpl
  Gnuplot script to visualize the map of material indices in a voxels file.

- ~/README.txt
  This file.

All files in the package (except the Windows executable and the reference paper in the documentation directory) are in plain text (i.e. ASCII) format. Except for this readme, all these files use the Unix new line convention, which differs from that of Windows. As a result, some Windows text editors, notably Microsoft Notepad, may not interpret correctly the new line mark and display scrambled text. This problem is solved by using an editor that can read both Unix and Windows formats, for instance the MS-DOS Editor available in most Windows systems (execute "edit"). Usually, compilers accept both new line conventions quietly.


References:

[1] R. Capote, R. Jeraj, C.-M. Ma, D. W. O. Rogers, F. Sanchez-Doblado, J. Sempau, J. Seuntjens and J. V. Siebers, Phase-Space Database for External Beam Radiotherapy, International Atomic Energy Agency Nuclear Data Section Technical Report INDC(NDS)-0484 (2006). Available from http://www-nds.iaea.org/publications/indc/indc-nds-0484.pdf.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
4. How to install it

NOTE: Hereafter is assumed that the installation is carried out on a GNU/Linux system. For concreteness, the gfortran Fortran compiler (freely available from http://gcc.gnu.org) is also assumed. For Windows systems, an executable binary file is included in the distribution (see above) and, therefore, steps b) and c) are usually unnecessary.

The installation is accomplished following these steps:

a) Copy all the penEasy files into a directory on your disc. We shall assume it has been named ~/penEasy/.

b) Copy the files penelope.f, pengeom.f, penvared.f and rita.f from the PENELOPE package into ~/penEasy/fortranCode/.

c) Compile and link penEasy.F :
   $ gfortran penEasy.F -o penEasy.x -O

This generates an executable file (penEasy.x), which should be moved to ~/penEasy/run/. The option '-O' serves to produce and optimized executable. This completes the installation.

Note that, although you must not include all the Fortran source files in the compilation command, those files should be present in the same directory where penEasy is located for the compilation to be successful.

Notice also that the compilation command to be executed differs in case you intend to use phase-space files in the IAEA binary format (see section 11 below). Please refer to that section for further details.

d) In order to be able to visualize the results in graphic format, the penEasy package includes a set of gnuplot scripts. Gnuplot is a function and data plotting software freely available from http://www.gnuplot.info. If you plan to use these scripts, you will also need to have gnuplot installed on your machine.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
5. How to use it

NOTE: In this section it is assumed that the installation has been performed on a GNU/Linux system.

To perform a simulation, take these steps:

a) It is advisable to create a new directory in your working area where all the files for your job will be stored. We shall suppose that this directory is named ~/mySimul/. Always keep the original, unaltered files in ~/penEasy/.

b) Copy the following files to ~/mySimul/.
     ~/penEasy/fortranCode/penEasy.x
     ~/penEasy/run/command.in
     ~/penEasy/run/penEasy.in
     ~/penEasy/gnuplotScripts/*

The geometry and material files required by PENELOPE need also be in ~/mySimul/.

c) Edit the configuration file penEasy.in, where the input data for penEasy is introduced. It is structured in sections, each one starting with a string of the type '[SECTION ...]'. Each source or tally has its own section. Detailed instructions for all sections are given within the file itself.

d) To run penEasy, execute this:
   $ penEasy.x < penEasy.in > penEasy.out

The main program reads from the standard input (the keyboard) and writes to the standard output (the screen), so it must be run redirecting these devices to the appropriate external files with the symbols '<' and '>', as shown in the previous command line.

e) The course of the simulation can be controlled by sending commands to the program while it is running. This is accomplished by modifying the file command.in. The available commands and their coding are briefly explained inside the file itself. For instance, the simulation can be terminated and the final results printed by setting the number of requested histories to zero at any time. The file command.in is scanned periodically.

f) Simulation results are written in a separate file for each tally. These files are named after the corresponding tally and have extension *.dat. The gnuplot scripts included in the distribution use these DAT files to display simulation results graphically. A summary of the execution is also printed to the output file (penEasy.out in the example above).

Note that the DAT files are regularly written at each update interval so that users can track the progress of the simulation.

All statistical uncertainties reported by penEasy are at two standard deviations (2 sigma).


>>>> NOTE: for advanced users on the use of TALLY >>>>

This note is useful in case you plan to modify an existing tally or wish to create your own. Notice that TALLY, which takes the arguments MODE (integer) and ARG (real*8), is called at various points of the penEasy main program where it is expected to perform different actions. The possible values of MODE identify the point at which TALLY is called so that the appropriate action, if any, can be carried out. The following list specifies the situation that corresponds to each mode and the value passed as ARG in each case. In the list, 'E' represents the particle's kinetic energy. By a 'history' we mean the simulation of the primary particle (or particles) produced by a single call to the SOURCE routine and the simulation of all the secondaries generated by it (or them). The required action for the particular case of tallies reporting quantities related with the deposited energy is also described; notice that, in general, the action to be performed depends on the nature of the tally.

Mode -99
Situation:
  A new particle (either primary or secondary) has been retrieved from the stack and its simulation is about to begin.
Required action:
  Its kinetic energy E should be subtracted from deposited energy counters.
Passed argument:
  -E. (This value should be *added* to energy counters.)

Mode -98
Situation:
  The particle has been absorbed.
Required action:
  Its kinetic energy E should be added to deposited energy counters.
Passed argument:
  E.
Comment:
  Note that KNOCK sets E=0 after absorption and, thus, scoring E in energy deposited counters seems irrelevant. However, in penEasy a particle can also be absorbed right after it enters a medium with an E smaller than Eabs, even if no interaction occurs in that medium---i.e., not as a result of a call to KNOCK.

Mode -97
Situation:
  A positron has been absorbed with E>0, i.e., not as a result of a KNOCK but because it entered a region with lower Eabs. Two photons, with an energy of mc^2=511 keV (approx) each, have been pushed to the stack.
Required action:
  2mc^2 should be added to deposited energy counters.
Passed argument:
  2mc^2

Mode -96
Situation:
  One or more particles have been added to the stack as a result of a particle splitting event.
Required action:
  E*(NumParticlesAdded) should be added to deposited energy counters. This mode is not used in case the splitting is produced by the PSF source, which uses mode=0 instead (see below).
Passed argument:
  E*(NumParticlesAdded)

Modes from -1 to -8
Situation:
  An interaction, simulated by KNOCK, has occurred. The type of interaction is returned by KNOCK through the variable ICOL, which takes on values between 1 and 8. The calling MODE is set to ICOL with reversed sign. ICOL labels are as follows:
    Electrons (KPAR=1) and positrons (KPAR=3):
      ICOL = 1 artificial soft event (hinge)
           = 2 hard elastic collision
           = 3 hard inelastic collision
           = 4 hard bremsstrahlung emission
           = 5 inner-shell ionization
           = 6 positron annihilation
           = 7 delta interaction
           = 8 'auxiliary' fictitious interactions
    Photons (KPAR=2):
      ICOL = 1 coherent (Rayleigh) scattering
           = 2 incoherent (Compton) scattering
           = 3 photoelectric absorption
           = 4 electron-positron pair production
           = 7 delta interaction
           = 8 'auxiliary' fictitious interactions
Required action:
  The energy DE lost by the particle should be added to deposited energy counters.
Passed argument:
  DE

Mode 0
Situation:
  A source (e.g., BIGS or PSF) has pushed a new particle to the stack.
Required action:
  The particle's kinetic energy E should be added to deposited energy counters.
Passed argument:
  E.
Comment:
  In the case of particle splitting of a PSF source only one call is made, although several particles will be added to the stack. The passed argument is then (E*NumSplitCopies), which accounts for all the energy deposited.

Mode 1
Situation:
  A new history is about to begin.
Required action:
  None specific.
Passed argument:
  N, the current history number.

Mode 2
Situation:
  Obsolete, not used anymore. Previously, used to signal that the history number N had been changed by the PSF source. Read ~/documentation/changeHistory.txt for further details.

Mode 3
Situation:
  The particle has been moved a distance DS by STEP without any interface crossing and it is about to interact.
Required action:
  None specific. Used by track length estimators to obtain the fluence.
Passed argument:
  DS, the distance travelled.

Mode 4
Situation:
  The particle has been moved a distance DSEF by STEP and an interface crossing has occurred at the end of the flight.
Required action:
  None specific. Used by track length estimators to obtain the fluence.
Passed argument:
  DSEF, the distance travelled.

Mode 5
Situation:
  A particle has been killed. This implies that its simulation is halted without any further action. No energy is deposited and, in the case of a positron, no annihilation radiation is created.
Required action:
  None specific.
Passed argument:
  None (a dummy value).

Mode 6
Situation:
  A history has been finished.
Required action:
  Perform bookkeeping procedures to store the average and variance of the quantities being calculated.
Passed argument:
  N, the history number just finished.

Mode 7
Situation:
  A source routine has created a new particle and has moved it (using STEP) up to the object. The particle has not been pushed to the stack yet.
Required action:
  None specific. Used to tally the current of particles that enter a detector. This information can be useful since, once in the stack, particles created by the source can not be distinguished from secondaries created during the simulation. (A PSF source can create particles that are not primaries.)
Passed argument:
  DSEF, the distance travelled from the point of birth up to the object.


Notice that in all MODE less or equal to zero the value of ARG should be added to the deposited energy counters.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
6. Example

NOTE: In this section it is assumed that the installation has been performed on a GNU/Linux system. Trivial changes (e.g., changing slash '/' by backslash '\' in directory trees) are required for Windows or Mac systems.

Consider a 10 MeV electron beam impinging on a semi-infinite (z>0) water phantom. The beam is modelled as emerging from a point source located at the Cartesian coordinates (0,0,-100), i.e., at 100 cm from the water surface. Electrons are emitted with an initial direction limited to a cone with its axis along the vector (0,0,1) and with an angular semiaperture alpha of 2.86241 deg. Since tg(alpha)=0.05, the field on the water surface is a circle with radius equal to 100*tg(alpha)=5 cm.

We wish to calculate the energy deposited per unit depth interval, that is, a depth-dose curve. More precisely, penEasy will report dE/(rho*dz) in eV*cm^2/g, where dE is the energy deposited in a bin, rho is the mass density of water and dz is the bin width along the z direction. The dose distribution in cylindrical coordinates (as a function of z and the radial distance from the beam axis) will also be computed. The depth interval from z=0 up to z=7 cm will be partitioned into 40 bins. For the cylindrical distribution, the radial interval [0,8] cm will be divided into 80 bins. Results will be reported per history, that is, per unit emitted electron.

Electrons, positrons and photons will be followed down to 100, 100 and 10 keV, respectively. The cutoffs for the production of secondary electrons and bremsstrahlung photons will be set to 100 and 10 keV, respectively. The simulation will be stopped after 100,000 histories, or when the average uncertainty in the depth-dose is 1% (2 sigma), whatever comes first. No variance reduction technique will be applied.

Now, follow these steps:

a) Create the directory ~/mySimul/ in your working area.

b) Copy these files to ~/mySimul/ :
     ~/penEasy/run/command.in
     ~/penEasy/run/penEasy.x
         (See installation procedure above; on a Windows system you may use the ready-made penEasy.exe provided with the distribution.)
     ~/penEasy/run/penEasy.in  (config file prepared for this example)
     ~/penEasy/run/phantom.geo (semiinfinite phantom, PENGEOM syntax)
     ~/penEasy/run/water.mat   (PENELOPE material file)
     ~/gnuplotScripts/tallyCylindricalDoseDistrib-rz.gpl
     ~/gnuplotScripts/tallyParticleTrackStructure.gpl
     ~/gnuplotScripts/tallySpatialDoseDistrib-1D.gpl

As in any other PENELOPE simulation, a material data file, water in our case, generated with the program MATERIAL (included in PENELOPE) is required. For convenience, this file is provided with penEasy for this example.

c) You may want to edit penEasy.in and check that the information introduced in this file reflects the definition of our depth-dose problem. Except for the ParticleTrackStructure, tallies other than those producing the requested distributions have been turned OFF to save CPU time.

A detailed description of the simulation parameters EABS,C1,C2,WCC,WCR and DSMAX can be found in the PENELOPE documentation.

d) Run the program:
   $ penEasy.x < penEasy.in > penEasy.out &
(on Unix systems, the final '&' puts the process in background)

When the simulation ends--it takes about one minute--the sought depth-dose data can be found in tallySpatialDoseDistrib-1D.dat. This file is updated every 100 s ('UPDATE INTERVAL' field in penEasy.in) whilst the program is running. Recall that the quoted uncertainties are at the 2 sigma level.

e) In case you installed gnuplot, you can visualize the depth dose curve by executing the corresponding script from the ~/mySimul/ directory:
   $ gnuplot tallySpatialDoseDistrib-1D.gpl
   (On Windows systems, enter 'wgnuplot' to invoke gnuplot.)

To visualize the cylindrical dose distribution use:
   $ gnuplot tallyCylindricalDoseDistrib-rz.gpl

You can also visualize a few particles tracks executing:
   $ gnuplot tallyParticleTrackStructure.gpl

f) You may send commands to penEasy by editing the file command.in. This file is scanned periodically while the program is running.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
7. Voxelized gometries

Three possible geometry models can be employed in penEasy: (i) quadrics; (ii) voxels, that is, homogeneous volume elements shaped as rectangular prisms, which can be obtained, e.g., from a CT scan; and (iii) a mixture of quadrics and voxels. Detailed instructions on how to select one of these models are given in the file penEasy.in.

To use a voxelized geometry users must provide a valid voxels file. The format of the voxels file is described in detail in the voxels.vox file provided under ~/voxSample/. Keep this file for future reference.

A tally that reports the absorbed dose in each voxel is also provided. As with any other tally, a detailed description of its operation is given in the documentation found in ~/documentation/.

Not all tallies are compatible with voxelized geometries. More precisely, fluence spectra are not reliable when reporting data that refers to the voxelized region. See the documentation of these tallies under ~/documentation/ for more information. To compute the dose distribution inside the voxelized region, the tally specifically tailored for voxels is recommended.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
8. Variance reduction

Some ill-conditioned problems take an unreasonably long simulation time to complete. In these cases it may be useful to have recourse to the so-called variance reduction (VR) techniques, intended to reduce the CPU time required to produce a given statistical uncertainty.

Several simple VR techniques are available in penEasy, some of which take advantage of the VR methods already present in the standard PENELOPE package. With penEasy, their application is controlled through the parameters introduced in the configuration file (the IN file). Please refer to this file for detailed instructions on the various VR options and on how to configure them.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
9. Parallel execution

Monte Carlo (MC) simulations can be parallelized relatively easily. Various strategies and tools have been developed to facilitate the implementation of parallel algorithms, e.g., MPI, openMP, PVM and openMOSIX. All these have been successfully used with MC, but they suffer from the drawback of requiring the modification of the sequential code and/or the installation of additional software, which may be inconvenient for some users.

An alternative solution for penEasy that does not have these limitations is provided by the package clonEasy, which can be freely downloaded from http://www.upc.es/inte/downloads/clonEasy.htm . Its principles and usage are described in [1].

In brief: the same MC job is distributed to several CPUs (the clones) but different random seeds are provided for each clone. After all the executions are done, the partial results are collected and averaged appropriately. To distribute the jobs, the secure shell (ssh) protocol is used. This protocol is usually embedded in Unix-like systems, including all GNU/Linux distributions that we are aware of. This means that any accessible Unix-like computer--e.g. connected to the Internet with a permanent IP address and with a valid user account--can be a clone for our parallel computation. This is achieved without the installation of any additional software.

Random seeds for the PENELOPE pseudo-random number generator (named RANECU) that initiate disjoint sequences, intended to prevent correlations between clones, can be easily obtained with the help of the seedsMLCG code. This code can also be freely downloaded from http://www.upc.es/inte/downloads/seedsMLCG.htm.

To feed different clones with different seeds it is necessary to exploit the feature of penEasy that allows the introduction of the initial seed values through an external file. See the instructions in penEasy.in for more details.


Reference:

[1] Andreu Badal and Josep Sempau, A package of Linux scripts for the parallelization of Monte Carlo simulations, Comput. Phys. Commun. 175 (2006) 440-450. Available from http://www.upc.es/inte/downloads.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
10. Restarting a simulation

A simulation may terminate before the quantity of interest has been obtained with the desired statistical accuracy. This circumstance can be the result of a power cut, of an underestimation of the uncetainty that was really required, or for other reasons.

With penEasy users may restart a simulation and resume the calculation at exactly the point at which it stopped. To this end, a 'snapshot' of the system state is taken periodically, as dictated by the value introduced in the 'INTERVAL BETWEEN DUMPS' field of the configuration (IN) file, and at the end of the simulation. The information obtained from the snapshot is stored in a so-called dump file, whose name is also set in the configuration file.

In a subsequent simulation, the dump file created before can be read to extract all the relevant information about the system state when the first simulation stopped. This includes the random seeds, number of histories, simulation time employed, values of the active tallies and of their uncertainties and state of the phase-space file (PSF) source (i.e., the last particle read, see the documentation for details), if this source model was active. Notice that, if the tally named PSF was active, the restarted simulation will append new particles to the previously created file. This is not the case for the file created by the tally called particle track, whose sole purpose is to display a few trajectories on the screen and not the sequential description of the electromagnetic shower.

The creation of the dump file is optional. If no name is provided, the file will not be created, thus possibly saving some CPU time at the cost of preventing a posterior restart. A simulation can be restarted as many times as desired, each time taking up from the point where the last execution ended.

More details on the syntax and precautions that should be taken when restarting a simulation are given in the help file contained under the ~/documentation/ directory.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
11. Phase-space files in IAEA format

penEasy includes both a source and a tally called phase-space file (PSF, see ~/documentation). A PSF is basically a record of the dynamical state of all particles arriving at (tally) or emerging from (source) a certain region in space. By dynamical state it is meant the type, energy, position, direction, etc. of those particles. The penEasy code allows users to store this information in two different formats, namely, in plain text (i.e., ASCII) or in binary format as specified by the International Atomic Energy Agency (IAEA), see [1].

While the plain text format has the advantages of being human-readable and compatible with all platforms, it also involves the drawback of requiring considerably more disc space than the alternative. Furthermore, the use of the binary format permits the use in penEasy of PSFs obtained from the IAEA database for linear accelerators and cobalt therapy units (see [2]), even when these PSFs where computed with MC codes other than PENELOPE/penEasy.

To facilitate the processing of IAEA-compliant PSFs, the Agency provides a package of routines written in the C++ programming language. penEasy invokes these routines if the user requests that PSFs be handled in IAEA format, but to do so a compilation procedure different from that presented in section 4 must be followed. The compilation command is somewhat more involved due, in part, to the fact that it is necessary to combine code in Fortran (penEasy) and C++ (IAEA routines). Note that the penEasy executable in the distribution has already been created following this procedure (see section 3) and, therefore, Windows users need not repeat the steps listed below.

The compilation procedure is described below in detail assuming a GNU/Linux system and the gfortran compiler. It is also assumed that steps a) and b) of section 4 have been successfully completed. Now do the following:

a) Copy the files
    iaea_config.h
    iaea_header.cpp
    iaea_header.h
    iaea_phsp.cpp
    iaea_phsp.h
    iaea_record.cpp
    iaea_record.h
    utilities.cpp
    utilities.h
from ~/IAEAcode/ to ~/fortranCode/. The rest of files included in ~/IAEAcode/ are given for completeness, but they must not be included in the penEasy compilation directory.

b) Compile and link penEasy.F :
   $ gfortran penEasy.F  *.cpp -o penEasy.x -O -lm -lstdc++ -DIAEAPSF

This generates the executable file penEasy.x. The option '-DIAEAPSF' ensures that the IAEA routines are invoked properly.

The gfortran compiler automatically invokes the C++ compiler for the *.cpp files. Other compiler environments may require a separate compilation of the Fortran and C++ codes. For example, with Intel Fortran for Linux the steps are as follows:

a) Copy the files as in item a) above.

b) Compile the IAEA routines:
   $ icl *.cpp -c -O2

c) Compile the penEasy routines and link with the IAEA ones:
   $ ifort penEasy.F  *.o -fpp -o penEasy.x -O2 -lm -lstdc++ -DIAEAPSF

This generates the executable penEasy.x. For the Intel compiler case, some issues deserve further comments:

- Depending on the system, the options '-lm' and '-lstdc++' must be present in both compilation commands. Note that '-DIAEAPSF' is needed only for the Fortran code.

- The Fortran files in the penEasy distribution all have extension '.F'. The capitalization of the 'F' is important because to implement calls to the IAEA routines penEasy uses what is known as preprocessor directives and the big 'F' tells the compiler that a preprocessor must be invoked before the compiler proper. Unfortunately, not all compilers use this same convention. For those that do not (we have found that Intel Fortran for Windows may not), it can be useful to force the entrance of the preprocessor by adding the option '-fpp' (or '/fpp' in Windows systems) to the compilation command. Alternatively, the file extension '.F' can be changed to '.fpp'. In any case, it is always a good idea to consult the compiler's manual before proceeding.

Note that in order to create or read PSFs in IAEA format it is necessary to select this option in the corresponding source or tally sections in the penEasy configuration file--usually named penEasy.in. Users are requested to provide a filename both for producing a PSF (with the PSF tally) or to read one (when using the PSF source). In both cases two files are involved in the process: one is the PSF header, in plain, human-readable ASCII format, and the other is the actual binary PSF. The two files bear filenames that are composed by the filename provided by the user and the terminations '.IAEAheader' and '.IAEAphsp' for the header and the binary PSF, respectively.


References:

[1] R. Capote, R. Jeraj, C.-M. Ma, D. W. O. Rogers, F. Sanchez-Doblado, J. Sempau, J. Seuntjens and J. V. Siebers, Phase-Space Database for External Beam Radiotherapy, International Atomic Energy Agency Nuclear Data Section Technical Report INDC(NDS)-0484 (2006). Available from http://www-nds.iaea.org/publications/indc/indc-nds-0484.pdf.

[2] http://www-nds.iaea.org/phsp.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
12. Known issues and limitations

- All tallies, with the exception of those reporting pulse height spectra, support variance reduction, i.e., they take the statistical weight into account when scoring quantities of interest. For pulse height spectra, however, some variance reduction techniques may bias the calculated distribution. This is not a limitation of penEasy or PENELOPE, but an unavoidable consequence of the nature of the tally.

- Not all tallies are compatible with voxelized geometries. See the section on voxelized geometries above.

- The use of PENELOPE with external static electromagnetic fields (penfield.f in the PENELOPE package) is not supported. See the PENELOPE manual for a description of the changes that should be introduced to adapt the code.

- Both the Phase-Space File (PSF) tally and the PSF source do not support photon polarization.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
13. Where can I find the latest update?

penEasy can be freely downloaded from http://www.upc.es/inte/downloads/penEasy.htm

It is strongly recommended that you get the latest versions of PENELOPE and penEasy before attempting a new simulation project.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
14. Comments and questions

Comments are most welcome, especially those regarding potential bugs of the code. You can send your comments or questions to josep.sempau[add_an_at_sign_here]upc.es. Please read the documentation carefully before sending comments. Since our human resources are extremely limited, failure to do so may result in not getting a reply to your email. Also, make sure that you have the latest version. Only the latest version is up-to-date. Old penEasy releases, intended for obsolete PENELOPE versions, are not maintained.



>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
15. Contributions

The following colleagues have made important contributions to the development of the code:

- Andreu Badal (e-mail: Andreu.Badal-Soler@fda.hhs.gov).
Dr. Badal contributed to the development of penVox, which handles voxelized geometries in penEasy.

- Immaculada Martinez-Rovira (e-mail: immaculada.martinez@upc.es).
Dr. Martinez-Rovira contributed to the implementation in penEasy of the IAEA format for PSFs.

Some other people have tested the code, given advise or pointed out bugs and potential improvements. See ~\documentation\changeHistory.txt for details.

Many thanks to all of them.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
16. Preferred citation

If penEasy is used for research conducting to publications the following bibliographical reference (given in BibTeX format for your convenience) should be cited:

@ARTICLE{sempau2011,
  AUTHOR  ={J. Sempau and A. Badal and L. Brualla},
  TITLE   ={{A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries--application to far-from-axis fields}},
  JOURNAL ={Med. Phys.},
  VOLUME  ={38},
  YEAR    ={2011},
  PAGES   ={5887 -- 5895},
  NOTE    ={Available at http://dx.doi.org/10.1118/1.3643029}
}


>>>> END OF FILE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
