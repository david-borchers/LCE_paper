Contents
========
This is code to go with the paper "A latent capture history model for digital aerial surveys". 

It contains source code for an R package called "twoplane" and a folder called "usercode". You will need to build the package before using it.

It is a standard R package, and has been tested (most recently) with R version 4.0.2. 

The package contains a small component written in C++, so it requires a C++ compiler toolchain to be installed. 

Building and installing in R
----------------------------

To build and install the package in R using devtools, you may use the following commands:

_Start R and set the working directory to the directory containing this file_

install.packages("devtools")

library(devtools)

build()

install()

_Now the package should be ready to load:_

library(twoplane)

Building and installing in RStudio
----------------------------------

Alternatively, to build and install the package in RStudio:
1. Open the project file: LCE_paper.Rproj
2. Open the Build tab
3. Under "More", select "Clean and Rebuild"

Usercode
--------

The usercode folder contains three files:

1. porpoiseFit.R contains code to load and fit to the semi-synthetic porpoise data used in the paper,

2. sim-longlag.r contains code to do the long-lag simulations, and

3. sim-shortlag.r contains code to do the long-lag simulations.


twoplane package
==============

Estimates animal density and related parameters from a mark-recapture line transect survey using two cameras, on which recaptures are not identified, there is animal movement and animals' availability is governed by a Markov model. 

RStudio
=======

If you are an RStudio user, you might want to use the R project file LCE_paper.Rproj as your RStudio project file.
