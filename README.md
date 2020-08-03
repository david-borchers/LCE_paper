Contents
========
This is code to go with the paper "A latent capture history model for digital aerial surveys". 

It contains an R package called "twoplane" and a folder called "usercode". You will need to build the package before using it, and because the likelihood is coded in C++ you will need a C++ compiler.

The usercode folder contains three files:

1. porpoiseFit.R contains code to load and fit to the semi-synthetic porpoise data used in the paper,

2. sim-longlag.r contains code to do the long-lag simulations, and

3. sim-shortlag.r contains code to do the long-lag simulations, and


twoplane package
==============

Estimates animal density and related parameters from a mark-recapture line transect survey using two cameras, on which recaptures are not identified, there is animal movement and animals' availability is governed by a Markov model. 
