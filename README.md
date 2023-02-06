# 2DChemistry

An atom-based simulation engine; designed to help secondary to early-tertiary students understand core chemsitry concepts.

## What this is for

This toy-model is meant to accomplish the following two design goals:

1. **2DChem can be opened in any browser on static and mobile devices, with one-time access to the internet.** This is accomplished by restricting the code to HTML5, client-side Javascript, and a minimum amount of imported JS libraries (at the time of writing, Chart.JS). 
2. **2DChem provides a semi-quantitative modelling of gas phase molecules, limited to what is generally considered high-school content.** This means zero quantum mechanics and that collision theory governs all chemical recations.

I expect to utilise this software as a part of topics on ideal gas law, reaction rates, activation energies, catalysis, etc.

In other words, this software has been written to succeed the pHET Reactions & rates simulation by the University of Colorado: https://phet.colorado.edu/en/simulations/reactions-and-rates. I wanted a faster engine that can tackle 10,000 atoms (on a decent PC) and approximate complex reaction chains such as hydrocarbon combustion.

## Approximations relative to the real world

1. Space is two-dimensional. 
2. Everything is basically a gas. (In other words, it will be inefficient to force this engine to model liquids and solids.)
3. Atoms are hard circles that undergo rigid-body collisions.
4. Molecules are like-wise rigid bodies with fixed bonds. There are no internal vibrational motions.
5. Energies barriers and transfers are one-tenth their real-world values. This vastly speeds up all process so they can be visualsed.
6. No long range interactions such as electrostatics just yet.

## Other Notes

The reason this engine doesn't do liquid-phase is that you ideally want a different MD engine with Lennard-Jones potentials and much smaller simulation speeds. 

WebAssembly and parallelisation are items on my wishlist in order to speed things up a little more.
