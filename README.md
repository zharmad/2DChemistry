# 2DChemistry

An atom-based simulation engine; designed to help secondary to early-tertiary students understand core chemistry concepts and visualise common gas-phase natural phenomena & industrial applications. A web demo is available at my GitHub Page: https://zharmad.github.io/2DChemistry/

Based on my preliminary testing (latest May 2024), the performance on Firefox is generally better than Chrome.

## What this is for

This toy-model is meant to accomplish the following design goals:

1. **2DChem can be opened in any browser on static and mobile devices, with one-time access to the internet.** This is accomplished by restricting the code to HTML5, client-side Javascript, and a minimum amount of imported JS libraries (at the time of writing, Chart.JS). 
2. **2DChem provides a semi-quantitative modelling of gas phase molecules, limited to what is generally considered high-school content.** This means that the simulation uses only classical mechanics and collision theory, without going into quantum mechanics details such as rotational-vibrational coupling. It also generally minimised research-standard methods such as Lennard-Jones, periodic boundary conditions, Monte-Carlo sampling, velocity-verlet integration, etc...
3. **2DChem does not hide away the fact that almost all "reactions" written at school represent a set of reaction pathways that involve intermediate products.** This detail is necessary to understand many extension concepts, *e.g.*, why reactions produce by-products and what particle/kinetic/collision theories each imply about observable phenomena.

In other words, this software has been written to succeed the pHET Reactions & rates simulation by the University of Colorado: https://phet.colorado.edu/en/simulations/reactions-and-rates. I wanted a faster engine that can tackle 10,000 atoms (on a decent PC) and approximate complex reaction chains such as hydrocarbon combustion, and I would strongly prefer that the toy model uses reality as the starting point with common examples. The expectation is for us to utilise this software as a part of topics on ideal gas law, reaction rates, activation energies, catalysis, etc.

## Example range of learning outcomes

### Abstract
- The **temperature** of a body of gas is proportional to the **kinetic energy** of its individual molecules. (Not directly velocity.) Thus, heavier gas consistently travels slower on average than lighter gas.
- **Le Chateliers principle:** changes in temperature will bias the formation of reactants and products in an equilibrium reaction. (NB: Pressure-changes is rather difficult to catch both in reality and in this toy model.)
- **Ideal gas law** Pressure corresponds to the frequency of molecules collisions with a given boundary, and is thus proportional to both gas density and temperature.
- The **reaction rate** is influenced by molecular geometry; *e.g.*, an atom to be transferred in a reaction must be able to physically access its expected destination.
- **Radicals** are short-lived molecules with unpaired electrons. These are often generated within a chain of chemical reactions as intermediate species.
- (Wishlist) **Catalysts** speed up reactions by providing new reaction pathways with lower energy barriers. They do NOT physically decrease the activation energies of existing reaction pathways.
- Decomposition reactions, in which a reactant separates into two products, can be catalysed by a photon with energy near the activation energy. (This is called photocatalysis.)

### Specific
- The **Ozone layer** is created and maintained by the absorption of UV light in the stratosphere by both oxygen and ozone molecules. Ozone layer depletion is enhanced by any molecule that can readily react with either ozone and/or oxygen radical.
- **Iodine gas** is purple because the decomposition reaction from di-iodine into iodine radicals has an activation energy roughly equal to the energy of a green photon. 
- (Wishlist) Nitro catalyses combustion reactions by providing reaction paths with significantly lower energy barriers.


## Approximations relative to the real world

1. Space is two-dimensional.
2. Everything is modelled like a gas with no long-range interactions. (In other words, it will be inefficient to force this engine to model liquids and solids.)
3. Atoms and molecules are hard circles that undergo rigid-body collisions.
4. Molecules are like-wise rigid bodies with fixed bonds. There are no internal vibrational motions.
5. Energies barriers and transfers are one-tenth their real-world values. This vastly speeds up all process so they can be visualised.
6. No long range interactions such as electrostatics just yet.

### Note about units.

The units of measurements in this toy model include Kelvins (K), nano/picometers (nm/pm), pico/femtoseconds (ps/fs), atomic-mass-units (amu), kilo-joules (kJ), moles (mol). Note that the reaction rates are many magnitudes faster because energy barriers are 0.1x their real life counterparts as well as the loss of one physical dimension.

## Current feature set

- A set of example gas compositions and reactions.
- Starting conditions for the above composition.
    - initial parameters: composition, area, density, and temperature.
- Dynamic control of gas parameters:
    - External world temperature for equilibration, and externally controlled area.
- Acceptable performance up to ~10,000 atoms on PC and on mobile.
- Basic recording of data over time:
    - Gas paramters: temperature, pressure, area, density, etc.
    - Composition.
    - Available as a dynamics Chart.JS graph, and a manual WIP export to CSV.
- Displays reaction diagrams for each step of the reaction pathway.

- (Partial) Carbon combustion and NOx decomposition pathways to showcase the underlyin complexity of reaction intermediates.
- (Partial) Customised starting compositions is now mostly done by manually creating the conditions and saving to clipboard.
- (Partial) Control: User interface to adjust what each graph shows.
- (Partial) Graph data can be exported to CSV, but the CSV is ugly.
- (Wishlist) Modify the reaction energies and watch responses. This will probably need a graph theory behind the scenes to ensure consistent energies.
- (Wishlist) NPT and NPE conditions that permit the gas to do work on a box wall or four. (Currently on μVT and μVE conditions).
- (Wishlist) Control: User interface to add and remove molecules so as to simulate industrial steps, e.g. during the conversion of sulfur dioxide pollutant into useful sulfuric acid.
- (Wishlist) Moddability: User interface to create molecules and reactions of their own.
- (Wishlist) Feature: Dedicate one wall for a surface catalyst. This enables Haber process catalysis and many others used in industry.

## How does the simulation work?

The toy model doesn't have all the possible reactions and molecules that can exist in reality. Instead, it is explicitly told the molecule types that are allowed to exist in the current simulation, and the reactions that might occur between them. This list of reactions lets the toy determine the necessary collision geometry and activation energy required for a particular reaction, and how to set the velocities and rotations of reaction products.[^1]

To seed the simulation box, the toy reads the gas composition mixture and places molecules one-by-one, randomly, in the available box. It will give up eventually if the box is too small to fit all the molecules. Each molecule is given a starting velocity and rotation corresponding to 2D Maxwell distributions - there is no magic here, kinetic theory specifies that every dimension of motion gets a random initial value. The Gaussian curve is used to determine this "random value", scaled by the molecule's mass or moment-of-inertia.

At the beginning of each simulation step, the integrator simply adds v×dt to the molecule positions. After this, all potential collisions and reactions are checked 
and then resolved. All collisions are simple, elastic collisions if a reaction does not occur. If a molecule collides with a box wall, it is optionally allowed to reset its velocities as part of heat exchange with the outside world.

When a reaction does, the simulation keeps track by deleting reactants and creating products as needed. In this toy model, energy is conserved by redistributing changes due to heats of reaction (ΔH) between the products. This means momentum is no longer conserved. (Note: quantum mechanics and internal motions are needed to conserve both energy and momentum in a reaction process.)

Every ~10 steps, the average velocity and rotation across all molecules is removed (center-of-mass motion), with its energy redistributed by rescaling individual molecule motions.[^2] Information about the simulation is collected regularly and shown in the side tabs.

[^1]: Although fun, it is not feasible to code in all possible reactions. If I want to include a comprehensive but not exhaustive list of intermediates in a methane-oxygen combustion, I would want around 20 molecules types and 50 equilbirium reactions. Add in nitrogen to make it methane-air and this can rise to 53 molecules types and 325 reactions ( see: https://www.cerfacs.fr/cantera/mechanisms/meth.php ). For the purposes of the accuracy, one would also want to look up the Arrhenius activation energies for every pair.

[^2]: This mitigates a common phenomenon for dense simulations where the system of molecules spontaneously begin rotating and moving like a 
single entity (the flying ice-cube phenomena). For this toy model, this happens because we don't respect the momentum and energy conservation laws during collisions; energy is supposed to be absorbed/released into the internal vibrations of molecules, but in this model it has to go into either the RE or KE.

There are industry-standard ways to redistribute this energy across the whole system such that it doesn't cause the "flying ice-cube" artefact. This isn't within scope for our purposes. With the basic center-of-mass motions removed, what you'll see instead in a hot, dense reacting gas are the next collective motions: either the gas swirls in two circles like a convection current, or (if the simulation box is narrow enough) dense aggregates moves back-and-forth like bouncing waves.

### Current directions and limiations.

The reason this engine doesn't do liquid-phase and solid phase is that you ideally want a different MD engine with Lennard-Jones potentials and much smaller time increments.

Still haven't split the configuration into JSON files, but the hooks have been set.

WebAssembly and parallelisation are items on my wishlist in order to speed things up a bit more. The current bottleneck is testing every pair of molecules to check if they are within collision range. This is particle-mesh Ewald territory for liquid-phase simulations with 10^6 atoms, but I think a simple parallelisation solution should be sufficient.

### Known bugs

1. Heat transfer via wall collisions isn't currently implemented correctly/rigorously. The equilibrium system temperature will be lower because faster moleulces experience more collisions that would reset their kinetic energies. Thus, an artificial scaling factor (1.38) is currently applied to the Maxwell 2D velocity distribution to keep the simulation from being too cool relative to the outside world.



