/*
    Physical units within the simulation.
    Time     = fs  = 1e-15 s
    Position = pm  = 1e-12 m
    ...
    Velocity = pm fs^-1 
    Mass     = amu ~- 1   g/mol
    Energy   = 1 kJ/mol = 1 amu pm^2 fs^-2
    
    = = = = 
    For an oxygen molecules in 2D, this translates ~ 0.3 pm/fs ~ 1.5 kJ/mol
    
    = = Notes on performance: = =
    The majority of time is spent on drawing circles for each atom. So there is a practical limit to the total number of atoms that can be included, rather than the complexity of the simulation.
    Thus, a large molecule like N2O4 wilil be more expensive.    
*/

var globalVars = {};

// World temperature. Used to determine initial velocities and heat-exchange after collisons.
globalVars.temperature    =  300;
globalVars.temperatureMin =   10;
globalVars.temperatureMax = 1000;

globalVars.bHeatExchange = true;

// Determines the zoom level between the simulatiuon and the visualisation.
// Defaults to 1 pixel = 10 pm = 0.1 Angs.
globalVars.lengthScale    =  10;
globalVars.lengthScaleMin =  10;
globalVars.lengthScaleMax =  50;

// Determines the length of shadows as the screen refreshes each frame.
globalVars.refreshAlpha = 0.4;

// Converts the default units to pm per fs.
globalVars.timeFactor = 1e-3;
globalVars.timeDelta    =  10.0;
globalVars.timeDeltaMin =  10.0;
globalVars.timeDeltaMax = 100.0;

// Simulation specific parameters.
globalVars.numMolecules    =  200;
globalVars.numMoleculesMin =  100;
globalVars.numMoleculesMax = 1000;
globalVars.statisticsUpdateInterval = 100;

globalVars.worldWidth = undefined;
globalVars.worldHeight = undefined;

globalVars.initialPreset = "nitrogen dioxide";
globalVars.bPresetsOverwriteParams = true; //Prevent the initial loading from overwriting HTML overrides.

// Molecular colouring section.
globalVars.moleculeColourScheme = "molecule"; //Choices: 'atom' or 'molecule'

// Summons a map instance. Use map.get(X) as syntax.
function get_html_variables() {

    const address = window.location.search  
    // Returns a URLSearchParams object instance
    const parameterList = new URLSearchParams(address) 
    
    let mapUserHTMLVars = new Map();
    parameterList.forEach((value, key) => { mapUserHTMLVars.set(key, value); });
    return mapUserHTMLVars;
}

function initial_setup_with_html_vars( mapUserHTMLVars ) {
    
    const p = mapUserHTMLVars.get( "initialPreset" );
    if ( undefined != p ) { globalVars.initialPreset = p };
    overwrite_global_values( globalVars.initialPreset );
    
    // Overwrite these values from user HTML given tags
    for (const key in globalVars) {
        let s = mapUserHTMLVars.get(key);
        if ( undefined != s ) { globalVars[key] = s; }    
    }
}

//NB: these numbers will need sanitisation with the existing minimum and maximum.


// Preset simulation variables go here.
globalVars.presets = {};

/*
    Notes:
        1. All reactants, products and expected intermediates should be defined prior to activation.
        2. Place all species that the user are allowed to add prior to the simulation at the beginning. In other words, short-lived intermediates go after the reactants, products, and other participating molecules.
        3. Component ratios don't need to be listed for intermediates. The undefined entries will just resolve to 0.0.
*/
// Noble Gas, i.e. hard spheres.
var temp = globalVars.presets[ "noble gas" ] = {};
temp.lengthScale = 20;
temp.timeDelta = 20;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 124;
temp.numComponentsShow = 5;
temp.componentIDs    = [ "He", "Ne", "Ar", "Kr", "Xe" ];
temp.componentRatios = [ 16, 8, 4, 2, 1 ];

temp = globalVars.presets[ "atmosphere" ] = {};
temp.lengthScale = 30;
temp.timeDelta = 20;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 400;
temp.numComponentsShow = 4;
temp.componentIDs    = [ "N₂", "O₂", "Ar", "H₂O" ];
temp.componentRatios = [ 0.78, 0.21, 0.01, 0.02 ];

temp = globalVars.presets[ "nitrogen dioxide" ] = {};
temp.lengthScale  = 20;
temp.timeDelta    = 40;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 160;
temp.numComponentsShow = 2;
temp.componentIDs    = [ "NO₂", "N₂O₄" ];
temp.componentRatios = [ 0.8, 0.2 ];

temp = globalVars.presets[ "hydrogen iodide equilibrium" ] = {};
temp.lengthScale  = 20;
temp.timeDelta    = 20;
temp.worldTemperature = 600;
temp.bDoHeatExchange = true;
temp.numMolecules = 300;
temp.numComponentsShow = 3;
temp.componentIDs    = [ "H₂", "I₂", "HI", "H•", "I•" ];
temp.componentRatios = [ 0.5, 0.5, 0.0 ];
temp.componentHidePlot = [ "H•", "I•" ];

temp = globalVars.presets[ "ozone layer formation" ] = {};
temp.lengthScale = 30;
temp.timeDelta = 40;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 400;
temp.numComponentsShow = 3;
temp.componentIDs    = [ "N₂", "O₂", "O₃", "O•" ];
temp.componentRatios = [ 0.78, 0.21, 0.01, 0.0 ];
temp.componentHidePlot = [ "N₂" ];

temp = globalVars.presets[ "combustion - H2 and O2" ] = {};
temp.lengthScale  = 30;
temp.timeDelta    = 20;
temp.worldTemperature = 800;
temp.bDoHeatExchange = true;
temp.numMolecules = 400;
temp.numComponentsShow = 3;
temp.componentIDs    = [ "H₂", "O₂", "H₂O", "O•", "H•", "OH•" ];
temp.componentRatios = [ 0.67, 0.33, 0.0 ];
temp.componentHidePlot = [ "O•", "H•", "OH•" ];

temp = globalVars.presets[ "combustion - H2 and O2 adv." ] = {};
temp.lengthScale  = 30;
temp.timeDelta    = 20;
temp.worldTemperature = 800;
temp.bDoHeatExchange = true;
temp.numMolecules = 400;
temp.numComponentsShow = 5;
temp.componentIDs    = [ "H₂", "O₂", "H₂O", "H₂O₂", "O₃", "O•", "H•", "OH•", "HO₂•" ];
temp.componentRatios = [ 0.67, 0.33, 0.0 ];
temp.componentHidePlot = [ "O•", "H•", "OH•", "HO₂•" ];

temp = globalVars.presets[ "combustion - hydrocarbon" ] = {};
temp.lengthScale  = 30;
temp.timeDelta    = 20;
temp.worldTemperature = 600;
temp.bDoHeatExchange = true;
temp.numMolecules = 500;
temp.numComponentsShow = 4;
temp.componentIDs    = [ "C₂H₆", "CH₄", "O₂", "Ar", "CO₂", "H₂O", "H₂", "O•", "H•", "OH•", "H₂O₂", "HO₂•", "CH₃•", "CO", "CH₂O", "CH₃OH", "C₂H₂", "C₂H₄", "CH₂CO" ];
temp.componentRatios = [ 0.1, 0.2, 0.65, 0.05 ];
temp.componentHidePlot = [ "O•", "H•", "OH•" ];

/*
    The objects listing potential reactions go here so that they can be loaded in a modular manner.
*/
globalVars.presetReactions = {}

/*
    Hydrogen iodide decomposition is one of the three main steps in a method to produce hydrogen and oxygen from water - klnown as the sulfur-iodine cycle.
    One fancy thing about this equilibrium is that the dissociation energy of I2 -> I+I is equivalent to green light at 578 nm.
    
    0. Get DeltaH from ANL database as usual: https://atct.anl.gov/Thermochemical%20Data/version%201.118/
        - H: 218, I: 107, I2: 62, HI: 26 
    1. We take the activation energies from this publication on the kinetics: Zhang et al. (2008), DOI: 10.1016/j.ijhydene.2007.10.025
    2. Some additional useful information of catalysed version can be found in, e.g.: Favuzza et al. (2011), DOI: 10.1016/j.apcatb.2011.03.032  
*/
globalVars.presetReactions[ "hydrogen iodide equilibrium" ] = [
    {
        reactantNames: [ "H•", "H•" ], productNames: [ "H₂" ],
        EActivation: 0.0, DeltaH: -43.6, lifetimeActivated: 1000,
    },
    {
        reactantNames: [ "I•", "I•" ], productNames: [ "I₂" ],
        EActivation: 0.0, DeltaH: -15.2, lifetimeActivated: 1000,
    },
    {
        reactantNames: [ "H•", "I•" ], productNames: [ "HI" ],
        EActivation: 0.0, DeltaH: -29.9, lifetimeActivated: 1000,
    },    
    {
        reactantNames: [ "H₂", "I•" ], productNames: [ "HI", "H•" ],
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 360, 360 ],
        productAngles:       [   0,   0 ],
        productAngleRanges:  [ 240, 360 ],                
        EActivation: 16.4, DeltaH: 13.7,
    },
    {
        reactantNames: [ "I₂", "H•" ], productNames: [ "I•", "HI" ],
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 360, 360 ],
        productAngles:       [   0, 180 ],
        productAngleRanges:  [ 360, 240 ],
        EActivation:  1.8, DeltaH: -14.7,
    }
    //The double transfer recation HI + HI <-> H2 + I2 is ignored in this system, as we don't have the corresponding reaction coded and it's also a relatively high barrier.
]


/*
    Ozone layer simple model
    General pathways: https://en.wikipedia.org/wiki/Ozone%E2%80%93oxygen_cycle
    Get DeltaH from ANL database: https://atct.anl.gov/Thermochemical%20Data/version%201.118/
    Model photochemical excitation separately by adding to the internal energy of the molecule (pure rotation here).
    NB: 200nm light gives 598 kJ/mol of energy -> convert to 59.8 in this model.
    Give oxygen 200nm light and ozone 260nm light. See: Section 2 of http://www.ccpo.odu.edu/SEES/ozone/class/Chap_5/index.htm   
    
    Activation energy source for ozone contributions: Sun et al. (2019), DOI: 10.1016/j.pecs.2019.02.002
    TODO: add things that destory the ozone layer.    
*/
globalVars.presetReactions[ "ozone layer formation" ] = [
    {
        reactantNames: [ "O•", "O•" ], productNames: [ "O₂" ],
        EActivation:   0, DeltaH: -49.8, lifetimeActivated: 1000,
    },
    {
        reactantNames: [ "O₂", "O•" ], productNames: [ "O₃" ],
        EActivation: 0, DeltaH: -10.75, lifetimeActivated: 1000,
    },
    {
        reactantNames: [ "O₃", "O•" ], productNames: [ "O₂", "O₂" ],
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 240, 360 ],
        productAngles:       [   0,   0 ],
        productAngleRanges:  [ 360, 360 ],                
        EActivation: 2.0, DeltaH: -39.1,
    }
]
/*
    Combustion reaction models. Datasets of activation energies has been taken from the San Diego mechanism hosted by UC San Diego Combustion Research Group. Available at: http://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html            
    Heat of formation dataset at 298 K are drawn from ANL: https://atct.anl.gov/Thermochemical%20Data/version%201.118/                
*/

/*
    Hydrogen and oxygen reaction chains. DeltaH_formation:
        H = 218 ; O = 249 ; H2O = -242 ; OH = 37            
    This process is simplified to eliminate additional paths such as H₂O₂ (peroxide). For fuller descriptions, see e.g.:
    1.  Dougherty and Rabitz (1980), DOI: 10.1063/1.439114            
*/
globalVars.presetReactions[ "combustion - H2 and O2" ] = [
    // Hydrogen direct decomposition and recombination.
    {
        reactantNames: ["H•", "H•"], productNames: ["H₂"],
        EActivation: 0.0, DeltaH: -43.6, lifetimeActivated: 1000,
    },
    // Oxygen direct decomposition and recombination.
    {
        reactantNames: ["O•", "O•"], productNames: ["O₂"],
        EActivation: 0.0, DeltaH: -49.8, lifetimeActivated: 1000,
    },
    // Water direct decomposition and recombination.
    {
        reactantNames: ["OH•", "H•"], productNames: ["H₂O"],
        EActivation: 0.0, DeltaH: -49.7, lifetimeActivated: 1000,
        reactantAngles:      [   0,   0 ], // Filled with 0.0 if not given.
        reactantAngleRanges: [ 240, 360 ], // Filled with 360 if not given. 
        productAngles:       [ 0 ],
        productAngleRanges:  [ 0 ],
    },
    // OH radical direct decomposition and recombination.
    {
        reactantNames: ["O•", "H•"], productNames: ["OH•"],
        EActivation: 0.0, DeltaH: -43.0, lifetimeActivated: 1000,
    },
    // Radical propagation 1: hydrogen and oxygen molecule
    {
        reactantNames: ["O₂", "H•"], productNames: ["OH•", "O•"],
        EActivation: 7.1, DeltaH: 6.8,
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 360, 360 ],
        productAngles:       [   0,   0 ],
        productAngleRanges:  [ 240, 360 ],
    },
    // Radical propagation 2: oxygen and hydrogen molecule
    {
        reactantNames: ["H₂", "O•"], productNames: ["OH•", "H•"],
        EActivation: 2.6, DeltaH: 0.6,
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 360, 360 ],
        productAngles:       [ 180,   0 ],
        productAngleRanges:  [ 240, 360 ],
    },
    //  Collision-based water formation 1. 
    {
        reactantNames: ["OH•", "H₂"], productNames: ["H₂O", "H•"],
        EActivation:  1.5, DeltaH: -6.1,
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 240, 360 ],
        productAngles:       [ 180,   0 ],
        productAngleRanges:  [ 240, 360 ],
    },
    //  Collision-based water formation 2. (Don't worry about collision symmetry just yet. This need a more advaned angle algorithm).
    {
        reactantNames: ["OH•", "OH•"], productNames: ["H₂O", "O•"],
        EActivation:  0.0, DeltaH: -6.7,
        reactantAngles:      [   0,   0 ], 
        reactantAngleRanges: [ 360, 360 ],
        productAngles:       [ 180,   0 ],
        productAngleRanges:  [ 240, 360 ],
    },    
    // Self reaction of hydrogen. Guesstimate as not experimentally measurable.
    {
        reactantNames: ["H₂", "H•"], productNames: ["H•", "H₂"],
        EActivation: 5.0, DeltaH: 0.0,
    }
    // Self reaction of oxygen not used -> goes to ozone.
]

// TODO: This and the carbon combustion requires an angle-based determination of readction mechanisms. Example: OH+OH resolves to HOOH and H2O + O depending on angle.

/*
    Add in peroxide and ozone pathways that are involved in combustion.
    1. http://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html1.  
    2. Sun et al. (2019), DOI: 10.1016/j.pecs.2019.02.002
*/
globalVars.presetReactions[ "combustion - H2 and O2 adv." ] = [
    
]