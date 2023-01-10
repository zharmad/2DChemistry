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

globalVars.initialPreset = 'nitrogen dioxide';
globalVars.bPresetsOverwriteParams = true; //Prevent the initial loading from overwriting HTML overrides.

// Molecular colouring section.
globalVars.moleculeColourScheme = 'molecule'; //Choices: 'atom' or 'molecule'

// Summons a map instance. Use map.get(X) as syntax.
function check_html_overrides() {

    const address = window.location.search  
    // Returns a URLSearchParams object instance
    const parameterList = new URLSearchParams(address) 
    
    let mapUserHTMLVars = new Map();
    parameterList.forEach((value, key) => { mapUserHTMLVars.set(key, value); });
  
    // Overwrite these values from user HTML given tags
    for (const key in globalVars) {
        let s = mapUserHTMLVars.get(key);
        if ( undefined != s ) { globalVars[key] = s; }    
    }
}

//NB: these numbers will need sanitisation with the existing minimum and maximum.


// Preset simulation variables go here.
globalVars.presets = {};

// Noble Gas, i.e. hard spheres 
var temp = globalVars.presets['noble gas'] = {};
temp.lengthScale = 20;
temp.timeDelta = 20;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 124;
temp.numComponentsShow = 5;
temp.componentIDs    = ['He', 'Ne', 'Ar', 'Kr', 'Xe'];
temp.componentRatios = [ 16, 8, 4, 2, 1 ];

temp = globalVars.presets['atmosphere'] = {};
temp.lengthScale = 30;
temp.timeDelta = 20;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 400;
temp.numComponentsShow = 4;
temp.componentIDs    = ['N₂', 'O₂', 'Ar', 'H₂O'];
temp.componentRatios = [ 0.78, 0.21, 0.01, 0.02 ];

temp = globalVars.presets['nitrogen dioxide'] = {};
temp.lengthScale  = 20;
temp.timeDelta    = 40;
temp.worldTemperature = 300;
temp.bDoHeatExchange = true;
temp.numMolecules = 160;
temp.numComponentsShow = 2;
temp.componentIDs    = ['NO₂', 'N₂O₄'];
temp.componentRatios = [ 0.8, 0.2 ];

temp = globalVars.presets['combustion - H2 and O2'] = {};
temp.lengthScale  = 30;
temp.timeDelta    = 20;
temp.worldTemperature = 800;
temp.bDoHeatExchange = false;
temp.numMolecules = 400;
temp.numComponentsShow = 3;
temp.componentIDs    = ['H₂', 'O₂', 'H₂O', 'O•', 'H•', 'OH•'];
temp.componentRatios = [ 0.67, 0.33, 0.0, 0.0, 0.0, 0.0 ];
