/*    
    Comes in a new different types. Raw counts and proportions.   
*/
    
class GasComposition {
    constructor(type) {
        if ( undefined === type ) { type = 'count'; }
        
        this.nComponents = 0;
        this.data = {};
        this.type = type;
        
        // the target total number of molecules to generate
        this.nTotal = undefined ;
        this.reset_nTotal();
    }

    reset_nTotal() {
        switch ( this.type ) {
            case 'count':
                this.nTotal = 0;
                break;
            case 'ratio':
                this.nTotal = 0.0;
                break;
            case 'random':
                this.nTotal = 0;
                break;
            default:
                throw `Unrecognised type ${this.type} within GasComposition constructor!`;
        }        
    }
       
    reset() {
        this.data = {};
        this.reset_nTotal();
    }   
        
    get_component_names()  { return Object.keys(this.data); }
    get_component_values() { return Object.values(this.data); }
    get_components() { return Object.entries(this.data); }
        
    add_component(name, val) {
        if ( undefined === this.data[name] ) {
            this.data[name] = val;
            this.nComponents++ ;
        } else {
            this.data[name] += val;
        }
        this.nTotal += val ;
    }
    
    remove_component(name) {
        if ( undefined === this.data[name] ) { return; }
        this.nTotal -= this.data[name] ;
        this.nComponents-- ;
        delete this.data[name];
    }
        
    //Variant with safety Check 
    set_component(name, val) {
        if ( undefined === this.data[name] ) { throw `Gas composition does not have component named ${name}!`; }
        const temp = this.data[name];
        this.data[name] = val;
        this.nTotal += ( val - temp ) ;
    }
    
    add_components_via_array(arrIDs, arrRatios) {
        const n = arrIDs.length;
        if ( n != arrRatios.length ) {throw "ERROR: The two input arrays to add_components_via_array() do not have the same length!"; }
        for ( let i = 0; i < n; i++ ) {
            this.add_component( arrIDs[i], arrRatios[i] );
        }
    }

    debug() {
        console.log( `= = = Gas composition contents: = = =` );
        Object.entries(this.data).forEach( ([key, val]) => {
            console.log(`${key} : ${val}` );
        });
    }
   
    check_against_library(molLib) {
        const arr = molLib.get_defined_molecules();
        Object.keys(this.data).forEach( name => {
            if ( !(name in arr) ) { throw `Molecules ${name} is not defined in the library!`; }
        });
    }
    //const arrMoles = moleculeLibrary.get_moltype_array();
    
    normalise() {
        if ( this.type == 'ratio' ) {
            this.nTotal = 0.0;
            Object.values(this.data).forEach( val => { this.nTotal += val; });
            if ( this.nTotal > 0.0 ) {
                //for ( key in Object.keys(this.data) ) { this.data[key] /= this.nTotal ; }
                Object.keys(this.data).forEach( key => { this.data[key] /= this.nTotal; });
            } else {
                Object.keys(this.data).forEach( key => { this.data[key] = 1.0 / this.nComponents; });
            }
            this.nTotal = 1.0;
        }
        // Else do nothing. Raw counts do not need normalisation.
    }
    
    // Use 0.1% tolerance since we're not using more than 1000 molecules.
    convert_ratio_to_count( nTotal ) {
        if ( this.type != 'ratio' ) { throw "Can only operate on composition ratios!"; }
        if ( (this.nTotal - 1.0) > 0.001 ) { this.normalise(); }
        this.type = 'count';
        this.nTotal = 0;
        Object.entries(this.data).forEach( ([name, val]) => {
            const c = Math.round( val * this.nTotal );
            this.nTotal += c ;
            this.data[name] = c;
        });
    }
    
    convert_count_to_ratio() {
        if ( this.type != 'count' ) { throw "Can only operate on explicit counts!"; }
        this.type = 'ratio';
        Object.entries(this.data).forEach(([name, val]) => {
            this.data[name] = val / this.nTotal ;
        });
        this.nTotal = 1.0;
    }
    
    export_array( bSort ) {
        const arr = [];
        // TODO
        return arr
    }
}

// Interface functions with overall program. 
// This one is about presets
function get_new_preset_gas_composition( type, args ) {
    if ( undefined === type ) { type = 'custom'; }
    if ( undefined === args ) { args = {}; }
    let gc = undefined;
    switch( type ) {
        // case 'demo':
            // let gc = new GasComposition('count');
            // gc.add_component( "He", 80 );
            // break;
        case 'noble gas':
            gc = new GasComposition('ratio');            
            gc.add_component( "He", 16 );
            gc.add_component( "Ne",  8 );
            gc.add_component( "Ar",  4 );
            gc.add_component( "Kr",  2 );
            gc.add_component( "Xe",  1 );            
            gc.normalise();
            break;
        case 'atmosphere':
            gc = new GasComposition('ratio');    
            gc.add_component( "N₂",  0.78 );
            gc.add_component( "O₂",  0.21 );
            gc.add_component( "Ar",  0.01 );
            gc.add_component( "H₂O", 0.02 );
            //gc.add_component( "CO₂", 0.0005 ); // 500ppm
            gc.normalise();          
            //if ( undefined === args["numTotal"] ) { args["numTotal"] = 200; }
            //gc.convert_ratio_to_count( args["numTotal"] );
            break;
        case 'nitrogen dioxide':
            gc = new GasComposition('ratio');
            gc.add_component( "NO₂",  0.8 );
            gc.add_component( "N₂O₄", 0.2 );
            break;
        case 'combustion - H2 and O2':
            gc = new GasComposition('ratio');
            gc.add_component( "H₂",  0.6 );
            gc.add_component( "O₂", 0.4 );
            gc.add_component( "H₂O", 0.0);
            gc.add_component( "O•", 0.0);
            gc.add_component( "H•", 0.0);
            gc.add_component( "OH•", 0.0);
            break;        
        case 'custom':
            gc = new GasComposition('count');
            if ( undefined === args["components"] ) { args["components"] = {}; }            
            Object.entries(args["components"]).forEach( ([name,v]) => {
                gc.add_component( name, v );
            });
            break;
        default:
            throw `Unrecognised gas composition preset ${type}!`;
    }
    return gc;
}