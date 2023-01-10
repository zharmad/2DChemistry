/*
    Molecule class that handles everything about a single entity in the simulation.
    = = = = =
    The basic premise is that all molecules will behave as rigid circlular bodies undergoing elastic 2D collisions, regardless of its interior geometry. Thus each molecule has an aggregate position and velocity vector associated with its center of mass, which simplifies collision detection.
    
    Assumes that a global tableOfElements is available.
    
    Default conversion rate to pixels:
    100 pm -> 10 px -> 1 Angs.
*/

// Declare all multi-atom molecules that are supported by this script. Internal offsets positions for drawing purposes.
// Atoms are stated in their draw order. The internal X/Y-axis are defined according to potential molecular dipoles.
/*
    This class defines the types of molecules that are supported. It handles the properties that are shared across all copies of such a molecule type:
    - n, atoms, offsets, mass, rotI, size, colours, molColours 
    
    The convention for the orientation of the molecule is this:
    - the x-axis defines axis of the major bidirectional decomposition for a given molecule.
    - the heavier atom will go first.
*/
class MoleculeLibrary {
      
    constructor() {        
        //if ( undefined === data ) { data = {}; }
        this.data = {} ;
        this.collisionRadiiFactor = 1.0;
        this.tableOfElements = new ElementList();
        this.tableOfElements.add_all_known_elements();
        this.tableOfElements.rescale_radii( 0.75 );
        //this.compute_all_derivative_properties();        
    }
    
    //Does not update all molecule types.    
    get_defined_molecules() { return Object.keys(this.data); }
    get_entry(name) { return this.data[name]; }
    get_molecule_color(name) { return this.data[name].molColour; }

    reset_library() {
        delete this.data;
        this.data = {};
        this.tableOfElements.reset_table();
        this.tableOfElements.rescale_radii( 0.75 );
    }

    add_entry(name, atoms, offsets) {        
        if ( undefined === offsets ) { throw "Not all arguments are present for add_entry() !"; }
        let n = atoms.length;
        if ( n != offsets.length ) { throw "The length of atoms and offset arrys are not the same!"; }        
        //Overwrite existing entries.
        if ( this.data[name] != undefined ) { delete this.data[name]; }        
        this.data[name] = { n, name, atoms, offsets };
        this.convert_offsets(name);
        this.compute_derivative_properties(name);
    }     
    
    convert_offsets(molName) {
        const d = this.data[molName];
        for ( let i = 0 ; i < d.n ; i++ ) {
            if ( !(d.offsets[i] instanceof Vector2D) ) {
                d.offsets[i] = new Vector2D( d.offsets[i][0], d.offsets[i][1] );                
            }
        }        
    }
    
    //Compute derivative molecular properties that will be used by each molecule.
    compute_derivative_properties(molName) {        
        const d = this.data[molName];        
        let m = 0.0 ; let I = 0.0 ; let size = 0.0 ; let cSum = [0,0,0] ;
        let c = []; let ar = [];
        d.nDegrees = (d.n > 1) ? 3 : 2;        
        for ( let i = 0 ; i < d.n ; i++ ) {               
            const elem = this.tableOfElements.get_by_name( d.atoms[i] );
            const radiusAtom  = elem.radius ;
            const r2 = d.offsets[i].norm2();
            m += elem.mass;
            size = Math.max( radiusAtom + Math.sqrt(r2), size );
            I += elem.mass * r2 ;
            for (let j = 0 ; j < 3 ; j++ ) { cSum[j] += elem.colourVec[j]; }
            c.push( elem.parse_colourVec() );
            ar.push( radiusAtom );
        }
        
        // Make monoatomic molecules have a null rotational kinetic energy.
        if ( 0.0 == I ) { I = null; }
        // Assign the newly computed properties back to the data entry.
        d.mass = m ; d.rotI = I ; d.size = size * this.collisionRadiiFactor ;
        for (let j = 0 ; j < 3 ; j++ ) { cSum[j] = Math.floor(cSum[j]/d.n) ; }
        d.molColour = `rgb(${cSum[0]},${cSum[1]},${cSum[2]})`;
        
        d.colours = c; d.radii = ar;
    }
    
    // compute_all_derivative_properties() {
        // //Process this initial library using the table of Elements information to compute the internal properties
        // for (const key in this.data ) {
            // this.convert_offsets(key);            
            // console.log(`Computing properties of ${key}`);            
            // this.compute_derivative_properties(key);
        // }
    // }

    set_molecule_colour( molName, strColour ) {
        this.data[molName].molColour = strColour ;
    }
    
    // Leave room for monoatomid, diatomic and polyatomic later.
    create_molecule( moltype, args ) {
        if ( undefined === args ) { args = {}; }
        var mol = undefined;
        if ( moltype.n > 1 ) {
            mol = new MoleculePolyatomic( moltype, args );
        } else {
            mol = new MoleculeMonoatomic( moltype, args );            
        }
        if ( undefined === mol.name ) { throw "Moltype argument not recognised! Did you mean create_molecule_by_name?"; }
        return mol;
    }
    
    create_molecule_by_name( molname, args ) {
        const moltype = this.get_entry( molname );
        return this.create_molecule( moltype, args );
    }
    
    add_all_known_molecule_types() {
        this.add_entry( 'He', ['He'], [[0,0]] );
        this.add_entry( 'Ne', ['Ne'], [[0,0]] );
        this.add_entry( 'Ar', ['Ar'], [[0,0]] );
        this.add_entry( 'Kr', ['Kr'], [[0,0]] );
        this.add_entry( 'Xe', ['Xe'], [[0,0]] );        
        this.add_entry( 'H₂', ['H','H'], [[32,0],[-32,0]] );
        this.add_entry( 'N₂', ['N','N'], [[71,0],[-71,0]] );
        this.add_entry( 'O₂', ['O','O'], [[66,0],[-66,0]] );
        this.add_entry( 'Cl₂', ['Cl','Cl'], [[102,0],[-102,0]] );
        this.add_entry( 'H₂O', ['H','H','O'], [[-52,76],[-52,-76],[7,0]] );
        this.set_molecule_colour('H₂O','#AFE4DE');
        this.add_entry( 'CO₂', ['O','O','C'], [[116,0],[-116,0],[0,0]] );

        //Nitrogen dioxide equilibrium. Borrow data from https://techiescientist.com/n2o4-lewis-structure/
        this.add_entry( 'NO₂', ['O','O','N'], [[14,110],[14,-110],[-32,0]] );
        this.set_molecule_colour('NO₂', 'brown');
        this.add_entry( 'N₂O₄', ['O','O','N','O','O','N'], [[134,112],[134,-112],[88,0],[-134,112],[-134,-112],[-88,0]] );
        this.set_molecule_colour('N₂O₄', 'white');

        //Hydrogen-oxygen combustion.
        this.add_entry( 'H•', ['H'], [[0,0]] );
        this.set_molecule_colour('H•', '#fbec5d'); //Maize yellow.
        this.add_entry( 'OH•', ['O','H'], [[0,6],[0,-91]] );
        this.add_entry( 'H₂O₂', ['H','H','O','O'], [[94,82],[-94,-82],[0,74],[0,-74]] );
        
        //Ozone layer equilibria. https://en.wikipedia.org/wiki/Ozone_layer
        this.add_entry( 'O•', ['O'], [[0,0]] );
        this.set_molecule_colour('O•', 'DarkRed');
        this.add_entry( 'O₃', ['O','O','O'], [[22,109],[22,-109],[-45,0]] );
        this.set_molecule_colour('O₃', 'LightCyan');
        this.add_entry( 'NO', ['O','N'], [[54,0],[-61,0]] );        
        // nitrate radical bond length is 124 pm, according to Wayne et al., 1991. DOI: 10.1016/0960-1686(91)90192-A
        this.add_entry( 'NO₃•', ['O','O','O','N'], [[21,107],[21,-107],[-41,0],[0,0]] );
        // Three reactions. O2 + hv -> 2O ; O + O2 <-> O3 ; O + O3 -> 2O2
        //this.add_entry( 'N₂O₅', ['O','O','N','O','N','O','O'], [ ] );
        //console.log( this.data['CO₂']);

        // Hydrogen sulfide oxidation and direct thermal decomposition
        // Important step to eliminate sulfur impurities in fuel and biofuel sources.
        // See, e.g.: 
        // th
        // this.add_entry('H₂S', ['H','H','S']);
    }
}

/*
    This class handles the properties that are unique to each molecule:
    - Position, velocity, rotation, colours, etc.
    
    This class will also handle reactions that create and destroy individual molecules.
*/
// This is the base class. Assume polyatomic with all 3 degrees of freedom.
class Molecule {
    constructor( moltype, args ) {
        if ( undefined === moltype ) { throw "ERROR: new Molecule() is not given an molecule type to create!"; }                
        this.name    = moltype.name;
        this.sync_molecular_properties(moltype);
        
        this.p  = undefined;
        this.v  = undefined;
        this.th = undefined;
        this.om = undefined;
        
        //Used by the simulation to skip this molecule, as it has reacted and is about to be deleted.
        this.bIgnore = false;
    }
    
    // Reference the matching molecular properties from the Library.
    sync_molecular_properties(moltype) {
        if ( undefined === moltype )  { throw "ERROR: sync_molecular_properties is not given an molecule type!"; }
        // References to molecular properties.
        this.mass   = moltype.mass;
        this.rotI   = moltype.rotI;
        this.size   = moltype.size;
        this.colour = moltype.molColour;
        
        // References to atomic properties.
        this.nAtoms      = moltype.n;
        this.nDegrees    = moltype.nDegrees;
        this.atomOffsets = moltype.offsets;
        this.atomColours = moltype.colours;
        this.atomRadii   = moltype.radii;
    }
    
    // Various shorthands for retrieval.
    get_size() { return this.size; }
    get_mass() { return this.mass; }
    get_rotI() { return this.rotI; }
    get_colour() { return this.colour; }
   
    //Report in kJ/mol rather than kg.m^2/s^2 
    get_kinetic_energy() { return 0.5 * this.mass * this.v.norm2(); }
    get_rotational_energy() { return 0.5 * this.rotI * this.om**2.0; }
    get_total_energy() { return this.get_kinetic_energy() + this.get_rotational_energy(); }
    get_energies() { return [ this.get_kinetic_energy(), this.get_rotational_energy() ]; }
    
    //momentum
    get_translational_momentum() { return this.v.scaled_copy( this.mass ) ; }
    get_angular_momentum( pRef ) {
        if ( undefined === pRef ) {
            return this.rotI * this.om;
        } else {
            return this.rotI * this.om + this.mass * Vector2D.cross( this.p.subtract(pRef), this.v ) ;
        }
    }   
    
    debug() {
        const entries = Object.entries(this);
        for (const x in entries) {
            var y = entries[x][1];
            if ( y instanceof Vector2D ) {
                console.log(`${entries[x][0]}: ${y[0]} ${y[1]}`);
            }  else {
                console.log(`${entries[x][0]}: ${y}`);
            }
        }
    }
    
    copy_pos_vel_from( mol ) {
        this.p.set_to( mol.p );
        this.v.set_to( mol.v );
        this.th = mol.th;
        this.om = mol.om;
    }
    
    // Draw each atom according to its offset and rotation. Most expensive part of the whole operation.
    // Suspect the arc and fill are the most expensive in total.
    draw(ctxLoc) {
        switch ( globalVars.moleculeColourScheme ) {
            case 'atom':
                for (let i = 0; i < this.nAtoms; i++) {
                    const off = this.atomOffsets[i].rotate( this.th );
                    const xPos = (this.p.vec[0] + off.vec[0])/globalVars.lengthScale ;
                    const yPos = (this.p.vec[1] + off.vec[1])/globalVars.lengthScale ;
                    const radius = this.atomRadii[i]/globalVars.lengthScale;                    
                    ctxLoc.beginPath();
                    ctxLoc.fillStyle = this.atomColours[i];
                    ctxLoc.arc( xPos, yPos, radius, 0, 2 * Math.PI );
                    ctxLoc.fill();
                    ctxLoc.lineWidth = 1;
                    ctxLoc.strokeStyle = '#221100';
                    ctxLoc.stroke();
                }
                break;
            case 'molecule':
                //Get better performance by using a single path that is closed by the fill command.
                ctxLoc.beginPath();                
                ctxLoc.fillStyle = this.colour;
                ctxLoc.lineWidth = 1;
                ctxLoc.strokeStyle = '#221100';
                for (let i = 0; i < this.nAtoms; i++) {
                    const off = this.atomOffsets[i].rotate( this.th );
                    const xPos = (this.p.vec[0] + off.vec[0])/globalVars.lengthScale ;
                    const yPos = (this.p.vec[1] + off.vec[1])/globalVars.lengthScale ;
                    const radius = this.atomRadii[i]/globalVars.lengthScale;
                    ctxLoc.moveTo( xPos + radius, yPos );
                    ctxLoc.arc( xPos, yPos, radius, 0, 2 * Math.PI );
                }
                ctxLoc.stroke();
                ctxLoc.fill();
                break;
            default:
                throw `Unknown colour scheme in global variables: ${globalVars.moleculeColourScheme} !`;
        }
    }

    //Dynamics functions.
    update_position( dt ) {
        this.p.sincr( dt, this.v );
        this.th = ( this.th + this.om * dt ) % ( 2.0 * Math.PI ) ;
    }
    
    rescale_velocities( s ) {
        this.v.scale( s );
        this.om *= s;
    }
    
    resample_speed( T ) {
        if ( undefined === T ) { throw "Resampling requires an input temperature!";}        
        const vNew = random_speed2D( T, this.mass );
        this.v.scale( globalVars.timeFactor * vNew / this.v.norm() ) ;
    }
    
    sample_velocity( T ) {
        if ( undefined === T ) { throw "Sampling requires an input temperature!";}        
        this.v = random_velocity2D( T, this.mass );
        this.v.scale( globalVars.timeFactor );
    }
    resample_velocity( T ) {
        if ( undefined === T ) { throw "Resampling requires an input temperature!";}        
        const vNew = random_velocity2D( T, this.mass );
        this.v.vec[0] = vNew.vec[0] * globalVars.timeFactor;
        this.v.vec[1] = vNew.vec[1] * globalVars.timeFactor;
    }
    
    set_theta( th ) { this.th = th };
    set_omega( om ) { this.om = om };
    resample_omega( T ) {
        if ( undefined === T ) { throw "Resampling requires an input temperature!";}
        //Cheat to just assume that rotations are exactly the same as another dimension of translational movement.        
        this.om = globalVars.timeFactor * random_speed1D( T, this.rotI ) ;
    }
}

//General class.
class MoleculePolyatomic extends Molecule {
    constructor( moltype, args ) {
        super( moltype, args );
        
        const bSample = ( 'bSample' in args ) ? args.bSample : false;
        if ( bSample ) {
            if ( !( 'T' in args ) ) { throw "Cannot get random properties without a temperature argument 'T'!"; }
            // Translational and rotational properties.        
            if ( 'p' in args ) { this.p = args.p } else { this.p = new Vector2D(0,0); }
            if ( 'v' in args ) { this.v = args.v } else { this.sample_velocity( args.T ); }
            if ( 'th' in args ) { this.th = args.th } else { this.th = Math.random() * 2.0 * Math.PI; }
            if ( 'om' in args ) { this.om = args.om } else { this.resample_omega( args.T ); }
        } else {
            // Translational and rotational properties.        
            if ( 'p' in args ) { this.p = args.p } else { this.p = new Vector2D(0,0); }
            if ( 'v' in args ) { this.v = args.v } else { this.v = new Vector2D(0,0); }        
            if ( 'th' in args ) { this.th = args.th } else { this.th = 0.0; }
            if ( 'om' in args ) { this.om = args.om } else { this.om = 0.0; }
        }
    }
}

// TODO: If we incorporate internal vibrational freedoms. Use this.
// class MoleculeDiatomic extends Molecule {}

// Subtype that ignores commands related to rotational freedoms
class MoleculeMonoatomic extends Molecule {
    constructor( moltype, args ) {
        super( moltype, args );
              
        this.th = null;
        this.om = null;
        
        const bSample = ( 'bSample' in args ) ? args.bSample : false;
        if ( bSample ) {
            if ( !( 'T' in args ) ) { throw "Cannot get random properties without a temperature argument 'T'!"; }
            // Translational and rotational properties.        
            if ( 'p' in args ) { this.p = args.p } else { this.p = new Vector2D(0,0); }
            if ( 'v' in args ) { this.v = args.v } else { this.sample_velocity( args.T ); }
        } else {
            // Translational and rotational properties.        
            if ( 'p' in args ) { this.p = args.p } else { this.p = new Vector2D(0,0); }
            if ( 'v' in args ) { this.v = args.v } else { this.v = new Vector2D(0,0); }        
        }
    }   
    
    get_rotI() { return null; }    
    //Report in kJ/mol rather than kg.m^2/s^2 
    get_rotational_energy() { return 0.0; }
    get_total_energy() { return this.get_kinetic_energy(); }
    get_energies() { return [ this.get_kinetic_energy(), 0.0 ]; }
    
    //momentum
    get_angular_momentum( pRef ) {
        if ( undefined === pRef ) {
            return 0.0 ;
        } else {
            return Vector2D.cross( this.p.subtract(pRef), this.v ) ;
        }
    }
        
    copy_pos_vel_from( mol ) {
        this.p.set_to( mol.p );
        this.v.set_to( mol.v );
    }
    
    //Dynamics functions.
    update_position( dt ) { this.p.sincr( dt, this.v ); }
    
    rescale_velocities( s ) { this.v.scale( s ); }    
    
    set_theta( th ) {}
    set_omega( om ) {}    
    resample_omega( T ) {}    
}
