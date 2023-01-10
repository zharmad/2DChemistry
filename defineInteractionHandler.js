/*
    This class handles a defined set of reactions between molecules.
    = = = =
    It is invoked with an array containing the list of molecules to be characterised.

    When a reaction check is considered successful, this handler returns the product molecules and their geometries for the simulation to then directly incorporate. This handler does not get to remove the rectants.
    
    This handler computes the velocities so as to try to preserve total kinetic and rotational energy.

    This requires defineMolecules.js to function. Does not need to know about molecule libraries alead of times, but will need to initialise them in order to specify in advance the relevant moltypes.
*/

class InteractionHandler {
    constructor( args ) {
        if ( undefined === args ) { args = {}; }
        //this.nEntries = 0;
        if ( undefined != args.species ) {
            this.initialise_species( args.species );
        } else {
            this.species = [];            
            this.dict = {};
            this.dictSelf = {};
        }
        //this.initialise( args.species );
        this.moleculeLibrary = args.moleculeLibrary;
        console.log("Gas interaction handler initialised.");
    }    

    set_molecule_library( ml ) { this.moleculeLibrary = ml; }
    
    // Create a bidirectional dictionary so that both orders point to the same object.
    initialise_species(arrNames) {
        console.log(arrNames);
        this.species = arrNames;
        this.dict = {};
        this.dictSelf = {};        
        arrNames.forEach( m => {
            this.dict[m] = {}
            this.dictSelf[m] = undefined;
        });
        arrNames.forEach( m1 => {
            arrNames.forEach( m2 => {
                this.dict[m1][m2] = this.dict[m2][m1] = undefined;
            });
        });
    }
    
    debug() {
        console.log(this);
        this.species.forEach( name => {
            console.log( `= = Defined molecule species: ${name}` );            
            if ( undefined != this.dictSelf[name] ) {
                console.log( `...Self reaction entry for ${name}:` );
                this.dictSelf[name].debug();
            }
            for ( const k in this.dict[name] ) {
                if ( undefined != this.dict[k][name] ) {
                    console.log( `...Collision entry between ${name} and ${k}` );
                    this.dict[k][name].debug();
                }
            }
        });
    }
    
    //Only self reactions have time dependence. Collision reactions are instant!
    update_time_dependence( dt ) {
        for ( const key in this.dictSelf ) {
            if ( undefined != this.dictSelf[key] ) { this.dictSelf[key].update_time_dependence( dt ); }
        }
    }
    
    /*
        The new function to parse an input reaction property. TODO convert to JSON or some other easy data format.
        Required arguments:
            - arrays: reactantNames, reactantAngles, reactantAngleRanges
                      productNames,  productAngles,  productAngleRanges
            - scalars: EActivation, DeltaH,
        Optional arguments:
            - scalars: lifetimeActivated, angleReactionOffset.
            - string: unitAngle
            - boolean: bDoForward, bDoReverse 
    */    
    parse_input_reaction( args ) {        
            
        //Initial settings and fil,ling out defaul definitions.
        var bDoForward = ( undefined === args.bDoForward ) ? true : args.bDoForward ;
        var bDoReverse = ( undefined === args.bDoReverse ) ? true : args.bDoReverse ;        
        if ( false == bDoForward && false == bDoReverse ) { return; } // Useful for bug-fixing? Do not throw.
        
        if ( undefined === args.unitAngle ) { args.unitAngle = 'degrees'; }
        if ( undefined === args.reactantAngles ) { args.reactantAngles = Array( args.reactantNames.length ).fill( 0.0 ); }        
        if ( undefined === args.productAngles ) { args.productAngles = Array( args.productNames.length ).fill( 0.0 ); }
        if ( undefined === args.reactantAngleRanges ) { args.reactantAngleRanges = Array( args.reactantNames.length ).fill( 360.0 ); }        
        if ( undefined === args.productAngleRanges ) { args.productAngleRanges = Array( args.productNames.length ).fill( 360.0 ); }
        
        // Linking to the molecular library for creating products during the simulation.
        args.moleculeLibrary = this.moleculeLibrary ;
        
        this.convert_angle_units( args );

        // Grab the moltype entries themselves here.
        args.reactants = []; const nAtomsReac = [];
        const nReactants = args.reactantNames.length;
        for ( var i=0;i<nReactants;i++ ) {
            args.reactants.push( this.moleculeLibrary.get_entry( args.reactantNames[i] ) );
            nAtomsReac.push( args.reactants[i].n );
        }

        args.products  = []; const nAtomsProd = [];
        const nProducts  = args.productNames.length;
        for ( var i=0;i<nProducts;i++ ) {
            args.products.push( this.moleculeLibrary.get_entry( args.productNames[i] ) );
            nAtomsProd.push( args.products[i].n );
        }
        //Complete determination of the type of reaction class to invoke.        
        switch ( nReactants ) {
            case 2:
                //Combination or general collision.
                switch( nProducts ) {
                    case 2:
                        // Classic two on two collision reaction.
                        if ( bDoForward ) {
                            this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ] = new reactionTransfer( args );
                            this.dict[ args.reactantNames[1] ][ args.reactantNames[0] ] = this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ];
                        }
                        if ( bDoReverse ) {
                            this.swap_reactant_and_products( args );
                            this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ] = new reactionTransfer( args );
                            this.dict[ args.reactantNames[1] ][ args.reactantNames[0] ] = this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ];
                        }
                        break;
                    case 1:
                        //Combination reaction.
                        if ( bDoForward ) {
                            this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ] = new reactionCombination( args );
                            this.dict[ args.reactantNames[1] ][ args.reactantNames[0] ] = this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ];
                        }
                        if ( bDoReverse ) {
                            this.swap_reactant_and_products( args );
                            this.dictSelf[ args.reactantNames[0] ] = new reactionDecomposition( args );
                        }
                        return;
                        break;
                    default:
                        throw `ERROR: Unsupported number of reactants and products! ${nReactants} -> ${nProducts}`;
                }
                break;
            case 1:
                //Decomposition or single-transformation.
                switch( nProducts ) {
                    case 2:
                        // Decomposition reaction.                       
                        if ( bDoForward ) {
                            this.dictSelf[ args.reactantNames[0] ] = new reactionDecomposition( args );                                
                        }
                        if ( bDoReverse ) {
                            this.swap_reactant_and_products( args );
                            this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ] = new reactionCombination( args );
                            this.dict[ args.reactantNames[1] ][ args.reactantNames[0] ] = this.dict[ args.reactantNames[0] ][ args.reactantNames[1] ];
                        }
                        return;
                        break;
                    case 1:
                        // One to one transformation reaction. Can be used to switch between stoichiometric forms.
                        throw "ERROR: one-to-one transformation reactions are not currently supported!";
                        break;
                    default:
                        throw `ERROR: Unsupported number of reactants and products! ${nReactants} -> ${nProducts}`;
                }
                break;
        }
        //console.log("...Added interaction between.");
    }
    
    convert_angle_units( args ) {
        if ( 'degrees' === args.unitAngle ) {
            const f = Math.PI/180.0;
            for(var i=0;i<args.reactantAngles.length; i++) { args.reactantAngles[i] *= f; }
            for(var i=0;i<args.reactantAngleRanges.length; i++) { args.reactantAngleRanges[i] = Math.cos( args.reactantAngleRanges[i] * f/2 ); }
            for(var i=0;i<args.productAngles.length; i++) { args.productAngles[i] *= f; }
            for(var i=0;i<args.productAngleRanges.length; i++) { args.productAngleRanges[i] = Math.cos( args.productAngleRanges[i] * f/2 ); }
            if( undefined != args.angleReactionOffset ) { args.angleReactionOffset *= f; }
            args.unitAngle == 'radians';
        }
    }
    
    swap_reactant_and_products( args ) {       
        let temp = undefined;
        
        temp = args.reactants;
        args.reactants = args.products;        
        args.products  = temp;

        
        temp = args.reactantNames;
        args.reactantNames  = args.productNames;
        args.productNames = temp;
        
        temp = args.reactantAngles;
        args.reactantAngles  = args.productAngles;
        args.productAngles = temp;                
        
        temp = args.reactantAngleRanges;
        args.reactantAngleRanges  = args.productAngleRanges;
        args.productAngleRanges = temp;
        
        args.EActivation -= args.DeltaH;
        args.DeltaH *= -1;
    }        
    
    remove_reaction( m1, m2 ) {
        delete this.dict[m1][m2];
        delete this.dict[m2][m1];
        //this.nEntries--;
    }
        
    // Entry point from Simulation object for all unimolecular processes.
    // Triggered for every molecule.
    // Return either null or an array containing new molecules.
    process_molecule_selfinteraction( mol ) {
        if ( undefined === this.dictSelf[mol.name] ) { return null; }
        //return entry.process_reaction( { mol } );
        return this.dictSelf[mol.name].process_reaction( new molSystem( [mol] ) );
    }

    // Entry point from Simulation object for all bimolecular interactions.
    // Triggered when two molecules comes within a specific distance.
    // Return either null or an array containing new molecules.
    process_molecule_encounter( mol1, mol2 ) {
        //Container for variables to be pased between different functions.
        const obj = new molSystem( [mol1, mol2] );
        //First confirm whether a collision will occur.        
        const bCollide = obj.confirm_molecule_collision();
        if ( false === bCollide ) { return null; }
        
        // Collision has definitely occurred. May want to shift molecules backwards and obtain more accurate contact point. THis is important for some edge cases.
        
        // Check whether there is a possiblity of reaction.
        let entry = this.dict[mol1.name][mol2.name];
        //if ( undefined === entry ) {throw `The names ${mol1.name} & ${mol2.name} are not part of the defined dictionary!`; }

        // Elastic collision exit point for molecules that has no possibility of reaction.        
        if ( undefined === entry ) {
            obj.resolve_collision();
            //console.log(`Collision between ${mol1.name} & ${mol2.name}, no possible reaction.`);
            return null;
        }
        
        // Determine if a reaction might occur due to this simple collision.
        // NB: All angle and energy checks should be specific to the reaction and handled within the object.        
        const ret = entry.process_reaction( obj );
        if ( null != ret ) {
            //console.log(`Collision between ${mol1.name} & ${mol2.name}, reaction success.`);
            return ret;
        } else {
            // Elastic collision exit point for molecules that could react but have failed above checks.
            //console.log(`Collision between ${mol1.name} & ${mol2.name}, reaction failed.`);
            obj.resolve_collision();
            return null;
        }
    }
                
}

// Class for a system of molecules defined as their centers of mass,
// all known to be within interaction distance.
// Useful shorthands for computing centers of mass.
// Meant to be extensible towards resolving multibody collisiosn.
class molSystem {
    constructor( arrMol ) {
        this.arrMol = arrMol;
        this.nMol = arrMol.length;
        this.mass = 0.0;
        this.p = undefined;
        this.v = undefined;
        this.pRel = undefined;
        this.vRel = undefined;
        
        this.pContact    = undefined;
        this.pContactRel = undefined;
        this.vecNormal   = undefined;
        
        //this.reset_component_energies();
    }

    reset_component_energies() {
        // Collision Energies.
        this.TEcom      = 0.0;
        // this.REcom = 0.0; Iw^2 -> sum m r^2 v^2/r^2 -> sum m v^2, which is just translational energy of molecular components as below.
        this.TEradials  = [];
        this.TEtangents = [];
        this.REradials  = [];
        this.REtangents = [];        
    }

    get mol()  { return this.arrMol[0]; }
    get mol1() { return this.arrMol[0]; }
    get mol2() { return this.arrMol[1]; }
    
    debug() {
        console.log(this);
        for (const k in this) {
            console.log(`${k}`, this[k]);
        }
    }
    
    check_integrity() {
        for (const mol of this.arrMol) {
            if ( Number.isNaN( mol.p.x ) || Number.isNaN( mol.p.y ) ) {
                mol.debug();
                throw "NaN position values have been detected in molSystem!";
            }            
            if ( Number.isNaN( mol.v.x ) || Number.isNaN( mol.v.y ) ) {
                mol.debug();
                throw "NaN velocity values have been detected in molSystem!";
            }            
        }
    }
    
    static distribute_angmoment_with_partial_energy_conservation( s1, s2 ) {
        //if( undefined === s1.pRel ) { s1.calculate_COM_configuration(); }
        if( undefined === s2.pRel ) { s2.calculate_COM_configuration(); }
        
        const angMomentOld = s1.get_angular_momentums_from_mols();
        const signAM = Math.sign( angMomentOld );                    
        //console.log( "old angmoment", angMomentOld );
        let f = 0.0;
        s2.arrMol.forEach( (mol, i) => {
            f += s2.pRel[i].norm() * Math.sqrt( 2.0 * mol.mass);
            if ( mol.rotI > 0.0 ) { f += Math.sqrt(mol.rotI); }
        });
        const KEnew = ( angMomentOld / f ) ** 2.0;                
        // The incrementation procedure isn't really energy conserving. Should try to fix...?
        s2.arrMol.forEach( (mol, i) => {            
            mol.v.sincr(  Math.sqrt( 2.0 * KEnew / mol.mass ), Vector2D.scalar_cross( signAM, s2.pRel[i].unit() ) );
            if ( mol.rotI > 0.0 ) { mol.om = signAM * Math.sqrt( KEnew / mol.rotI ); }
        });
        //console.log( "new angmoment", s2.get_angular_momentums_from_mols() );
        //console.log( s1.get_total_energy_from_mols(), s2.get_total_energy_from_mols() );
    }
    
    // Rescales all motions in system 2 so as to match the starting energy of system 1 with a DeltaH.
    static rescale_velocities_by_energy_conservation( s1, s2, DeltaH ) {
        if ( undefined === DeltaH ) { DeltaH = 0.0; }
        const EOld = s1.get_total_energy_from_mols();
        const ENew = s2.get_total_energy_from_mols();
        const EExcess = EOld - DeltaH - ENew;
        if ( ENew + EExcess < 0) { console.log(ENew, EExcess); throw "Negative output energy found! Aborting."; }                    
        let ratio = Math.sqrt( (ENew + EExcess) / ENew );
        s2.arrMol.forEach( (mol, i) => { mol.rescale_velocities( ratio ); });
        //console.log( EOld, s2.get_total_energy_from_mols(), DeltaH );
    }
    
    calculate_COM_configuration() {
        this.p = new Vector2D( 0,0 );
        this.v = new Vector2D( 0,0 );
        this.pRel = [], this.vRel = [];
        this.arrMol.forEach( (mol, i) => {
            const mass = mol.get_mass();
            this.mass += mass;
            this.p.sincr( mass, mol.p );
            this.v.sincr( mass, mol.v );
            this.pRel.push( mol.p.copy() );
            this.vRel.push( mol.v.copy() );                                    
        });        
        this.p.scale( 1.0 / this.mass );// COM Position in world frame
        this.v.scale( 1.0 / this.mass );
        // relative positions in COM frame
        this.pRel.forEach( p => { p.decr( this.p ); });
        this.vRel.forEach( v => { v.decr( this.v ); });
        if ( undefined != this.pContact ) {
            this.pContactRel = this.pContact.subtract( this.p );
        }
    }

    calculate_angular_momentums() {
        //this.AngMomentumsExt = []; this.AngMomentumsInt = [];
        this.AngMomentumsExt = 0.0; this.AngMomentumsInt = 0.0;
        this.arrMol.forEach( (mol, i) => {        
            // this.AngMomentumsExt.push( mol.mass * Vector2D.cross( this.pRel[i], this.vRel[i] ) );
            // this.AngMomentumsInt.push( mol.get_angular_momentum() );
            this.AngMomentumsExt += mol.mass * Vector2D.cross( this.pRel[i], this.vRel[i] );
            this.AngMomentumsInt += mol.get_angular_momentum();
        });
    };

    get_angular_momentums_from_mols() {
        if( undefined === this.p ) { this.calculate_COM_configuration(); }
        let ret = 0.0;
        this.arrMol.forEach( (mol, i) => { ret += mol.get_angular_momentum( this.p ); });
        return ret;
    }

    calculate_component_energies() {        
        this.reset_component_energies();
        // Compute the breakdown of all energies.
        // This allows for a more detailed breakdown of all energies, extracting components that are relevant for a reaction.

        this.TEcom = 0.5 * this.mass * this.v.norm2();        
        this.arrMol.forEach( (mol, i) => {
            const sRadial = Vector2D.dot( this.vRel[i], this.vecNormal );
            this.TEradials.push( 0.5* mol.mass * sRadial*sRadial );            
            const sTangent = Vector2D.cross( this.vRel[i], this.vecNormal );
            this.TEtangents.push( 0.5* mol.mass * sTangent*sTangent );
            const RE = mol.get_rotational_energy();
            if ( RE == 0.0 ) {
                this.REradials.push( 0.0 ); this.REtangents.push( 0.0 );
            } else {
                const frac = Vector2D.cross( this.pContact.subtract( this.pRel[i] ).unit(), this.vecNormal );
                this.REradials.push( frac * RE ); this.REtangents.push( (1.0-frac) * RE );
            }
        });
    }

    get_rotational_energy() {
        let ret = 0.0;
        for ( let i = 0; i < this.nMol; i++ ) { ret += this.REradials[i] + this.REtangents[i]; }
        return ret;
    }

    get_radial_components_sum() {
        let ret = 0.0;
        for ( let i = 0; i < this.nMol; i++ ) { ret += this.TEradials[i] + this.REradials[i]; }
        return ret;
    }

    get_total_energy_from_mols() {
        let ETot = 0.0;
        this.arrMol.forEach( ( mol, i ) => { ETot += mol.get_total_energy(); } );
        return ETot;
    }   

    // This one is used when the COM energies have already been calculated.
    get_total_energy( bIncludeTEcom ) {
        if ( undefined === bIncludeTEcom ) { bIncludeTEcom = true; }
        let ETot = 0.0;       
        if ( bIncludeTEcom ) { ETot += this.TEcom; }        
        this.arrMol.forEach( ( mol, i ) => {
            ETot += this.TEradials[i] + this.TEtangents[i] + this.REradials[i] + this.REtangents[i];
            //console.log(this.TEradials[i], this.TEtangents[i], this.REradials[i], this.REtangents[i]);
        });
        return ETot;
    }
    
    // Need to account for A+B->C+D and  A+B->C
    // Just do the two body case for now.
    flip_product_orientation(reactantNames) {
        if ( this.arrMol[0].name != reactantNames[0] ) {
            return true;
        } else {
            return false;
        }
    }
    
    swap_molecules(i,j) {
        this.reset_component_energies();
        this.vecNormal.scale( -1.0 );
        
        var temp = this.arrMol[i];
        this.arrMol[i] = this.arrMol[j];
        this.arrMol[j] = temp;        
        temp = this.pRel[i];
        this.pRel[i] = this.pRel[j];
        this.pRel[j] = temp;
        temp = this.vRel[i];
        this.vRel[i] = this.vRel[j];
        this.vRel[j] = temp;        
    }
    
    /* More detailed collision checks. Either returns false or returns the position of contact. */
    // Must be called after creation.
    confirm_molecule_collision() {
        const mol1 = this.mol1, mol2 = this.mol2 ;
        const n1 = mol1.nAtoms, n2 = mol2.nAtoms;
        if ( 1 == n1 && 1 == n2 ) {
            this.pContact  = Vector2D.weighted_avg( mol2.get_size(), mol1.p, mol1.get_size(), mol2.p );
            this.vecNormal = mol2.p.subtract( mol1.p ).unit();
            return true;
        }
        const p1 = mol1.p, p2 = mol2.p, pa1 = new Vector2D(0,0), pa2 = new Vector2D(0,0);
        const o1 = [], o2 = [], r1 = [], r2 = [];
        mol1.atomOffsets.forEach( x => { o1.push( x.rotate(mol1.th) ) } );
        mol2.atomOffsets.forEach( x => { o2.push( x.rotate(mol2.th) ) } );
        mol1.atomRadii.forEach( x => { r1.push( x ) } );
        mol2.atomRadii.forEach( x => { r2.push( x ) } );        
        for (let i = 0; i < n1; i++) {
            pa1.x = p1.x+o1[i].x; pa1.y = p1.y+o1[i].y;
            for (let j = 0; j < n2; j++) {
                pa2.x = p2.x+o2[j].x; pa2.y = p2.y+o2[j].y;
                const pDiff = pa2.subtract( pa1 );
                const ratio = pDiff.norm() / ( r1[i] + r2[j] );
                if ( ratio < 1.0 ) {
                    this.pContact = Vector2D.weighted_avg( r2[j], pa1, r1[i], pa2 );
                    this.vecNormal = pDiff.unit() ;
                    //if ( ratio < 0.9) { console.log(`Warning, severe overlap between molecules have been detected between ${mol1.name} amd ${mol2.name} !`); }
                    return true;
                }
            }
        }
        return false;
    }

    /* General rigid body collision solver. Used when there is no reaction, but requires a pre-defined contant point and normal. */
    // This integrator has a tendency to add slight amounts of energy for monoatomic gases, and lose it for diatomic gases. The drift is very slight and within expectations of numerical integrators.
    // Based on maths of MyPhysicsLab. At: https://www.myphysicslab.com/engine2D/collision-en.html
    calculate_collision_impulse( elasticity, sep1P, sep2P ) {
        const mol1 = this.mol1, mol2 = this.mol2;
        const mass1 = mol1.get_mass(), rotI1 = mol1.get_rotI() ;
        const mass2 = mol2.get_mass(), rotI2 = mol2.get_rotI() ;
        const vInit12 = new Vector2D( 0, 0); //Relative initial volcity at point of contact.
        vInit12.x = mol1.v.x - mol1.om * sep1P.y - mol2.v.x + mol2.om * sep2P.y;
        vInit12.y = mol1.v.y + mol1.om * sep1P.x - mol2.v.y - mol2.om * sep2P.x;
        const vecN = this.vecNormal;
        //const vecN = mol2.p.subtract( mol1.p ).unit();
        const sRel = Vector2D.dot( vInit12, vecN ); // Should be >0.

        //const f = mass1 * mass2 / ( mass1 + mass2 ) ;
        //const impulse = f * sRel * (1 + elasticity) ;
        const w1 = (rotI1 != null ) ? sep1P.cross(vecN)**2.0/rotI1 : 0.0;
        const w2 = (rotI2 != null ) ? sep2P.cross(vecN)**2.0/rotI2 : 0.0;
        const f = 1.0 / ( 1.0/mass1 + 1.0/mass2 + w1 + w2 ) ;
        const impulse = f * sRel * ( 1.0 + elasticity ) ;
        return impulse;
    }
    
    /*
        Notes: These elastic collisions model favours redistribution of energy towards rotational energy.
        - instead of a 2:1 ratio...
        - systems vary between 3:2 ~ 5:4.
    */
    resolve_collision() {
        if ( undefined === this.pContact || undefined === this.vecNormal ) {
            const ret = this.confirm_molecule_collision();
            if ( !ret ) { throw "ERROR: Collisions resolution has been called where no collision is apparent!"; }
        }
        const mol1 = this.mol1, mol2 = this.mol2;        
        //console.log(`Elastic collision solver between ${mol1.name} and ${mol2.name}`);
        //console.log( mol1.get_total_energy(), mol2.get_total_energy() );
        
        const mass1 = mol1.get_mass(), rotI1 = mol1.get_rotI();
        const mass2 = mol2.get_mass(), rotI2 = mol2.get_rotI();
        const sep1P   = this.pContact.subtract( mol1.p ), sep2P = this.pContact.subtract( mol2.p );
        const impulse = this.calculate_collision_impulse( 1.0, sep1P, sep2P );
        mol1.v.sincr( -impulse / mass1, this.vecNormal  );
        mol2.v.sincr(  impulse / mass2, this.vecNormal  );
        if ( rotI1 != null ) { mol1.om -= impulse * sep1P.cross( this.vecNormal ) / rotI1 }
        if ( rotI2 != null ) { mol2.om += impulse * sep2P.cross( this.vecNormal ) / rotI2 }        
        
        // We want to move the molecules away so as to ensure that they defintiely no longer overlap.
        // This addresses the initial overlap problems, although the current form results in a little jitter.
        const dp  = mol2.p.subtract( mol1.p ); // Not directly opposite to vecN in general but are in case of circles.
        const sep = dp.norm(); 
        const pmove = mol1.get_size() + mol2.get_size() - sep;
        const f = 1.0 / ( 1.0/mass1 + 1.0/mass2 ) ;
        mol1.p.sincr( -1.0 * f / mass1 / sep * pmove, dp );
        mol2.p.sincr(  1.0 * f / mass2 / sep * pmove, dp );
        
        //console.log( mol1.get_total_energy(), mol2.get_total_energy() );
        //throw "pause";
    }    
}

// General class for handling reactions.
// Contains a self-consistency check to make sure all keys have defined values.
class reaction {
    constructor( args ) {
        //Set by the host handler.
        //Core variables required by all reactions.
        this.moleculeLibrary = args.moleculeLibrary ;        
        this.EActivation = args.EActivation;
        this.DeltaH = args.DeltaH;
        this.reactants = args.reactants;
        this.products  = args.products;
        this.productAngles = args.productAngles;
        this.angleReactionOffset = ( undefined === args.angleReactionOffset ) ? 0.0 : args.angleReactionOffset;
        // Variables helpful for reporting, side-effect of specificaiton..
        this.reactantNames = args.reactantNames;
        this.productNames  = args.productNames;

        // Parameters for decomposition reactions.
        // this.lifetimeActivated = args.lifetimeActivated ;
        // this.probRemain = 0.0 ;
        
        // Parameters for combination and general reactions where the impact angle matters.
        // this.reactantAngles = args.reactantAngles;
        // this.reactantAngleRanges = args.reactantAngleRanges;
        // Housekeeping parameters that should not be passed to the reaction class.
        // this.productAngleRanges = args.productAngleRanges;                
    }
    
    assert() {
        const keys = Object.keys(this);
        keys.forEach( k => {
            if (undefined === this[k]) { throw `Reaction object has not be properly set up! Key ${k} is undefined.` ;}
        });
    }

    debug() {
        console.log("= = Debug report of reaction instance:", this );
        for (const k in this) {  console.log( `  ${k}`, this[k] ); }
    }
    // General function to position two product molecules so that the center of mass remains the same.
    // theta here is either the orientation of the old molecule, or based on the diffference vector betwen the old COMs.
    // the division angle and new orientation angles are defined per reaction.
    assign_new_positions_2mol( mNew1, mNew2, pCOM, theta ) {
        //Slightly larger to mitigate recollision issues.
        const s    = 1.01 * (mNew1.size + mNew2.size);
        const mSum = mNew1.mass + mNew2.mass;
        const delta = new Vector2D( Math.cos( theta + this.angleReactionOffset ), Math.sin( theta + this.angleReactionOffset ) );
        mNew1.p.set_to( pCOM );
        mNew2.p.set_to( pCOM );        
        mNew1.p.sincr( s*mNew2.mass/mSum, delta );
        mNew2.p.sincr( -s*mNew1.mass/mSum, delta );
        if ( mNew1.rotI > 0 ) { mNew1.th = theta + this.angleReactionOffset + this.productAngles[0]; }
        if ( mNew2.rotI > 0 ) { mNew2.th = theta + this.angleReactionOffset + this.productAngles[1]; }
    }
    
    //Important for anything that depends on simulation time.
    update_time_dependence( dt ) {
        //Do nothing.
    };
    
    // Dummy functions that may be available depending on the type of reaction. Written here as documentation of expected usage
    // check_angles( mol1, mol2, pContact ) { return false; }
    check_angles_1mol( obj ) {
        const mol1 = obj.mol1, pContact = obj.pContact; 
        const dot1 = Vector2D.dot( pContact.subtract(mol1.p).unit(), Vector2D.UnitVector( mol1.th + this.reactantAngles[0] ) );
        if ( dot1 < this.reactantAngleRanges[0] ) { return false; }
        return true;
    }        
    
    check_angles_2mol( obj ) {
        const mol1 = obj.mol1, mol2 = obj.mol2, pContact = obj.pContact; 
        const dot1 = Vector2D.dot( pContact.subtract(mol1.p).unit(), Vector2D.UnitVector( mol1.th + this.reactantAngles[0] ) );
        if ( dot1 < this.reactantAngleRanges[0] ) { return false; }
        const dot2 = Vector2D.dot( pContact.subtract(mol2.p).unit(), Vector2D.UnitVector( mol2.th + this.reactantAngles[1] ) );
        if ( dot2 < this.reactantAngleRanges[1] ) { return false; }
        return true;
    }
    
    check_angles_general( obj ) {
        const mol1 = obj.mol1, mol2 = obj.mol2, pContact = obj.pContact; 
        if ( mol1.rotI > 0.0 ) {
            const dot1 = Vector2D.dot( pContact.subtract(mol1.p).unit(), Vector2D.UnitVector( mol1.th + this.reactantAngles[0] ) );
            if ( dot1 < this.reactantAngleRanges[0] ) { return false; }
        }
        if ( mol2.rotI > 0.0 ) {
            const dot2 = Vector2D.dot( pContact.subtract(mol2.p).unit(), Vector2D.UnitVector( mol2.th + this.reactantAngles[1] ) );
            if ( dot2 < this.reactantAngleRanges[1] ) { return false; }
        }
        return true;
    }    
    
    // get_available_energy( args ) { return undefined; }
    // check_energies( ETotAvail ) { return false; }
    // Overall pseudocode template for process reaction. Return null whenever reaction fails.
    process_reaction( args ) {
        // Check angles
        if ( false === this.check_angles( args ) ) { return null; }
        // Check energy
        var ETotAvail = this.get_available_energy( args );
        if ( ETotAvail < this.EActivation ) { return null; }
        // Check lifetime for decomposition reactions.
        if ( Math.random() > Math.exp( -args.dt / this.lifetimeActivated ) ) {return null; }        
        // All checks passed. Generate products.
        arrMolNew = [];
        //Add products.
        return arrMolNew;
    }
}

/* = = Subclasses for Decomposition reactions = = */
class reactionDecomposition extends reaction {
    // Decompose a polyatomic molecule into two polyatomic species each with rotational energy.
    // Hack and ignore the rigid-body details, but adopt equipartition instead.
    // Note that the order of molecules in the reaction matters. The first atom in the diatomic decomposes to the first atom in the product.
    constructor( args ) {
        super( args );
        this.lifetimeActivated = args.lifetimeActivated;
        this.probRemain = 0.0 ; // Start with a value that will instantly decomopse anything.
        this.assert();
    }

    update_time_dependence( dt ) { this.probRemain = Math.exp( -dt / this.lifetimeActivated ); }    
    //check_angles( obj ) { return true; }
    get_available_energy( obj ) { return obj.mol.get_rotational_energy() ; }
    
    process_reaction( obj ) {
        const mol = obj.mol;
        //Hack to borrow energy from the translational freedoms.
        var RE = mol.get_rotational_energy() ;
        if ( RE > this.EActivation && Math.random() < this.probRemain ) {
            //Successful reaction. Compute product kinetic properties and pass back to the simulation.
            // this.deltaH ; 
            const mNew1 = this.moleculeLibrary.create_molecule( this.products[0], { v:mol.v.copy() });
            const mNew2 = this.moleculeLibrary.create_molecule( this.products[1], { v:mol.v.copy() });
            this.assign_new_positions_2mol( mNew1, mNew2, mol.p, mol.th );
            
            const sysNew = new molSystem( [mNew1, mNew2] );            
            // Conserve angular momentum as well as linear momentum.            
            molSystem.distribute_angmoment_with_partial_energy_conservation( obj, sysNew );
            
            // Rescale energies so as to conserve energy from breaking the bond.
            // This breaks memomentum conservation.
            molSystem.rescale_velocities_by_energy_conservation( obj, sysNew, this.DeltaH );
            //console.log(`Decomposition reaction occurred: ${mol.name} -> ${mNew1.name} + ${mNew2.name}`);
            //sysNew.check_integrity();

            return [ mNew1, mNew2 ];

        } else {
            return null;
        }    
    }
}


/* = = Subclasses for Combination reactions = =
    Note that the two directionas are currently not equivalent.
    Normally, the excess energy terms are converted into the molecule's internal freedoms.
    They are then often re-emitted as photon(s) when electrons quickly hop back down from their initial excited states.
    We're currently hiding all this as if they don't exist.
    //const EExcess = ( obj.TETangent + obj.REinitial - mNew.get_rotational_energy() ) + (obj.TERadial - this.DeltaH);
*/
class reactionCombination extends reaction {
    // Combine two polyatomic molecules into one polyatomic species. Yes to initial angle check.
    // Note that the order of molecules in the reaction matters. The first partner in the collision should for the basis of theta = 0.
    constructor( args ) {
        super( args );
        this.reactantAngles = args.reactantAngles;
        this.reactantAngleRanges = args.reactantAngleRanges;
        this.assert();
    }
    
    process_reaction( obj ) {
        if ( false === this.check_angles_general( obj ) ) { return null; }
        const mol1 = obj.mol1, mol2 = obj.mol2 ;
        //Get pre-defined args.dp, args.mSum and args.pCOM
        obj.calculate_COM_configuration();
        obj.calculate_component_energies();
        const EAvail = obj.get_radial_components_sum();        
        if ( EAvail > this.EActivation ) {
            // Successful reaction. Compute product kinetic properties and pass back to the simulation.            
            const thSwitch = obj.flip_product_orientation( this.reactantNames ) ? Math.PI : 0.0 ;            
            const theta = Math.atan2( obj.pRel[1][1], obj.pRel[1][0] ) + thSwitch + this.angleReactionOffset;
            // Enact conservation of linear momentum here.
            const mNew = this.moleculeLibrary.create_molecule( this.products[0], { p: obj.p, v: obj.v, th: theta });
            // Enact conservation of angular momentum here. Sometimes this is nowhere near energy conserved!!
            obj.calculate_angular_momentums();
            mNew.om = ( obj.AngMomentumsExt + obj.AngMomentumsInt ) / mNew.rotI ;
            
            let TEnew = mNew.get_kinetic_energy(), REnew = mNew.get_rotational_energy() ;
            const EExcess = obj.get_total_energy( false ) - mNew.get_rotational_energy() - this.DeltaH;            
            let ERatio = ( REnew + TEnew + EExcess )/ (REnew + TEnew);            
            if (ERatio < 0) { console.log(REnew, TEnew, EExcess, ERatio); throw "Negative output energy found! Aborting."; }
            // Rescale energies so as to conserve energy and re-break conservation of momentum.
            mNew.rescale_velocities( Math.sqrt(ERatio) ); 
            // TEnew *= ERatio ; REnew *= ERatio;
            // mNew.v.scale( Math.sqrt( 2.0 * TEnew / (mNew.mass * mNew.v.norm2()) ) ); // 2 of 3 degrees of freedom to KE.
            // mNew.om = Math.sign( mNew.om ) * Math.sqrt( 2.0 * REnew / mNew.rotI ); // 1 of 3 degrees of freedom.
            // let Ediff = mNew.get_total_energy() - mol1.get_total_energy() - mol2.get_total_energy();
            // console.log( Ediff, this.DeltaH , ERatio ) ;

            //console.log(`Combination reaction occurred: ${mol1.name} + ${mol2.name} -> ${mNew.name}`);
            //const sysNew = new molSystem([ mNew ] );
            //sysNew.check_integrity();
            
            return [ mNew ];
        } else {
            //console.log(`Energy check fail. Energy: ${EAvail}`);
            return null;
        }    
    }
}

/*
    = = Subclass for transfer reactions where an atom hops from one molecule to another. = =
    Formulate this interaction process as an elastic collision, followed by the direct transfer of a single atom between hte two masses.
    
*/
// Example: NO2 + O3 <-> NO3 + O2.
class reactionTransfer extends reaction {
    constructor( args ) {
        super( args );
        this.reactantAngles = args.reactantAngles;
        this.reactantAngleRanges = args.reactantAngleRanges;
        this.assert();
    }
            
    process_reaction( obj ) {
        if ( false === this.check_angles_general( obj ) ) { return null; }
        //Get pre-defined args.dp, args.mSum and args.pCOM
        obj.calculate_COM_configuration();
        obj.calculate_component_energies();
        const EAvail = obj.get_radial_components_sum();
        const ETotalOld = obj.get_total_energy( true );
        if ( EAvail > this.EActivation ) {            
            // Successful reaction.
            //console.log(`Successful reaction between ${obj.mol1.name} and ${obj.mol2.name}`);
            // Check if the reactants need to be flipped so-as to generate a new product.
            if ( obj.flip_product_orientation( this.reactantNames ) ) {
                obj.swap_molecules( 0, 1 );
            }
            const mol1 = obj.mol1, mol2 = obj.mol2 ;
            
            // Compute momentum to be transferred via a delta-v directly between molecules.
            const dp = mol2.p.subtract( mol1.p ).unit();
            const dv = dp.scaled_copy( Vector2D.dot( mol1.v, dp ) - Vector2D.dot( mol2.v, dp ) );
            
            // Now enact an elastic collision.
            obj.resolve_collision();
            
            //NB: It's likely that some rotational momentunma nd energy might be lost when one partner becomes a monoatomic molecule.
            const mNew1 = this.moleculeLibrary.create_molecule( this.products[0], { p: obj.p.copy(), om: mol1.om });
            const mNew2 = this.moleculeLibrary.create_molecule( this.products[1], { p: obj.p.copy(), om: mol2.om });
            //console.log(`Products ${mNew1.name} and ${mNew2.name} have been created.`);
            this.assign_new_positions_2mol( mNew1, mNew2, obj.p, Math.atan2( obj.pRel[1][1], obj.pRel[1][0] ) );
            // mNew1.p.sincr( 1.01 * mNew1.get_size(), obj.pRel[1].unit() );
            // mNew2.p.sincr( 1.01 * mNew2.get_size(), obj.pRel[0].unit() );
            // const theta = Math.atan2( obj.pRel[1][1], obj.pRel[1][0] ) + this.angleReactionOffset;
            // mNew1.set_theta( this.productAngles[0] + theta );
            // mNew2.set_theta( this.productAngles[1] + theta );
            
            const sysNew = new molSystem( [ mNew1, mNew2 ]);
            //console.log( ETotalOld, sysNew.get_total_energy_from_mols() );
            if ( Number.isNaN ( sysNew.get_total_energy_from_mols() ) ) { throw "NaN detected!"; }
            //Enact COM-aligned momentum transfer post elastic-collision. Assume that mass is lost by the first molecule and equal mass is gained by the second.
            if ( true ) {
                const mDiff = mol1.mass - mNew1.mass;
                mNew1.v.set_to( Vector2D.weighted_sum( mol1.mass, mol1.v, -mDiff, dv ) );
                mNew1.v.scale( 1.0/mNew1.mass );
                mNew2.v.set_to( Vector2D.weighted_sum( mol2.mass, mol2.v, mDiff, dv ) );
                mNew2.v.scale( 1.0/mNew2.mass );
            } else {
                mNew1.v.set_to( mol1.v ) ; mNew2.v.set_to( mol2.v );
            }
            
            // Rescale energies so as to conserve energy instead of momentum.
            // console.log( ETotalOld, sysNew.get_total_energy_from_mols() );
            // if ( Number.isNaN ( sysNew.get_total_energy_from_mols() ) ) { throw "NaN detected!"; }            
            molSystem.rescale_velocities_by_energy_conservation( obj, sysNew, this.DeltaH );
            //console.log( ETotalOld, '-', this.DeltaH, '->', sysNew.get_total_energy_from_mols() );            
            //console.log(`Transfer reaction occurred: ${mol1.name} + ${mol2.name} -> ${mNew1.name} + ${mNew2.name}`);
            //sysNew.check_integrity();
            
            
            return [ mNew1, mNew2 ];
        } else {
            //console.log(`Energy check fail. Energy: ${EAvail}`);
            return null;
        }    
    }    
}




/*
    Interface function with overall program to setup example reactions.
    
    All energy values are divided by 10 relative to literature.
    This greatly speeds up reaction rates and helps lower the barriers relative to the lower dimensionality of the model.
    Notably, the internal vibrational energies are not available to us to help overcome activation barriers.    
    
    Usage:
    (1) Define an authoritative list of all species that can exist in the simulation.
    (2) Define every potential reaction, their energies, required impact geometries, and success probabilities.
    (3) Code will parse the arguments and attempt to determine the optimal reaction functions and where to slot in.
    
    Notes:
    (1) The forward and reverse reactions are created simultaneously, unless overriden by bDoReverse or bDoForward.
        - This is particularly useful for performance when the final product requires an unreasonable amount of energy to break apart,
          e.g. O + O -> O2. The reverse reaction has a DeltaH of ~ 500 kJ/mol.
    (2) Angles are defined relative to the vector from COM_1 to COM_2, i.e. p12.
    (2a) Decomposition reactions rely purely on the internal energy, i.e. the rotational energy in this simplified model.
    (2b) Whereas collision reactions compute the available energy from the collision encounter. Only the perpendicular components of the kinetic and rotational energy to the point of impact will be available for reaction.
    (2c) However, ALL energy will be redistributed to the products as neede to conserve energy. Conservation of momentum is sometimes used to compute initial vectors that will thenm be rescaled as needed.
    (3) lifetimeActivated is useful to describe an expected lifetime of an activated species in decomposition reactions:
        - It helps a decomposition reaction to not proceed immediately.
        - P(decay) = 1.0 - e ^ -( dt / t_lifetime ).
        - this defined in femtosecond, and multiplied by the dt parameter.
*/
function get_new_preset_gas_reactions( args ) {
    const type = ( undefined === args.type ) ? 'custom' : args.type;
    if ( undefined === args.moleculeLibrary ) { throw "Creation requires a molecule library to be given!"; }

    let gr = undefined;
    //let gr = new InteractionHandler();
    //gr.set_
    switch( type ) {
        case'none':
        case 'noble gas':
            gr = new InteractionHandler({ species: ["He","Ne","Ar","Kr","Xe"], moleculeLibrary: args.moleculeLibrary } );
            break;
        case 'atmosphere':
            gr = new InteractionHandler({ species: ["N₂","O₂","Ar","H₂O"], moleculeLibrary: args.moleculeLibrary } );
            break;
        case 'nitrogen dioxide':
            // DeltaH for the dissociation is ~13 kcal/mol with Ea of ~14 kcal/mol, or 54 & 59 kJ/mol.
            // The experimental self dissociation rate constant is ~ 1e6 per second at 298K. This is roughly 39 kJ/mol
            // Source: Ornellas et al., 2003. DOI: 10.1063/1.459256
            gr = new InteractionHandler({ species: ["NO₂","N₂O₄"], moleculeLibrary: args.moleculeLibrary } );
            gr.parse_input_reaction({
                //Participants.
                reactantNames: ["NO₂", "NO₂"],                
                productNames: ["N₂O₄"],
                //Impact geometry.
                reactantAngles:      [   0, 180 ], // Filled with 0.0 is not given.
                reactantAngleRanges: [ 180, 180 ], // Filled with 360 if not given. 
                productAngles:       [ 0 ],
                productAngleRanges:  [ 0 ],
                EActivation: 0.5, DeltaH: -3.9,
                lifetimeActivated: 2000,
                //Note: all arguments below reflect default setting
                //angleReactionOffset: 0.0,                
                //unitAngle: 'degrees',
                //bDoForward: true, bDoReverse: true,
            });  
            break;
            
        // Heat of formation dataset at 298 K: https://atct.anl.gov/Thermochemical%20Data/version%201.118/        
        // Full list: H = 216 ; O = 246 ; H2O = -239 ; OH = 37
        // Arrhenius activation energies have not been searched, and thus have been guessed.
        case 'combustion - H2 and O2':
            gr = new InteractionHandler( { species: ["H₂","O₂","H₂O","O•","H•","OH•"], moleculeLibrary: args.moleculeLibrary } );
            //Radical creation and elimination
            // H2 <-> 2H* DeltaH = 2* 216 kJ/mol.
            gr.parse_input_reaction({
                reactantNames: ["H•", "H•"], productNames: ["H₂"],
                EActivation: 0, DeltaH: -43.2, lifetimeActivated: 1000,
            });
            // O2 <-> 2O* DeltaH = 2 * 246 kJ/mol.
            gr.parse_input_reaction({
                reactantNames: ["O•", "O•"], productNames: ["O₂"],
                EActivation: 0, DeltaH: -49.2, lifetimeActivated: 1000,
                bDoReverse: false,
            });
            // OH + O <-> H2O.  DeltaH = -522 kJ/mol
            gr.parse_input_reaction({
                reactantNames: ["OH•", "H•"], productNames: ["H₂O"],
                EActivation: 0.5, DeltaH: -52.2, lifetimeActivated: 1000,
                reactantAngles:      [   0,   0 ], // Filled with 0.0 is not given.
                reactantAngleRanges: [ 240, 360 ], // Filled with 360 if not given. 
                productAngles:       [ 0 ],
                productAngleRanges:  [ 0 ],
                bDoReverse: false,
            });
            //  O2 + H <-> OH + O. Propagation collisions.
            // DeltaH ~ 67 kJ/mol            
            gr.parse_input_reaction({
                reactantNames: ["O₂", "H•"], productNames: ["OH•", "O•"],
                EActivation: 7.0, DeltaH: 6.7,
                reactantAngles:      [   0,   0 ], 
                reactantAngleRanges: [ 360, 360 ],
                productAngles:       [   0,   0 ],
                productAngleRanges:  [ 240, 360 ],
            });            
            //  H2 + O <-> OH + H.  DeltaH = 7kJ/mol
            gr.parse_input_reaction({
                reactantNames: ["H₂", "O•"], productNames: ["OH•", "H•"],
                EActivation: 1.0, DeltaH: 0.7,
                reactantAngles:      [   0,   0 ], 
                reactantAngleRanges: [ 360, 360 ],
                productAngles:       [ 180,   0 ],
                productAngleRanges:  [ 240, 360 ],
            });                        
            //  OH + H2 <-> H2O + H.  DeltaH = -60 kJ/mol 
            gr.parse_input_reaction({
                reactantNames: ["OH•", "H₂"], productNames: ["H₂O", "H•"],
                EActivation:  6.0, DeltaH: -6.0,
                reactantAngles:      [   0,   0 ], 
                reactantAngleRanges: [ 240, 360 ],
                productAngles:       [ 180,   0 ],
                productAngleRanges:  [ 240, 360 ],
                bDoReverse: false,
            });
            // Self transmission probability.  H2 + H <-> H + H2.  DeltaH = 0kJ/mol 
            gr.parse_input_reaction({
                reactantNames: ["H₂", "H•"], productNames: ["H•", "H₂"],
                EActivation: 10.8, DeltaH: 0.0,
            });
            
            break;
        case 'custom':
            gr = new InteractionHandler([]);        
            break;
        default:
            throw `Unrecognised gas composition preset ${type}!`;
    }
    return gr;
}