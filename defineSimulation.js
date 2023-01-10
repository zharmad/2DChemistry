// Simulation Class to Handle overall simulation operations.
class Simulation {
    constructor() {
        //Molecule Handler
        this.nMolecules = 0;
        this.molecules = [];
        this.moletypeNames   = []; // This one is an array to keep track of order.
        this.moletypeCounts  = []; // All are arrays to make it easy to syn up with chart.js
        this.moletypeColours = [];
        this.nDegrees = 0;
        this.nMoleculesTarget = 0;
        
        //Objects to be initialised.
        this.moleculeLibrary = undefined;
        this.gasComp = undefined;
        this.gasReactions = undefined;
        
        //Simulation parameters. dt converts discrete timesteps into actual time in femtoseconds.
        this.timestep    = 0;
        this.timeElapsed = 0.0;
        this.dt          = 10.0 ;
        this.statsUpdateInterval  = 100;
        this.systemZeroInterval    = 100;
        //this.dt          = 50; //Basically in ~ picoseconds. average velocity of 2D oxygen molecule is ~300 ms^-1, or 30 px per timestep.       
        this.temperature = 300;
        this.bSet = false;
        this.bHeatExchange = true ;
        
        this.timeFactor = 0.001; // Convert time units from picoseconds, by default to femtoseconds.
        this.lengthScale  = 10; // Convert spatial units. Currently used to convert pm to pixels.        
        
        //Graphical interface
        this.graphicalContext = null ;
        this.xBounds = new Vector2D(0, 1);
        this.yBounds = new Vector2D(0, 1);
        this.refreshAlpha = 0.5 ;
        this.bDrawMolecules = true;
        
        // Accounting features.
        // Designed to be exportable directly to Chart.JS. Constructor does not know about its contents.
        // Default to current composition for now.
        this.chartDoughnutGr = undefined;        
        this.chartBarGr = undefined;       
        this.chartLineGr1 = undefined;
        this.chartLineGr2 = undefined;
        
        // TODO: Hookup with the dynamic component registry and create a single trajectory data object for use in plotting and export.
        /*
            For performance, dataframe syntax follows that of chart.js graphs data objects so that it can be directly sorted into Chart.js for plotting.
            User should decide how to format for specific purposes.
            - Line graphs: { label: 'T', data: [ [x1, y1], [x2, y2], ...], backgroundColor: 'black' }
            - Bar graphs: { label: 'KE', data: [ y1, y2, ...], backgroundColor: 'black' }
            - Pie graphs: {label: 'Count', data: [y1, y2, ...], backgroundColor: [c1, c2,...], hoverOffset: 4, borderWidth: 1}
            NB: Univariate graphs have a separate labels array to denote the X-values.
            NB: Bivariate graphs have a separate labels array to denote the legend.
        */
        this.dataFrame = undefined;
    }
    
    reset() {
        this.nMolecules = 0;
        this.molecules = [];
        this.nDegrees = 0;
        this.timestep    = 0;
        this.timeElapsed = 0.0;
        this.bSet = false;
        
        this.reset_data_frame();
    }
    
    // Internal dimensions will be in picometers. So will convert.
    setup_graphical_context(ctx, width, height, alpha) {
        if ( undefined === ctx ) { throw "Cannot setup graphical context without an argument!"; }
        if ( undefined === width ) { throw "Cannot setup graphical context without a width/height dimension!"; }
        if ( undefined === height ) { throw "Cannot setup graphical context without a width/height dimension!"; }
        if ( undefined === alpha ) { alpha = 0.4; }
        this.graphicalContext = ctx ;
        this.refreshAlpha = alpha ;
        this.set_world_boundaries( width, height );
    }
    
    // NB: this one function should never be called while the simulation is running.
    set_target_number_of_molecules( nMol ) { this.nMoleculesTarget = nMol; }
    
    set_world_temperature( T ) { this.temperature = T; }
    get_world_temperature() { return this.temperature; }
    set_bool_heat_exchange( bool ) { this.bHeatExchange = bool; }
    get_bool_heat_exchange() { return this.bHeatExchange; }
    set_bool_draw_molecules( bool ) { this.bDrawMolecules = bool; }
    set_world_length_scale( x ) { this.scaleLength = x ; }
    get_world_length_scale() { return this.scaleLength; }
    set_world_time_factor( x ) { this.timeFactor = x ; }
    get_world_time_factor() { return this.timeFactor; }  
    
    // Updating overall simulation speed will also affect some reaction properties.
    set_world_time_delta( x ) {        
        this.dt = x ;
        if ( undefined != this.gasReactions ) {
            this.gasReactions.update_time_dependence( x );
        }
    }
    
    get_world_time_delta() { return this.dt; }      
    set_statistics_update_interval( x ) { this.statsUpdateInterval = x; }
    set_system_zero_interval( x ) { this.systemZeroInterval = x; }
    set_world_boundaries( xMax, yMax ) {
        this.xBounds[1] = xMax * this.scaleLength;
        this.yBounds[1] = yMax * this.scaleLength;
    }
    
    update_values_from_globals() {
        this.set_world_time_factor( globalVars.timeFactor );
        this.set_world_time_delta( globalVars.timeDelta );
        this.set_world_length_scale( globalVars.lengthScale );
        this.set_world_boundaries( globalVars.worldWidth, globalVars.worldHeight );
        this.set_world_temperature( globalVars.temperature );
        this.set_bool_heat_exchange( globalVars.bHeatExchange );
        this.set_target_number_of_molecules( globalVars.numMolecules );
        this.set_statistics_update_interval( globalVars.statisticsUpdateInterval );
        
    }
    
    set_molecule_library( ml ) {
        this.moleculeLibrary = ml ;
    }
    
    set_gas_composition( gc ) {
        if ( this.gasComp != undefined ) { delete this.gasComp; }
        this.gasComp = gc;
        if ( 'count' === gc.type ) { this.nMoleculesTarget = gc.nTotal; }
    }

    set_gas_reactions( gr ) {
        if ( this.gasReactions != undefined ) { delete this.gasReactions; }
        this.gasReactions = gr;
        this.gasReactions.set_molecule_library( this.moleculeLibrary );
    }    
    
    get_random_position( margin ) {
        return random_2DUniform( this.xBounds.x + margin, this.xBounds.y - margin, this.yBounds.x + margin, this.yBounds.y - margin);
    }
    // mass is defined in amu. velocity as pm per fs.
    // KE is thus measured as 1 kJ/mol = 1000 kg m^2 / ( s^2 mol) = 1 amu pm^2 / fs^2.

    // Below are the only set of functions that should be allowed to update the molecule counts in the simulation.    
    // Second use case if an argument is given with 0 or less molecules to be created. Create the molecule type to be tracked later.
    // This is useful for pformation of products and intermediates.
    create_and_add_molecules(moletypeName, n) {
        if ( undefined === n ) { throw "Function create_and_add_molecules requires two arguments: molecule type and number of molecules to create!"; }
        const moltype = this.moleculeLibrary.get_entry( moletypeName );
        const j = this.moletypeNames.indexOf( moletypeName );
        if ( -1 === j ) {
            this.moletypeNames.push( moletypeName );
            this.moletypeCounts.push( 0 );
            this.moletypeColours.push( moltype.molColour );
            j = this.moletypeCounts.length - 1;
        }
        if ( n <= 0 ) { return; }
        
        let pInit = undefined, mol = undefined;
        for ( let i = 0; i < n; i++ ) {
            pInit = this.get_random_position( moltype.size );
            mol = this.moleculeLibrary.create_molecule( moltype, { bSample: true, T: this.temperature, p: pInit } );
            this.molecules.push( mol );
        }
        this.nMolecules += n; this.nDegrees += mol.nDegrees * n;
        this.moletypeCounts[j] += n;
    }
    
    // Add one molecule.
    add_molecule( mol ) {
        this.molecules.push( mol );
        this.nMolecules++; this.nDegrees += mol.nDegrees;
        const j = this.moletypeNames.indexOf( mol.name );
        this.moletypeCounts[j]++;
    }        
    // Move references to molecules around. Does update global statistics.
    remove_molecule( mol ) {
        const i = this.molecules.indexOf( mol );
        delete this.molecules[i];
        this.molecules[i] = this.molecules[this.nMolecules-1];
        this.molecules.pop();        
        this.nMolecules--; this.nDegrees -= mol.nDegrees;
        const j = this.moletypeNames.indexOf( mol.name );        
        this.moletypeCounts[j]--;
    }
    
    /* Generic initialiser that takes up the most memory */
    // Args are for specifying the size of the library and what molecules to build.
    initialise_molecules_libraries(args) {
        this.moleculeLibrary.reset_library();
        this.moleculeLibrary.add_all_known_molecule_types();
    }
    
    initialise_moletype_accounting( arrNames ) {
        const n = arrNames.length;
        this.moletypeNames   = arrNames;
        this.moletypeCounts  = [];
        this.moletypeColours = [];
        let colour = undefined;
        for ( let i = 0; i < n; i++ ) {
            this.moletypeCounts.push( 0 );
            colour = this.moleculeLibrary.get_entry( arrNames[i] ).molColour;
            this.moletypeColours.push( colour );
        }
    }
    
    /* The main builder  */
    build_simulation() {
        
        if ( undefined === this.gasComp ) { throw "Simulation instance has no defined composition to sample from!"; }
        this.molecules = [];
        this.moletypeNames = [];
        this.nMolecules = 0, this.nDegrees = 0;
        
        const obj = this.gasComp.get_components();
        this.initialise_moletype_accounting( this.gasComp.get_component_names() );
        switch ( this.gasComp.type ) {
            case 'count':
                obj.forEach( ([name,val]) => {
                    this.create_and_add_molecules( name, val );
                });
                break;
            case 'ratio':
                obj.forEach( ([name,val]) => {
                    const nCopies = Math.round( val * this.nMoleculesTarget );
                    this.create_and_add_molecules( name, nCopies );
                });
                break;
            default:
                throw `Simulation instance does not understand the type of Gas-composition: ${this.gasComp.type} !`;
        }
        
        //Create and manage the data frames.
        this.reset_data_frame();
        this.moletypeNames.forEach( name => {
            var colour = this.moleculeLibrary.get_molecule_color(name);
            // Switch from white to black for visual purposes.
            switch( colour ) {
                case 'rgb(255,255,255)':
                case 'white':
                case '#FFFFFF':
                    colour = 'rgb(0,0,0)';
                    break;
                default:
                    //Do nothing.
            }
            this.create_data_frame_entry( name, name, colour );
        });        

        //Synchronise and set up graph data.
        this.bSet = true;
        console.log(`Build complete. Created ${this.nMolecules} molecules.`);
    }   

    // This is meant to be all-encompassing so that it should pass a self-check.
    regenerate_simulation(args) {
        this.update_values_from_globals();
        //this.initialise_molecules_libraries(args);
        
        // Make sure the gas composition is standardised.
        this.gasComp.normalise();    
        
        this.build_simulation();
        this.draw_background(1.0);
        this.draw();
        
        this.timestep    = 0;
        this.timeElapsed = 0.0;
        //this.reset_line_graph();
        
        //Graph inital velocities.
        this.sync_all_graphs();
        this.update_all_graphs();
        
        this.run_self_check();
    }

    run_self_check() {
        if ( !this.bSet ) {throw "Simulation is not set!";}
        if ( 0 === this.nMolecules || [] === this.molecules) {throw "There are no molecules in the simulation!";}
        for (const mol of this.molecules) {
            if ( Number.isNaN( mol.p.x ) || Number.isNaN( mol.p.y ) ) {throw "A molecule has invalid positions!";}
            if ( Number.isNaN( mol.v.x ) || Number.isNaN( mol.v.y ) ) {throw "A molecule has invalid velocities!";}
        }
        if ( null === this.graphicalContext ) { throw "The graphical context has not been given!";}
        return true;
    }

    // This is responsible for wiping previous images.
    draw_background( alpha ) {
        if ( alpha === undefined ) { alpha = this.refreshAlpha; }
        const ctxLoc = this.graphicalContext ;
        ctxLoc.fillStyle = `rgba(255, 255, 255, ${alpha})`;
        ctxLoc.fillRect(0, 0, widthSim, heightSim);
        ctxLoc.lineWidth = 2;
        ctxLoc.strokeStyle = '#221100';
        ctxLoc.strokeRect(0, 0, widthSim, heightSim);        
        //Draw on the back of canvas.
        //context.globalCompositeOperation = 'destination-over';
        //Draw on the front again
        //context.globalCompositeOperation = 'source-over';
    }
    
    draw() {
        this.draw_background();
        for (const mol of this.molecules) {
            mol.draw(this.graphicalContext);
        }
    }
    
    // New strategy is to put all of the circles of the same colour in a single path.
    // This is best done with a molecule colouring scheme.    
    draw_all_new() {
        
        this.draw_background();
        //Collect every atom grouped by molecule colour.
        const xPos = {}, yPos = {}, rads = {}, colours = {};
        for ( const name of this.moletypeNames ) {
            xPos[name] = []; yPos[name] = []; rads[name] = [];
            colours[name] = this.moleculeLibrary.get_entry(name).molColour;
        }
        for (const mol of this.molecules) {            
            for (let i = 0; i < mol.nAtoms; i++) {
                const off = mol.atomOffsets[i].rotate( mol.th );                
                xPos[mol.name].push( (mol.p.x + off.x) / globalVars.lengthScale );
                yPos[mol.name].push( (mol.p.y + off.y) / globalVars.lengthScale );
                rads[mol.name].push( mol.atomRadii[i] / globalVars.lengthScale );
            }
        }
        
        //Draw one path for each colour.
        const ctxLoc = this.graphicalContext;        
        for ( const name of this.moletypeNames ) {
            ctxLoc.beginPath();                
            ctxLoc.fillStyle = colours[name];
            ctxLoc.lineWidth = 1;
            ctxLoc.strokeStyle = '#221100';
            const nCircles = rads[name].length;
            for ( let i = 0; i < nCircles; i++ ) {                
                ctxLoc.moveTo( xPos[name][i] + rads[name][i], yPos[name][i] );
                ctxLoc.arc( xPos[name][i], yPos[name][i], rads[name][i], 0, 2 * Math.PI );
            }
            ctxLoc.stroke();
            ctxLoc.fill();
        }
    }
    
    resolve_molecule_changes( arrAdd, arrDel ) {
        //console.log(`DEBUG at timestep ${this.timestep}: Molecules reacted!`);
        arrDel.forEach( m => { this.remove_molecule(m); });
        arrAdd.forEach( m => { this.add_molecule(m); });
    }
    
    // This O(n^2) step takes the most time. 32 of 80 seconds on last check for ~2000 molecule system.
    // TODO: Try emscripten -> Webassembly this piece of code.
    detect_potential_collisions() {
        const nMol = this.nMolecules;
        const xPos = new Float32Array(nMol);
        const yPos = new Float32Array(nMol);
        const sizes = new Float32Array(nMol);
        //Assign positions once and for all to stop pointer chasing.
        for ( let i = 0; i < nMol; i++ ) {
            const p = this.molecules[i].p;
            xPos[i] = p.vec[0];
            yPos[i] = p.vec[1];
            sizes[i] = this.molecules[i].size;
        }
        
        const molPairs = [];
        for (let i = 0; i < nMol-1; i++) {            
            for (let j = i + 1; j < nMol; j++ ) {
                var sepSq = (xPos[j]-xPos[i])*(xPos[j]-xPos[i]) + (yPos[j]-yPos[i])*(yPos[j]-yPos[i]);
                if ( sepSq < (sizes[j]+sizes[i])*(sizes[j]+sizes[i]) ) {
                    molPairs.push( [ this.molecules[i], this.molecules[j] ] );
                }
            }
        }        
        return molPairs;
    }
        
    resolve_all_potential_collisions( molPairs ) {    
        var ret = undefined, arrDel = [], arrAdd = [];
        const nPairs = molPairs.length;
        for( var i = 0; i < nPairs; i++ ) {            
            const m1 = molPairs[i][0]; const m2 = molPairs[i][1];            
            if ( m1.bIgnore || m2.bIgnore ) { continue; }
                
            // Offload each encounter to the interaction handler.
            var ret = this.gasReactions.process_molecule_encounter(m1, m2);
            
            // When a reaction has occurred, take note of the molecules to be created and deleted.
            if ( null != ret ) {
                m1.bIgnore = true; m2.bIgnore = true;
                arrDel.push( m1 ); arrDel.push( m2 );
                ret.forEach( m => { arrAdd.push( m ); });
            }
        }
        
        /* Account for additions and subtractions here. */
        if ( arrDel.length > 0 ) { this.resolve_molecule_changes( arrAdd, arrDel ); }
    }
    
    step() {
        this.timestep ++;
        const dt = this.dt;
        this.timeElapsed += dt * this.timeFactor ;
        //console.log( this.timestep, this.timeFactor );
        let nMol = this.nMolecules;

        //Resolve any spontaneous decomposition reactions
        //Placeholder for reaction list.
        var ret = undefined, arrDel = [], arrAdd = [];
        for (const mol of this.molecules) {
            ret = this.gasReactions.process_molecule_selfinteraction(mol);
            if ( null != ret ) {
                arrDel.push( mol );
                ret.forEach( m => { arrAdd.push( m ); });
            }
        } 
        /* Account for additions and subtractions here. */
        if ( arrDel.length > 0 ) { this.resolve_molecule_changes( arrAdd, arrDel ); }        
                
        // Simple movement
        //this.draw();
        if( this.bDrawMolecules ) { this.draw_all_new(); }        
        for (const mol of this.molecules) {
            mol.update_position( dt );
        }
        
        // Detect and resolve collisions
        const molPairs = this.detect_potential_collisions();
        if ( molPairs.length > 0 ) { this.resolve_all_potential_collisions( molPairs ); }
        
        //Run self check.
        //console.log(`Checking integrity at step ${this.timestep}`);        
        // for ( const mol of this.molecules ) {
            // if ( Number.isNaN( mol.v.x ) || Number.isNaN( mol.v.y ) ) {
                // mol.debug();
                // throw "NaN velocity values have been detected!";
            // }
        // }
        
        // Resolve any collision with the walls.
        for (const mol of this.molecules) {
            //Shift molecules            
            this.process_wallBounce(mol);
        }       
        
        // inform downstream updaters
        if ( this.timestep % this.statsUpdateInterval == 1 ) {
            this.push_data_frame();
            this.update_all_graphs();
        }
        
        //Zero all angular momentum to reduce ice cude phenomenon. This now breaks strict energy conservation of the system.
        if ( this.timestep % this.systemZeroInterval == 1 ) {
            this.zero_total_momentum();
            // const sysP = this.get_linear_momentum(); 
            // sysP.scale( 1.0/this.nMolecules );
            // console.log(`Reporting at time step ${this.timestep} ( ${this.nMolecules} molecules )` );
            // console.log("System linear momentum per molecule:", sysP.x, sysP.y );
            // console.log("System angular momentum per molecule:", this.get_angular_momentum() / this.nMolecules );
        }        
    }

    debug() {
        for (const mol of this.molecules) {
            var n = [mol.p[0], mol.p[1], mol.v[0], mol.v[1], mol.th, mol.om]
            var out = n.map(n => parseFloat(n.toPrecision(3)));
            console.log( mol.name, out);
        }
    }

    // TO-DO: switch to the rigid body collider with the wall.
    // const vel1PInit = vInit1 + scalar_cross( om, sep1P );
    // const w = (rotI != null ) ? sep1P.cross(vecN)**2.0/rotI : 0.0;
    // const f = 1.0 / ( 1.0/mass + w ) ;
    // const impulse = f * vel1PInit.dot(vecN) * (1 + elasticity) ;
    process_wallBounce(mol) {
        const xBounds = this.xBounds, yBounds = this.yBounds;
        const s = mol.size, p = mol.p, v = mol.v ;
        let bCollide = false ;
        if ( p.x - s < xBounds.x ) {
            v.x = -v.x;
            p.x += 2.0*( s + xBounds.x - p.x );
            bCollide = true;
        }
        if ( p.x + s > xBounds.y) {
            v.x = -v.x;
            p.x += 2.0*( xBounds.y -p.x - s );
            bCollide = true;
        }

        if ( p.y - s < yBounds.x) {
            v.y = -v.y;
            p.y += 2.0*( s + yBounds.x - p.y );
            bCollide = true;
        }        
        if ( p.y + s > yBounds.y) {
            v.y = -v.y;
            p.y += 2.0*( yBounds.y - p.y - s );
            bCollide = true;
        }
        
        /*
        Resample energies upon contact with the outside world. The default loses energy over time because faster molecules collide with the walls more often. This results in net energy tranfer to the outside. A correction constant is thus determined numerically from atmospheric samples.
        A values of 1.38 leads to the following:
            - 200 molecules in 303nm^2 box gives 292K.
            - 500 molecules in 303nm^2 box gives 300K.
        This confirms that the constant will depend on density, i.e. collision rate between molecules.
        */
        if ( this.bHeatExchange && bCollide ) {
            // mol.resample_speed( this.temperature * 1.26 );
            // mol.resample_omega( this.temperature * 1.26 );
            mol.resample_speed( this.temperature * 1.38 );
            mol.resample_omega( this.temperature );            
        }
    }
    
    /* General analysis functions*/
    get_volume() {
        return (this.xBounds.y-this.xBounds.x)*(this.yBounds.y-this.yBounds.x);
    }
    get_total_energy() {
        let ETot = 0.0; 
        for (const mol of this.molecules) {
            ETot += mol.get_total_energy();
        }
        return ETot;        
    }
    get_total_kinetic_energy() {
        let KE = 0.0;
        for (const mol of this.molecules) {
            KE += mol.get_kinetic_energy();
        }
        return KE;
    }
    get_total_rotational_energy() {
        let RE = 0.0;
        for (const mol of this.molecules) {
            RE += mol.get_rotational_energy();
        }
        return RE;
    }        
    get_mean_velocity() {
        let v = 0.0;
        for (const mol of this.molecules) {
            v += mol.v.norm();
        }
        return v/this.nMolecules;
    }       
    get_temperature() {
        // Note: the minus 3 comes from the constraints on setting the center of mass and rotation to zero.
        const totE = this.get_total_energy();
        return totE / ( 0.5 * (this.nDegrees - 3) * 8.314 * this.timeFactor**2.0 * 1000 ) ;
        //const totE = this.get_total_kinetic_energy();
        //return totE / ( this.nMolecules * 8.314 * this.timeFactor**2.0 * 1000 ) ;
    }
    get_centre_of_mass() {
        let pCent = new Vector2D(0,0), mTot = 0.0;
        for (const mol of this.molecules) {
            pCent.sincr( mol.mass, mol.p );
            mTot += mol.mass;
        }
        pCent.scale( 1.0 / mTot );
        return pCent;
    }
    get_angular_momentum() {
        let pCent = this.get_centre_of_mass();
        let L = 0.0;
        for ( const mol of this.molecules ) { L += mol.get_angular_momentum( pCent ); }
        return L;
    }
    get_linear_momentum() {
        let P = new Vector2D( 0, 0 );
        for ( const mol of this.molecules ) { P.sincr( mol.mass, mol.v); }
        return P;
    }

    // Conduct an iso-kinetic energy redistribution of total momentums and energies.
    // This transfer the system's rotational and linear energies back into the kinetic energies of individual molecules.
    // Note that the rotational energies if indviidual atoms have remained untouched for now. This may need to be redistributed as well.
    zero_total_momentum() {
        //if ( undefined === bLinear ) { bLinear = true; }
        //if ( undefined === bAngular ) { bAngular = true; }
        let pCent = this.get_centre_of_mass();
        let PTot = new Vector2D( 0, 0 );
        let MTot = 0, LTot = 0.0, ITot = 0.0, ETotOld = 0.0;
        for (const mol of this.molecules) {
            MTot += mol.mass ;
            PTot.sincr( mol.mass, mol.v );
            ITot += mol.mass * mol.p.subtract( pCent ).norm2();            
            LTot += mol.get_angular_momentum( pCent );
            ETotOld += mol.get_total_energy();
        }
        const sysOm = LTot / ITot; const sysV = PTot.scaled_copy( 1.0 / MTot );
        let ETotNew = 0.0;
        //console.log(`Zeroing system angular momentum about ( ${pCent.x}, ${pCent.y} ) with angular momentum ${LTot} and rotation ${sysOm} ...`);
        for (const mol of this.molecules) {
            const vRect = Vector2D.scalar_cross( -sysOm, mol.p.subtract( pCent ) );
            vRect.decr( sysV );
            mol.v.incr( vRect );
            ETotNew += mol.get_total_energy();
        }
        // Energy redistribution
        const ratio = Math.sqrt( ETotOld/ETotNew );
        for (const mol of this.molecules) {
            mol.v.scale( ratio );
            mol.set_omega( mol.om*ratio );
        };
        
        //debugging
        // const PFinal = this.get_linear_momentum();
        // const LFinal = this.get_angular_momentum();
        // const ETotFinal = this.get_total_energy();        
        // console.log( "Old momentum and energies:", PTot.x, PTot.y, LTot, ETotOld);
        // console.log( "New momentum and energies:", PFinal.x, PFinal.y, LFinal, ETotFinal );
    }
    
    /* Reaction handling functions */
    test_decomposition( mol ) {
        const totKE = mol.get_kinetic_energy() + mol.get_rotational_energy();        
        console.log( totKE );
        if ( totKE > 12.6 ) {
            // Create decomposed molecules and delete this one.
        }
    }
    
    /*
        All functions for handling data
        - Line graphs: { label: 'T', data: [ [x1, y1], [x2, y2], ...], backgroundColor: 'black' }
        - Bar graphs: { label: 'KE', data: [ y1, y2, ...], backgroundColor: 'black' }
        - Pie graphs: {label: 'Count', data: [y1, y2, ...], backgroundColor: [c1, c2,...], hoverOffset: 4, borderWidth: 1}
    */
    
    reset_data_frame() {
        delete this.dataFrame;
        this.dataFrame = {};
        this.create_data_frame_entry( 'temperature', 'T' );
    }
    
    create_data_frame_entry( key, label, BGColour ) {
        if ( undefined === label ) { label = key; }
        if ( undefined === BGColour ) { BGColour = 'black'; }
        this.dataFrame[key] = { label: label, data: [], backgroundColor: BGColour }        
    }
    
    // 
    push_data_frame() {
        const t = this.timeElapsed;
        this.dataFrame['temperature'].data.push( [ t, this.get_temperature() ] );
        const n = this.moletypeNames.length;       
        let name = undefined, count = undefined;
        for ( let i = 0; i < n; i++ ) {
            name  = this.moletypeNames[i];
            count = this.moletypeCounts[i];
            this.dataFrame[name].data.push( [ t, count ] );
        }
        this.dataFrame['temperature'].data.push( [ t, this.get_temperature() ] );
    }
        
    /* Analysis functions for graphing in Chart.JS */
    // The synchronisation functions simply point the data objects within the charts to their equivalents within the simulation so that they are automatically updated with the timestep
    // This is so that when the times come to update, we can just call update_graph().
    
    // Usage scenario: user passes an array of 
    sync_doughnut_graph() {
        const d = this.chartDoughnutGr.data ;        
        d.labels = this.moletypeNames ;
        d.datasets[0].data = this.moletypeCounts ;
        d.datasets[0].backgroundColor = this.moletypeColours ;
    }
    
    sync_bar_graph( arr ) {
        if ( undefined === arr ) { arr = ['velocities']; }
        const d = this.chartBarGr.data ; const n = arr.length;
        for ( let i = 0; i < n; i++ ) {
            switch ( arr[i] ) {
                case 'velocities':
                    d.datasets[i] = {
                        label: 'velocities (pm/fs)',
                        data: [],
                        backgroundColor: 'rgb(72,  96, 128)'
                    }
                    break;
                default:
                    throw `Bar graph type ${arr[i]} is not recognised!`;
            }
        }
    }
  
    //Temperature
    sync_line_graph_1() {
        this.chartLineGr1.data.datasets = [ this.dataFrame['temperature'] ]
    }
    
    //Copunts of individual molecule types
    sync_line_graph_2( arrEntries ) {
        if ( undefined === arrEntries ) { arrEntries = this.moletypeNames; }
        //this.chartLineGr1.data.labels = arrLabels;        
        const d = this.chartLineGr2.data.datasets = [];
        const n = arrEntries.length;
        let entry = undefined;
        for ( let i = 0; i < n; i++ ) {
            entry = this.dataFrame[ arrEntries[i] ];
            if ( undefined === entry ) { throw `Unrecognised entry for simulation data frame! ${arrEntries[i]}`; }
            d.push( entry );
        }
    }
    
    sync_all_graphs() {
        this.sync_doughnut_graph(); //Molecule composition for now.        
        this.sync_bar_graph(); //Atom velocities for now.
        this.sync_line_graph_1();
        this.sync_line_graph_2();
    }    
    
    inventory_bar_graph() {
        //this.reset_bar_graph();
        const arrVal = [], nBins = 20;
        //var maxVal = 0.0;
        var meanVal = 0.0;
        var temp = undefined ;     
        for (const mol of this.molecules) {
            //temp = mol.get_kinetic_energy() + mol.get_rotational_energy();
            //temp = mol.get_kinetic_energy();
            temp = mol.v.norm();
            arrVal.push( temp );
            //maxVal = Math.max( maxVal, temp );
            meanVal += temp;
        }
        const maxVal = 3.0* meanVal / this.nMolecules;
        /* Bucket into bins and plot */
        const d = this.chartBarGr.data;
        d.labels = Array.apply(null, Array(nBins)).map(function (x, i) { const a = i * maxVal / nBins ; return a.toFixed(2); }) 
        d.datasets[0].data = Array.apply(null, Array(nBins)).map(function (x, i) { return 0; }) 
        arrVal.forEach( v => {
            temp = Math.floor( nBins * v / ( maxVal * 1.001 ) );
            d.datasets[0].data[temp]++;
        });
        if ( this.chartBarGr != undefined ) { this.chartBarGr.update(); }
    }    
    
    
    update_doughnut_graph() {
        this.chartDoughnutGr.update();
    }
    
    update_line_graphs() {
        if ( this.chartLineGr1.bUpdate ) { this.chartLineGr1.update(); }
        if ( this.chartLineGr2.bUpdate ) { this.chartLineGr2.update(); }
    }
    
    update_all_graphs() {
        this.update_doughnut_graph();
        this.inventory_bar_graph();
        this.update_line_graphs();
        //this.update_line_graph( 0, this.timeElapsed, this.get_temperature() );
        //this.update_line_graph( 0, this.timeElapsed, this.get_total_kinetic_energy() );
        //this.update_line_graph( 1, this.timeElapsed, this.get_total_rotational_energy() );
    }
}