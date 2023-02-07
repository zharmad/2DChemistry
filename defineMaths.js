// functions to generate random numbers 

function random(min, max) {
    const num = Math.random()*(max-min)+min;
    return num;
}

function randomInt(min, max) {
    const num = Math.floor(Math.random() * (max - min + 1)) + min;
    return num;
}

function random_1DGaussian(s) {
    if ( undefined === s ) { s = 1.0; }
    let a = random_chi2D() * s ;
    let b = 2.0 * Math.PI * Math.random() ;
    return a * Math.cos(b);
}

// Taken from https://stackoverflow.com/questions/12556685/is-there-a-javascript-implementation-of-the-inverse-error-function-akin-to-matl
// Quick approximation from Abramowitz and Stegun (1964). Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables.
function erfinv(x) {
  // maximum relative error = .00013
  const a  = 0.147 ;
  //if (0 == x) { return 0 }
  const b = 2/(Math.PI * a) + Math.log(1-x**2)/2 ;
  const sqrt1 = Math.sqrt( b**2 - Math.log(1-x**2)/a ) ;
  const sqrt2 = Math.sqrt( sqrt1 - b ) ;
  return sqrt2 * Math.sign(x) ;
}

function gaussian( x, m, s ) { return Math.exp( -0.5*(x-m)*(x-m)/(s*s) ); }
/* Note: something seems wrong in the sampling of these velocities. The temperatuer returned is never quite the input. */

// These factors are to be used to scale the random distrubtions so as to replicate the relevant chi distributions.
// Apply to Random_2DGaussian to generate the relevant Maxwell velocity disturbtion for a particle.
// R = 8314.46
// Standard.
function scale_factor_Maxwell(T, mass) { return Math.sqrt( 8314.46 * T / mass ); }

// Suipplementary functions that compute the mean velocity of a particle from the PDF
// 1D = sqrt( 2/pi* RT/mass )
//function scale_factor_Maxwell1D(T, mass) { return Math.sqrt( 5293.14963255936 * T / mass ); }
// 2D ~= sqrt( 5/pi* RT/mass ). This is numerically fitted instead.
//function scale_factor_Maxwell2D(T, mass) { return Math.sqrt( 13060.317403376966 * T / mass ); }
// 3D = sqrt( 8/pi* RT/mass )
//function scale_factor_Maxwell3D(T, mass) { return Math.sqrt( 21172.59853023744 * T / mass ); }

// Returns the scalar speed of a particle in 2D motion.
function random_chi2D() { return Math.sqrt( -2.0 * Math.log( 1.0 - Math.random() ) ); }
function random_chi1D() { return erfinv( 2.0 * Math.random() - 1.0 ) * 1.4142; }

// Returns speed in pm ps^1. Mass in amu. Temperature in Kelvin. 1.4142 is a squareroot 2 for two degrees of freedom.
function random_speed2D( T, mass ) { return random_chi2D() * scale_factor_Maxwell( T, mass ); }
function random_velocity2D( T, mass ) { return random_2DGaussian( scale_factor_Maxwell( T, mass ) ); }

function random_speed1D( T, mass ) { return random_chi1D() * scale_factor_Maxwell( T, mass ); }
// function random_speed1D( T, mass ) {
    // let a = random_chi2D() * scale_factor_Maxwell( T, mass );
    // let b = 2.0 * Math.PI * Math.random() ;
    // return a * Math.cos(b);
// }

//random_chi1D needs a tabular implemention.
//  return Math.sqrt(2)*erfinv(2.0*Math.random()-1);

// Box Muller transform, but eliminate outliers to reduce edge cases.
// https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
// This is the 2D chi distribution equivalent.
function random_2DGaussian(s) {    
    if ( undefined === s ) { s = 1.0; }
    let a = random_chi2D() * s;
    //a = Math.min( 4.0, a ) * s ;
    let b = 2.0 * Math.PI * Math.random() ;
    return new Vector2D( a * Math.cos(b), a * Math.sin(b) );
}


function random_2DUniform(xmin, xmax, ymin, ymax) {
    let x = Math.random()*(xmax-xmin)+xmin;
    let y = Math.random()*(ymax-ymin)+ymin;
    return new Vector2D(x,y);
}

// function to generate random color
function randomRGB() {
    return `rgb(${random(0, 255)},${random(0, 255)},${random(0, 255)})`;
}

function array_average( v ) {
    const n = v.length;
    let tot = 0.0;
    for ( let i = 0; i < n; i++ ) {
        tot += v[i]
    }
    return tot/n;
}

// 2D Vector functions
class Vector2D {
    constructor( x, y ) { this.vec = [x,y] ; }
    
    // Using a getter/setter to make a couple of useful aliases
    get x() { return this.vec[0]; }
    set x(a) { this.vec[0] = a; }
    get y() { return this.vec[1]; }
    set y(a) { this.vec[1] = a; }    
    get 0() { return this.vec[0]; }
    set 0(a) { this.vec[0] = a; }
    get 1() { return this.vec[1]; }
    set 1(a) { this.vec[1] = a; }         
    
    // Static functions. Always returns a new vector.
    static weighted_sum( f1, v1, f2, v2 ) {
        return new Vector2D( f1*v1.vec[0]+f2*v2.vec[0], f1*v1.vec[1]+f2*v2.vec[1] );
    }
    static weighted_avg( f1, v1, f2, v2 ) {
        const fTot = f1+f2;
        return new Vector2D( (f1*v1.vec[0]+f2*v2.vec[0])/fTot, (f1*v1.vec[1]+f2*v2.vec[1])/fTot );
    }    
    static scalar_cross( s, v ) { return new Vector2D( -s*v.vec[1], s*v.vec[0] ); }
    static cross( v1, v2 ) { return v1.vec[0] * v2.vec[1] - v1.vec[1] * v2.vec[0]; }
    static dot( v1, v2 ) { return v1.vec[0] * v2.vec[0] + v1.vec[1] * v2.vec[1]; }
    static dist( v1, v2 ) { return Math.sqrt( Vector2D.dist2( v1, v2 ) ); }
    static dist2( v1, v2 ) { return (v2.vec[0]-v1.vec[0])*(v2.vec[0]-v1.vec[0]) + (v2.vec[1]-v1.vec[1])*(v2.vec[1]-v1.vec[1]) };
    static atan2( unit ) { return Math.atan2( unit.vec[1], unit.vec[0] ); }
    static rotate( v, th ) {
        let x = v.vec[0] * Math.cos(th) - v.vec[1] * Math.sin(th) ;
        let y = v.vec[0] * Math.sin(th) + v.vec[1] * Math.cos(th) ;
        return new Vector2D( x, y );
    }
    static duplicate( v ) { return new Vector2D( v.vec[0], v.vec[1] ); }
    
    static UnitVector( theta ){ return new Vector2D( Math.cos(theta), Math.sin(theta) ); }
    
    // Object functions. 
    debug() { console.log(`Values: [${this.vec[0]},${this.vec[1]}]`); }
    
    copy() { return Vector2D.duplicate( this ); }
    set_to( v ) { this.vec[0] = v.vec[0] ; this.vec[1] = v.vec[1]; }

    scaled_copy( s ) { return new Vector2D( s*this.vec[0], s*this.vec[1] ); }
    scale( s ) { this.vec[0] *= s; this.vec[1] *= s; }    
    incr( v ) { this.vec[0] += v.vec[0] ; this.vec[1] += v.vec[1]; }
    decr( v ) { this.vec[0] -= v.vec[0] ; this.vec[1] -= v.vec[1]; }        
    rotate( th ) {
        let x = this.vec[0] * Math.cos(th) - this.vec[1] * Math.sin(th) ;
        let y = this.vec[0] * Math.sin(th) + this.vec[1] * Math.cos(th) ;
        this.vec[0] = x; this.vec[1] = y;
    }

    sincr( s, v ) {
        this.vec[0] += s * v.vec[0];
        this.vec[1] += s * v.vec[1];
    }
    
    add( v ) {
        return new Vector2D(this.vec[0] + v.vec[0], this.vec[1] + v.vec[1]);
    }
    subtract( v ) {
        return new Vector2D(this.vec[0] - v.vec[0], this.vec[1] - v.vec[1]);
    }
    
    norm2() { return this.vec[0] * this.vec[0] + this.vec[1] * this.vec[1]; }
    norm() { return Math.sqrt( this.norm2() ); }
    
    unit() {
        const d = this.norm();
        return new Vector2D(this.vec[0]/d, this.vec[1]/d);
    }
    
    dot( v ) {
        return this.vec[0]*v.vec[0]+this.vec[1]*v.vec[1];
    }
    cross( v ) {
        return this.vec[0]*v.vec[1]-this.vec[1]*v.vec[0];
    }

    dist2( v ) {
        return (this.vec[0]-v.vec[0])**2.0 + (this.vec[1]-v.vec[1])**2.0 ;
    }
    dist( v ) {
        return Math.sqrt( this.dist2( v ) );
    }
    
    component_along( v ) {
        let ret = Vector2D.duplicate( v );
        ret.scale( this.dot(v) / v.norm2() );
        return ret;
    }
}

//fetch("distCheck.wasm").then(bytes=> bytes.arrayBuffer()).then(mod=> WebAssembly.compile(mod)).then(module=> {return new WebAssembly.Instance(module)}).then( instance => { check_collision_wasm = instance.exports.check_collision;});


/*
    2D Maxwell-velocity distribution CDF lookup table.
*/
// function CDF_chi2(x) {
    // return 1.0-Math.exp(-0.5*x*x);
// }

// class LookupTableSampler {
    // constructor(func, min, max, step) {
        // this.func = func ;
        // this.min = min ;
        // this.max = max ;
        // this.step = step ;
        
        // this.x = [];
        // this.y = [];
        // this.construct_lookup_table();
    // }
    
    // construct_lookup_table() {
        // for (let i = this.min; i <= this.max ; i += this.step) {
            // this.x.push(i);
            // this.y.push( this.func(i) );
        // }
    // }

    // sample() {
        // const yval = Math.random();
        
    // }
// }
// lookupMW = new LookupTableSampler( CDF_chi2, 0.0, 4.0, 0.1 );
