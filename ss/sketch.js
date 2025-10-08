// Simple 2D DEM code in JavaScript
// Improved version with better structure and comments using World class

// World object - contains all simulation data
let world;

// system variables
let ngw = 5; // number of grains in width
let ngh = 5; // number of grains in height
let rmin = 0.5E-3; // minimum radius
let rmax = 1.0E-3; // maximum radius
let density = 2700; // density
let scale = 45000; // scale for drawing
let vRand;

// target confining pressure
let pressure = 1.0;  // target confining pressure

// non-dimensional parameters
let kappa = 1000
let viscoRate = 0.98;
let InertialNumber = 1.0E-2;

// force law parameters
let kn = kappa * pressure // contact stiffness
let kt = kn;
let mu = 0.2; // friction coefficient

// precomputed variables
let ng = ngw * ngh; // total number of grains
let meanDiameter = rmin + rmax; // mean diameter
let meanMass = Math.PI * meanDiameter * meanDiameter * 0.25 * density;

// Calculate wall mass based on grain mass 
let wallMass = ngw * Math.PI * rmax * rmax * density;

// plot variables
let coordNumber = [];
let solidFrac = [];
let Vsolid;

// contact search variables
let dmax = 0.95 * rmin; // maximum distance for neighbor search
let nStepsUpsdate = 25; // number of steps for neighbor search update

// processing
let fnmax = 0.0; // maximum force (display purpose)
let nc = 0; // number of contacts

// loading
let stepBiaxStarts = 4000;
let axialVelocity = (ngh * meanDiameter) * InertialNumber * Math.sqrt(pressure / meanMass);

// pressure display variables
let currentPressureRight = 0.0;  // current pressure on right wall
let currentPressureTop = 0.0;    // current pressure on top wall
console.log("axialVelocity: " + axialVelocity);


// time flow
let dt = Math.PI * Math.sqrt(meanMass / kn) / 20 // time step
console.log("dt: " + dt);
let istep = 0; // step counter

function createGrains() {
    let grains = [];
    vRand = rmin / (200 * dt); // 5 * axialVelocity

    for (let iy = 0; iy < ngh; iy++) {
        for (let ix = 0; ix < ngw; ix++) {
            // position
            let x = rmax + ix * 2 * rmax;
            let y = rmax + iy * 2 * rmax;
            let pos = createVector(x, y);
            // radius
            let radius = random(rmin, rmax);
            // velocity
            let vel = createVector(random(-vRand, vRand), random(-vRand, vRand));

            // Create grain with material properties
            grains.push(new Grain(pos, vel, radius, density));
        }
    }

    let dt_critical = Math.PI * Math.sqrt(Math.PI * rmin * rmin * 1.0 * density / kn);
    dt = dt_critical / 10;

    return grains;
}


function graph(xcanvas, ycanvas, w, h, data, title) {
    stroke(255, 0, 0); // red
    line(xcanvas, ycanvas, xcanvas, ycanvas - h);
    line(xcanvas, ycanvas, xcanvas + w, ycanvas);

    if (data.length < 2) return
    let ymin = data[0];
    let ymax = data[0];
    for (let i = 1; i < data.length; i++) {
        if (data[i] < ymin) ymin = data[i];
        if (data[i] > ymax) ymax = data[i];
    }
    let Dx = w / (data.length - 1);
    let sc = h / (ymax - ymin);
    stroke(0, 0, 255);
    for (let i = 0; i < data.length - 1; i++) {
        line(xcanvas + Dx * (i - 1), ycanvas - sc * (data[i - 1] - ymin), xcanvas + Dx * i, ycanvas - sc * (data[i] - ymin));
    }
    fill(0);
    text(title + ' min: ' + round(ymin, 3) + ' max: ' + round(ymax, 3), xcanvas + 3, ycanvas - h + 10);
}

// P5.js setup and draw functions
function setup() {
    createCanvas(1000, 500);
    background(240);

    // Create grains
    let grains = createGrains();

    // Create World with all components
    world = new World(grains, dt, wallMass, kn);

    // Initialize neighbor list
    world.updateNeighborList(dmax);

    // Get solid volume
    Vsolid = world.getGrainVolume();
}

function draw() {
    background(240);

    // STEP 1: Reset forces first
    world.resetForces();

    // STEP 2: Compute all contact forces (grain-grain and grain-wall)
    world.computeWallForces(kn, viscoRate);
    fnmax = world.computeAllGrainForces(kn, kt, mu, viscoRate, dt);

    // STEP 3: Apply force-based wall control AFTER computing contact forces
    let rightWallArea = world.topWall.position.y - world.bottomWall.position.y; // height of right wall
    let topWallArea = world.rightWall.position.x - world.leftWall.position.x;   // width of top wall

    if (istep < stepBiaxStarts) {
        // ISOTROPIC COMPRESSION - both walls compress to target pressure
        world.controlWallWithPressure(world.rightWall, pressure, rightWallArea, 0.1);
        world.controlWallWithPressure(world.topWall, pressure, topWallArea, 0.1);
    } else {
        // BIAXIAL COMPRESSION - maintain pressure on right wall, compress top wall
        world.controlWallWithPressure(world.rightWall, pressure, rightWallArea, 0.1);
        world.controlWallWithVelocity(world.topWall, -axialVelocity);
    }

    // STEP 4: Update accelerations on grains (from contact forces)
    for (let i = 0; i < world.grains.length; i++) {
        world.grains[i].acc = p5.Vector.div(world.grains[i].force, world.grains[i].mass);
    }

    // STEP 5: Integrate motion (both grains and walls)
    world.grainTimestep();
    world.updateWalls(dt);

    // Prevent walls from moving outward beyond initial positions
    if (world.rightWall.position.x > ngw * 2 * rmax) {
        world.rightWall.position.x = ngw * 2 * rmax;
    }
    if (world.topWall.position.y > ngh * 2 * rmax) {
        world.topWall.position.y = ngh * 2 * rmax;
    }

    // Calculate current boundary pressures for display
    currentPressureRight = rightWallArea > 0 ? Math.abs(world.getRightWallForce()) / rightWallArea : 0;
    currentPressureTop = topWallArea > 0 ? Math.abs(world.getTopWallForce()) / topWallArea : 0;

    // Debug output every 1000 steps
    if (istep % 1000 == 0) {
        console.log(`Step ${istep}: Right P=${currentPressureRight.toFixed(3)}, Top P=${currentPressureTop.toFixed(3)}, Target=${pressure}`);
    }

    istep++; // increment step counter

    // Update neighbor list periodically
    if (istep % nStepsUpsdate == 0) {
        world.updateNeighborList(dmax);
    }

    // DISPLAY STUFF
    world.drawParticles(scale, vRand + abs(axialVelocity));
    world.drawWalls(scale);
    counts = world.drawContacts(fnmax, rmin);

    // compute and store data for plots using World data
    solidFrac.push(Vsolid / world.getGrainVolume());
    coordNumber.push(counts / ng);
    // draw graphs
    graph(550, 250, 500, 200, coordNumber, 'Coordination Number');
    graph(550, 500, 590, 200, solidFrac, 'Solid Fraction');

    // display step count
    fill(0);
    text('Step: ' + istep, 900, 20);
    text('Top pressure: ' + round(currentPressureTop, 3) + ' (Pa)', 700, 20);
    text('Right pressure: ' + round(currentPressureRight, 3) + ' (Pa)', 700, 40);



}
