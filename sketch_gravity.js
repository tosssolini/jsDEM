// Simple 2D DEM code in JavaScript
// Improved version with better structure and comments using World class

// World object - contains all simulation data
let world;

// Scale for drawing
let scale = 20000;

// system variables
let ngw = 10; // number of grains in width
let ngh = 10; // number of grains in height
let rmin = 0.5E-3; // minimum radius
let rmax = 1.0E-3; // maximum radius
let vRand;
let ng = ngw * ngh; // total number of grains
let meanDiameter = rmin + rmax; // mean diameter

// target confining pressure
let pressure = 1.0;  // target confining pressure

// non-dimensional parameters
let kappa = 100000
let viscoRate = 0.98;
let InertialNumber = 1.0E-2;
let damping = 0.01; // damping coefficient

// Grain parameters
let kn = kappa * pressure // contact stiffness
let kt = kn;
let mu = 0.2; // friction coefficient
let density = 2700; // density

// Wall parameters
let knWall = kn; // wall stiffness
let mWall = ngw * Math.PI * rmax * rmax * density; // wall mass based on grain mass

// Timestep first guess based on mean mass
let meanMass = Math.PI * meanDiameter * meanDiameter * 0.25 * density;
let dt = Math.PI * Math.sqrt(meanMass / kn) / 20 // time step

// contact search variables
let dmax = 0.95 * rmin; // maximum distance for neighbor search
let nStepsUpsdate = 1; // number of steps for neighbor search update

// other parameter declarations
let fnmax = 0.0; // maximum force (display purpose)
let nc = 0; // number of contacts
let istep = 0; // step counter
let coordNumber = [];
let solidFrac = [];
let Vsolid;

function createGrains() {
    // Create grains with random positions and velocities on a grid without overlap
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

    return grains;
}

function createWalls(xmin, xmax, ymin, ymax) {
    let walls = [];
    walls.push(new Wall(createVector(xmin, ymin), createVector(1, 0), Infinity, knWall)); // left
    walls.push(new Wall(createVector(xmax, ymin), createVector(-1, 0), mWall, knWall)); // right
    walls.push(new Wall(createVector(xmin, ymin), createVector(0, 1), Infinity, knWall)); // bottom
    walls.push(new Wall(createVector(xmin, ymax), createVector(0, -1), mWall, knWall)); // top
    return walls;
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

    // Create walls
    let walls = createWalls(0, ngw * 2 * rmax, 0, ngh * 2 * rmax);

    // Create World with all components
    world = new World(grains, walls, dt, kn, kt, mu, viscoRate, damping);

    // Initialize neighbor list
    world.updateNeighborList(dmax);

    // Get solid volume
    Vsolid = world.getGrainVolume();
    console.log('Solid volume: ' + Vsolid);
}

function draw() {
    background(240);

    // Reset forces first
    world.resetForces();

    // Compute all contact forces (grain-grain and grain-wall)
    world.computeWallForces();
    fnmax = world.computeGrainForces();

    // Apply gravity to all grains
    world.applyAcceleration(createVector(0, -10));

    // Integrate motion (both grains and walls)
    world.grainsTimestep();
    world.wallsTimestep(dt);

    // Update neighbor list periodically
    if (istep % nStepsUpsdate == 0) {
        world.updateNeighborList(dmax);
    }

    // increment step counter
    istep++;

    // DISPLAY STUFF
    world.drawParticles(scale, 5 * vRand);
    world.drawWalls(scale);
    counts = world.drawContacts(scale, fnmax, rmin);

    // compute and store data for plots using World data
    solidFrac.push(Vsolid / world.getVolume());
    coordNumber.push(2 * counts / ng);
    // draw graphs
    graph(550, 250, 500, 200, coordNumber, 'Coordination Number');
    graph(550, 500, 590, 200, solidFrac, 'Solid Fraction');

    // display step count
    fill(0);
    text('Step: ' + istep, 900, 20);



}
