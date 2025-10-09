///////////////////////////////////////////////////////////////////////////
// Simple 2D DEM code in JavaScript
// Improved version with better structure and comments using World class
// Grain assembly with 4 walls. 
// Packing, isotropic compression, biaxial compression.
// 2025-01-09 OTo
///////////////////////////////////////////////////////////////////////////

// World object - contains all simulation data
let world;

// Drawing parameters
let colorBy = 'pressure'; // 'velocity' or 'pressure'
let maxParam = 5; // max value for color scale
let scale; // scale for drawing (calculated automatically)

// system variables
let ngw = 5; // number of grains in width
let ngh = 5; // number of grains in height
let rmin = 0.3E-3; // minimum radius
let rmax = 1.0E-3; // maximum radius

let ng = ngw * ngh; // total number of grains
let meanDiameter = rmin + rmax; // mean diameter

// packing
let targetSolidFraction = 0.7;
let wallVelocity = 0.1;

// target confining pressure
let targetPressure = 100.0;  // target confining pressure

// biaxial loading
let targetAxialStrain = 0.1; // target axial strain
let axialVelocity = 0.01; // axial wall velocity 
let initialCellHeight = 0; // to track axial strain

// non-dimensional parameters
let kappa = 1000
let viscoRate = 0.2; // typical for granular materials (0.1-0.3)
let inertiaNumber = 1.0E-3; // currently not used!!
let damping = 0.1; // global damping coefficient

// Grain parameters
let kn = kappa * targetPressure // contact stiffness
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
let nStepsUpsdate = 5; // number of steps for neighbor search update

// Data and plot management
let dataManager;
let plotManager;

// Reference dimensions for strain calculation
let referenceHeight = null;
let referenceWidth = null;

// Values used to control simulation phases
let istep = 0; // step counter
let currentPhase = 0; // simulation phase
let currentSolidFraction;
let currentAxialStrain = 0; // current axial strain
let currentPressure = 0.0;  // current pressure on walls
let currentKineticEnergy; // current kinetic energy
let simulationComplete = false; // flag to indicate simulation is finished

// other parameter declarations
let fnmax = 0.0; // maximum force (display purpose)

function createGrains() {
    // Create grains with random positions and velocities on a grid without overlap
    let grains = [];
    let vRand = rmin / (200 * dt);

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

function createWalls(xmin, xmax, ymin, ymax) {
    let walls = [];
    walls.push(new Wall(createVector(xmin, ymin), createVector(1, 0), Infinity, knWall)); // left
    walls.push(new Wall(createVector(xmax, ymin), createVector(-1, 0), mWall, knWall)); // right
    walls.push(new Wall(createVector(xmin, ymin), createVector(0, 1), Infinity, knWall)); // bottom
    walls.push(new Wall(createVector(xmin, ymax), createVector(0, -1), mWall, knWall)); // top
    return walls;
}


// P5.js setup and draw functions
function setup() {
    createCanvas(1400, 800);
    background(240);

    dataManager = new DataManager();
    plotManager = new PlotManager(dataManager);

    // Calculate world dimensions
    let worldWidth = ngw * 2 * rmax;   // Physical world width in meters
    let worldHeight = ngh * 2 * rmax;  // Physical world height in meters

    // Calculate scale to fit within display bounds
    let maxDisplayWidth = 0.5 * width;   // 500 pixels
    let maxDisplayHeight = 0.9 * height; // 450 pixels

    let scaleX = maxDisplayWidth / worldWidth;   // pixels per meter (width constraint)
    let scaleY = maxDisplayHeight / worldHeight;  // pixels per meter (height constraint)
    scale = min(scaleX, scaleY);  // Use smaller scale to fit both dimensions

    console.log(`World size: ${worldWidth.toExponential(2)}m x ${worldHeight.toExponential(2)}m`);
    console.log(`Display area: ${maxDisplayWidth}px x ${maxDisplayHeight}px`);
    console.log(`Auto-calculated scale: ${scale.toExponential(2)} pixels/meter`);

    // Create grains
    let grains = createGrains();

    // Create walls
    let walls = createWalls(0, ngw * 2 * rmax, 0, ngh * 2 * rmax);

    // Create World with all components
    world = new World(grains, walls, dt, kn, kt, mu, viscoRate, damping);

    // Initialize neighbor list
    world.updateNeighborList(dmax);
}

function draw() {
    background(0);

    // Reset forces first
    world.resetForces();

    // Compute all contact forces (grain-grain and grain-wall)
    world.computeWallForces();
    fnmax = world.computeGrainForces();

    // Calculate current kinetic energy
    currentKineticEnergy = world.calculateKineticEnergy();
    currentSolidFraction = world.calculateSolidFraction();

    // Apply gravity to all grains
    //world.applyAcceleration(createVector(0, -10));

    // Apply boundary conditions...start with packing to a target solid fraction by moving walls inward
    if (currentSolidFraction < targetSolidFraction && currentPhase == 0) {
        // set zero friction during packing
        for (let grain of world.grains) {
            grain.mu = 0.0;
            grain.viscoRate = 1.0; // no dissipation during packing
            grain.damping = 10; // add some global damping
        }
        // move top and right walls inward
        world.walls[1].moveWall(wallVelocity * dt); // move right wall left
        world.walls[3].moveWall(wallVelocity * dt); // move top wall down
    }

    if (currentSolidFraction >= targetSolidFraction && currentPhase == 0) {
        console.log('Packing phase complete at step ' + istep + ' with solid fraction ' + currentSolidFraction);
        currentPhase = 1;
        // Capture reference dimensions for relaxation phase strain calculation
        referenceHeight = world.getCellHeight();
        referenceWidth = world.getCellWidth();
    }

    // Relaxation until kinetic energy is low

    if (currentPhase == 1 && currentKineticEnergy > 1E-6) {
        // do nothing, just let it relax

    } else if (currentPhase == 1 && currentKineticEnergy <= 1E-6) {
        console.log('Relaxation phase complete at step ' + istep + ' with KE ' + currentKineticEnergy.toExponential(3));
        currentPhase = 2;
        // Capture reference dimensions for isotropic compression phase strain calculation
        referenceHeight = world.getCellHeight();
        referenceWidth = world.getCellWidth();
        // set friction back to target value
        for (let grain of world.grains) {
            grain.mu = mu;
            grain.viscoRate = viscoRate; // restore dissipation
            grain.damping = damping; // restore global damping
        }
    }

    // Isotropic compression phase

    if (currentPhase == 2) {
        // Ramp up confining pressure to target value
        if (currentPressure < targetPressure) {
            currentPressure += targetPressure / 1000; // ramp up over 5000 steps
        }
        let forceRightWall = currentPressure * world.getCellHeight();
        let forceTopWall = currentPressure * world.getCellWidth();
        world.wallTimestep(1, forceRightWall); // right wall
        world.wallTimestep(3, forceTopWall); // top wall

        // Transition to Phase 3 when target pressure is reached AND kinetic energy is low
        if (currentPressure >= targetPressure && currentKineticEnergy <= 1E-6) {
            console.log('Isotropic compression complete at step ' + istep + ' with pressure ' + currentPressure.toFixed(3) + ' and KE ' + currentKineticEnergy.toExponential(3));
            currentPhase = 3;
            // Capture reference dimensions for biaxial compression phase strain calculation
            referenceHeight = world.getCellHeight();
            referenceWidth = world.getCellWidth();

        } else if (currentPressure >= targetPressure && currentKineticEnergy > 1E-6) {
            // Pressure reached but still need to wait for equilibrium
            console.log('Pressure reached, waiting for equilibrium. Current KE: ' + currentKineticEnergy.toExponential(3) + ', Target: ' + 1E-6.toExponential(3));
        }
    }

    if (currentPhase == 3) {
        // Biaxial compression: move top wall down while maintaining pressure on right wall

        // Calculate current axial strain using hybrid approach
        let strainData = world.calculateStrainFromReference(referenceHeight, referenceWidth);
        currentAxialStrain = Math.abs(strainData.strainY); // Axial strain (vertical compression)

        if (currentAxialStrain < targetAxialStrain && !simulationComplete) {
            // Continue axial compression: move top wall down
            world.walls[3].moveWall(axialVelocity * dt); // move top wall down

            // Maintain confining pressure on right wall
            let cellHeight = world.getCellHeight();
            let forceRightWall = targetPressure * cellHeight;
            world.wallTimestep(1, forceRightWall); // right wall pressure control
        } else if (currentAxialStrain >= targetAxialStrain && !simulationComplete) {
            // Target strain reached - complete simulation
            simulationComplete = true;
            console.log('Simulation complete! Target axial strain reached at step ' + istep);
            console.log('Final axial strain: ' + (currentAxialStrain * 100).toFixed(2) + '%');

            // Auto-save data when simulation completes
            let timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, -5);
            dataManager.exportToCSV(`output/simulation_complete_${timestamp}.csv`);
            console.log('Data automatically exported to CSV file');
        }

        // Data collection handled by dataManager
    }

    // Only run physics if simulation is not complete
    if (!simulationComplete) {
        world.grainsTimestep();

        // Update neighbor list periodically
        if (istep % nStepsUpsdate == 0) {
            world.updateNeighborList(dmax);
        }

        // increment step counter
        istep++;
    }

    // DISPLAY STUFF
    world.drawParticles(scale, colorBy, maxParam);
    world.drawWalls(scale);
    world.drawContacts(scale, fnmax, rmin);

    // Collect and display data
    dataManager.collect(world, istep, currentPhase, referenceHeight, referenceWidth);
    plotManager.drawAll(currentPhase);

    // Display simulation information
    fill(255); // white text
    textSize(16);

    // Step number
    text('Step: ' + istep, 20, 30);

    // Phase information
    let phaseText = '';
    if (currentPhase == 0) {
        phaseText = 'Phase 0: Packing to Target Solid Fraction';
    } else if (currentPhase == 1) {
        phaseText = 'Phase 1: Relaxation';
    } else if (currentPhase == 2) {
        phaseText = 'Phase 2: Isotropic Compression';
    } else if (currentPhase == 3) {
        phaseText = simulationComplete ? 'Phase 3: Biaxial Compression - COMPLETE' : 'Phase 3: Biaxial Compression';
    }
    text(phaseText, 20, 55);

    // Show completion status
    if (simulationComplete) {
        fill(0, 255, 0); // green text for completion
        text('SIMULATION COMPLETE - Data Saved', 20, 180);
        fill(255); // back to white
    }

    // Phase-specific information
    if (currentPhase == 2) {
        text('Applied Pressure: ' + currentPressure.toFixed(3) + ' Pa', 20, 80);
        text('Target Pressure: ' + targetPressure.toFixed(1) + ' Pa', 20, 105);
        text('Solid Fraction: ' + currentSolidFraction.toFixed(3), 20, 130);
    } else if (currentPhase == 3) {
        text('Confining Pressure: ' + targetPressure.toFixed(1) + ' Pa', 20, 80);
        text('Axial Strain: ' + (currentAxialStrain * 100).toFixed(2) + '%', 20, 105);
        text('Target Strain: ' + (targetAxialStrain * 100).toFixed(1) + '%', 20, 130);
        text('Solid Fraction: ' + currentSolidFraction.toFixed(3), 20, 155);
    } else {
        text('Solid Fraction: ' + currentSolidFraction.toFixed(3), 20, 80);
    }
}

function keyPressed() {
    if (key === 's' || key === 'S') {
        let timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, -5);
        dataManager.exportToCSV(`output/manual_save_${timestamp}.csv`);
        console.log('Data exported to CSV file');
    }
}
