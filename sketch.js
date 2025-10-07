class Grain {
    constructor(mass, radius, xpos, ypos) {
        this.mass = mass;
        this.inertia = 0.5 * mass * radius * radius;
        this.radius = radius;

        this.x = xpos;
        this.y = ypos;
        this.rot = 0.0;
        this.vx = 0.0;
        this.vy = 0.0;
        this.vrot = 0.0;
        this.ax = 0.0;
        this.ay = 0.0;
        this.arot = 0.0;

        // resultant force and torque
        this.fx = 0.0;
        this.fy = 0.0;
        this.frot = 0.0;
    }

    draw() {
        // grat level of velocity
        let vel = Math.sqrt(this.vx * this.vx + this.vy * this.vy);
        let M = max(vRand, 0.8 * abs(topv));
        let colorVel = map(vel, 0.0, M, 255, 0);
        fill(colorVel);

        push();
        translate(this.x * scale, height - (this.y * scale));
        rotate(-this.rot);
        stroke(0);
        let r = this.radius * scale;
        circle(0, 0, 2 * r);
        line(0, 0, r, 0);
        pop();
    }
}

class Neighbor {
    constructor(i, j) {
        this.i = i;
        this.j = j;
        this.touch = false;
        this.fn = 0.0;
        this.ft = 0.0;
    }
}

// system variables
let grain = [];
let ngw = 5; // number of grains in width
let ngh = 5; // number of grains in height
let rmin = 0.5E-3; // minimum radius
let rmax = 1.0E-3; // maximum radius
let density = 2700; // density
let scale = 45000; // scale for drawing
let xmax = 0.0;
let ymax = 0.0;
let ymax0 = 0.0;
let vRand;

// precomputed variables
let ng = ngw * ngh; // total number of grains
let meanDiameter = rmin + rmax; // mean diameter
let meanMass = Math.PI * meanDiameter * meanDiameter * 0.25 * density;

// plot variables
let coordNumber = [];
let solidFrac = [];
let Vsolid = 0.0;
let rwallpos = [];
let twallpos = [];
let rwallf = [];
let twallf = [];

// kinematic of the walls
let topf = 0.0;
let topa = 0.0;
let topv = 0.0;
let rightf = 0.0;
let righta = 0.0;
let rightv = 0.0;

// neighbor array
let neighbor = [];
let dmax = 0.95 * rmin; // maximum distance for neighbor search
let nStepsUpsdate = 25; // number of steps for neighbor search update

// processing
let fnmax = 0.0; // maximum force (display purpose)
let nc = 0; // number of contacts

// non-dimentional parameters
let kappa = 1000
let viscoRate = 0.98;
let InertialNumber = 1.0E-2;

// loading
let stepBiaxStarts = 4000;
let wallMass = 0.0;
let pressure = 1.0;  // target confining pressure
let axialVelocity = (ngh * meanDiameter) * InertialNumber * Math.sqrt(pressure / meanMass);

// pressure control parameters
let pressureController = 0.001; // pressure control gain (adjust for stability)
let currentPressureRight = 0.0;  // current pressure on right wall
let currentPressureTop = 0.0;    // current pressure on top wall
console.log("axialVelocity: " + axialVelocity);

// force law parameters
let kn = kappa * pressure // contact stiffness
let kt = kn;
let mu = 0.1; // friction coefficient
let fcohesion = 0.0; // cohesion force

// time flow
let dt = Math.PI * Math.sqrt(meanMass / kn) / 20 // time step
console.log("dt: " + dt);
let istep = 0; // step counter

function putOnGrid() {
    Vsolid = 0.0;
    for (let iy = 0; iy < ngh; iy++) {
        for (let ix = 0; ix < ngw; ix++) {
            let x = rmax + ix * 2 * rmax;
            let y = rmax + iy * 2 * rmax;
            let radius = random(rmin, rmax);
            Vsolid += Math.PI * radius * radius;
            let mass = Math.PI * radius * radius * density;
            let i = ix + iy * ngw;
            grain.push(new Grain(mass, radius, x, y));
        }
    }

    // random velocity in any direction, and compute initial kinetic energy
    vRand = rmin / (200 * dt); // 5 * axialVelocity
    for (let i = 0; i < grain.length; i++) {
        grain[i].vx = random(-vRand, vRand);
        grain[i].vy = random(-vRand, vRand);
    }

    xmax = ngw * 2 * rmax;
    ymax = ngh * 2 * rmax;
    let dt_critical = Math.PI * Math.sqrt(Math.PI * rmin * rmin * 1.0 * density / kn);
    dt = dt_critical / 10;
    //console.log(dt);

    wallMass = ngw * Math.PI * rmax * rmax * density * 1000;
}

function computeForces(k) {
    let i = neighbor[k].i;
    let j = neighbor[k].j;

    // branch from j to i
    let bx = grain[i].x - grain[j].x;
    let by = grain[i].y - grain[j].y;
    let b = Math.sqrt(bx * bx + by * by);

    // distance between the grains
    let dn = b - (grain[i].radius + grain[j].radius);

    if (dn < 0.0) {
        neighbor[k].touch = true;
        let nx = bx / b;
        let ny = by / b;
        let tx = -ny;
        let ty = nx;

        let vrelx = grain[i].vx - grain[j].vx + grain[i].radius * grain[i].vrot * ny + grain[j].radius * grain[j].vrot * ny;
        let vrely = grain[i].vy - grain[j].vy - grain[i].radius * grain[i].vrot * nx - grain[j].radius * grain[j].vrot * nx;
        let vn = vrelx * nx + vrely * ny;
        let vt = vrelx * tx + vrely * ty;

        let meff = (grain[i].mass * grain[j].mass) / (grain[i].mass + grain[j].mass);
        let visco = viscoRate * 2.0 * Math.sqrt(kn * meff);
        let fn = -kn * dn - visco * vn + fcohesion;
        neighbor[k].fn = fn;

        // let ft = -visco * vt;
        let ft = neighbor[k].ft - kt * vt * dt;
        let threshold = mu * abs(fn);
        if (ft > threshold) {
            ft = threshold;
        } else if (ft < -threshold) {
            ft = -threshold;
        }
        neighbor[k].ft = ft;

        // (fx,fy) = projection of n and t forces into the global frame
        let fx = fn * nx + ft * tx;
        let fy = fn * ny + ft * ty;

        grain[i].fx += fx;
        grain[i].fy += fy;
        grain[i].frot -= -ft * (grain[i].radius + 0.5 * dn);
        grain[i].p += fn;

        grain[j].fx -= fx;
        grain[j].fy -= fy;
        grain[j].frot -= ft * (grain[j].radius + 0.5 * dn);
        grain[j].p += fn;
    } else {
        neighbor[k].touch = false;
        neighbor[k].fn = 0.0;
        neighbor[k].ft = 0.0;
    }
}

function computeForceWallLeft(i) {
    let dn = grain[i].x - grain[i].radius;
    if (dn < 0.0) {
        let vn = grain[i].vx;
        let visco = viscoRate * 2.0 * Math.sqrt(kn * grain[i].mass);
        let fn = -kn * dn - visco * vn;
        grain[i].fx += fn;
        grain[i].p += fn;
    }
}

function computeForceWallRight(i) {
    let dn = (xmax - grain[i].x) - grain[i].radius;
    if (dn < 0.0) {
        let vn = -grain[i].vx - rightv;
        let visco = viscoRate * 2.0 * Math.sqrt(kn * grain[i].mass);
        let fn = -kn * dn - visco * vn;
        grain[i].fx -= fn;
        grain[i].p += fn;
        rightf += fn;
    }
}

function computeForceWallBottom(i) {
    let dn = grain[i].y - grain[i].radius;
    if (dn < 0.0) {
        let vn = grain[i].vy;
        let visco = viscoRate * 2.0 * Math.sqrt(kn * grain[i].mass);
        let fn = -kn * dn - visco * vn;
        grain[i].fy += fn;
        grain[i].p += fn;
    }
}

function computeForceWallTop(i) {
    let dn = (ymax - grain[i].y) - grain[i].radius;
    if (dn < 0.0) {
        let vn = -grain[i].vy - topv;
        let visco = viscoRate * 2.0 * Math.sqrt(kn * grain[i].mass);
        let fn = -kn * dn - visco * vn;
        grain[i].fy -= fn;
        grain[i].p += fn;
        topf += fn;
    }
}

function updateNeighborList() {
    // store the current neighbor list
    let tmp = [];
    for (let k = 0; k < neighbor.length; k++) {
        tmp.push(neighbor[k]);
    }
    // clear the current neighbor list
    neighbor = [];
    neighbor.length = 0;

    // rebuild the nieghbor list 
    for (let i = 0; i < grain.length; i++) {
        for (let j = i + 1; j < grain.length; j++) {

            // branch from j to i
            let bx = grain[i].x - grain[j].x;
            let by = grain[i].y - grain[j].y;
            let b = sqrt(bx * bx + by * by);

            // distance between the grains
            let dn = b - (grain[i].radius + grain[j].radius);
            if (dn <= dmax) {
                neighbor.push(new Neighbor(i, j));
            }
        }
    }

    // retrieve the previous pairs (that held fn and ft values)
    // This is important because ft is updated incrementally
    if (tmp.length > 0) {
        let kprev = 0;
        let nk = neighbor.length;
        let nkprev = tmp.length;
        for (let k = 0; k < nk; k++) {
            while (kprev < nkprev && tmp[kprev].i < neighbor[k].i && tmp[kprev].j < neighbor[k].j) {
                ++kprev;
            }
            if (kprev == nkprev) break;

            while (kprev < nkprev && tmp[kprev].i == neighbor[k].i && tmp[kprev].j < neighbor[k].j) {
                ++kprev;
            }
            if (kprev == nkprev) break;

            if (tmp[kprev].i == neighbor[k].i && tmp[kprev].j == neighbor[k].j) {
                neighbor[k].fn = tmp[kprev].fn;
                neighbor[k].ft = tmp[kprev].ft;
                neighbor[k].touch = tmp[kprev].touch;
            }
        }
    }
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
    putOnGrid();
    updateNeighborList();
}

function draw() {
    background(240);

    // integration for particles (Euler scheme)
    for (let i = 0; i < grain.length; i++) {
        //update position
        grain[i].x += grain[i].vx * dt
        grain[i].y += grain[i].vy * dt
        grain[i].rot += grain[i].vrot * dt

        //update velocity
        grain[i].vx += grain[i].ax * dt
        grain[i].vy += grain[i].ay * dt
        grain[i].vrot += grain[i].arot * dt

        //artificial and brutal dissipation of energy
        grain[i].vx *= 0.9999;
        grain[i].vy *= 0.9999;
        grain[i].vrot *= 0.9999;

        // reset resultant forces
        grain[i].fx = 0.0;
        grain[i].fy = 0.0;
        grain[i].frot = 0.0;
        grain[i].p = 0.0;
    }

    // Reset wall forces FIRST
    topf = 0.0;
    rightf = 0.0;

    // Compute wall forces
    for (let i = 0; i < grain.length; i++) {
        computeForceWallLeft(i);
        computeForceWallRight(i);
        computeForceWallTop(i);
        computeForceWallBottom(i);
    }

    // Compute inter-particle forces
    fnmax = 0.0;
    for (let k = 0; k < neighbor.length; k++) {
        computeForces(k);
        if (neighbor[k].fn > fnmax) {
            fnmax = neighbor[k].fn;
        }
    }

    // Calculate current boundary pressures AFTER computing forces
    let rightWallArea = ymax;
    let topWallArea = xmax;
    currentPressureRight = Math.abs(rightf) / rightWallArea;
    currentPressureTop = Math.abs(topf) / topWallArea;

    let speed = 1
    // Pressure control for walls
    if (istep < stepBiaxStarts) {
        // ISOTROPIC COMPRESSION - both walls compress to target pressure

        // Right wall pressure control (increase gain significantly)
        let rightPressureError = pressure - currentPressureRight;
        rightv = -speed * rightPressureError * meanDiameter; // 100x faster!

        // Top wall pressure control  
        let topPressureError = pressure - currentPressureTop;
        topv = -speed * topPressureError * meanDiameter; // 100x faster!

    } else {
        // BIAXIAL COMPRESSION - maintain pressure on right wall, compress top wall

        // Right wall: maintain confining pressure
        let rightPressureError = pressure - currentPressureRight;
        rightv = -speed * rightPressureError * meanDiameter;

        // Top wall: constant displacement rate (deviatoric loading)
        topv = -axialVelocity;
    }

    // Integrate wall positions
    xmax += rightv * dt;
    ymax += topv * dt;

    // Prevent walls from moving outward beyond initial positions
    if (xmax > ngw * 2 * rmax) xmax = ngw * 2 * rmax;
    if (ymax > ngh * 2 * rmax) ymax = ngh * 2 * rmax;

    // Update accelerations on grains
    for (let i = 0; i < grain.length; i++) {
        grain[i].ax = grain[i].fx / grain[i].mass;
        grain[i].ay = grain[i].fy / grain[i].mass;
        grain[i].arot = grain[i].frot / grain[i].inertia;
    }

    istep++; // increment step counter
    // Update neighbor list periodically
    if (istep % nStepsUpsdate == 0) {
        updateNeighborList();
    }

    // DISPLAY STUFF

    // Draw top and right wall
    stroke(255, 0, 0);
    line(0, height - (ymax * scale), xmax * scale, height - (ymax * scale));
    line(xmax * scale, height - (ymax * scale), xmax * scale, height);

    // Draw particles
    for (let i = 0; i < grain.length; i++) {
        grain[i].draw();
    }

    // neighbours
    stroke(255, 100, 0);
    counts = 0;
    for (let k = 0; k < neighbor.length; k++) {
        if (neighbor[k].touch) {
            let w = map(neighbor[k].fn, 0, fnmax, 0, rmin * scale);
            strokeWeight(w);
            let xi = grain[neighbor[k].i].x * scale;
            let yi = height - (grain[neighbor[k].i].y * scale);
            let xj = grain[neighbor[k].j].x * scale;
            let yj = height - (grain[neighbor[k].j].y * scale);
            line(xi, yi, xj, yj);
            counts++;
        }
    }
    strokeWeight(1);

    // compute and store data for plots
    solidFrac.push(Vsolid / (xmax * ymax));
    coordNumber.push(counts / ng);
    rwallpos.push(xmax);
    twallpos.push(ymax);
    rwallf.push(rightf);
    twallf.push(topf);
    // draw graphs
    graph(550, 250, 500, 200, twallpos, 'Top position');
    graph(550, 500, 590, 200, rwallpos, 'Right position');

    // display step count
    fill(0);
    text('Step: ' + istep, 900, 20);
    text('Top pressure: ' + round(currentPressureTop, 3) + ' (Pa)', 700, 20);
    text('Right pressure: ' + round(currentPressureRight, 3) + ' (Pa)', 700, 40);



}
