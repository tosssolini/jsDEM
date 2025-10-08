class World {
    constructor(grains, walls, dt, kn, kt, mu, viscoRate, damping) {
        this.grains = grains;
        this.leftWall = walls[0];
        this.rightWall = walls[1];
        this.bottomWall = walls[2];
        this.topWall = walls[3];
        this.walls = [walls[0], walls[1], walls[2], walls[3]];
        this.dt = dt;
        this.kn = kn;
        this.kt = kt;
        this.mu = mu;
        this.viscoRate = viscoRate;
        this.damping = damping;
        this.neighbors = [];  // Move from global variable

    }

    applyAcceleration(acceleration) {
        // Apply acceleration to all grains as body force
        // acceleration is a p5.Vector
        for (let grain of this.grains) {
            grain.force.add(p5.Vector.mult(acceleration, grain.mass));
        }
    }

    grainsTimestep() {
        // Update positions and velocities of all grains
        for (let grain of this.grains) {
            grain.updatePosition(this.dt, this.damping);
        }
    }

    // Update movable walls (replaces wallsTimeStep function)
    wallsTimestep() {
        this.walls.forEach(wall => {
            if (wall.mass < Infinity) {  // Only update movable walls
                wall.updatePosition(this.dt);
            }
        });
    }


    resetForces() {
        // Reset forces and pressures on grains and walls
        for (let grain of this.grains) {
            grain.resetForce();
        }
        // Reset all walls
        for (let wall of this.walls) {
            wall.resetForce();
        }
    }

    updateNeighborList(dmax) {
        // store the current neighbor list
        let tmp = [];
        for (let k = 0; k < this.neighbors.length; k++) {
            tmp.push(this.neighbors[k]);
        }
        // clear the current neighbor list
        this.neighbors = [];
        this.neighbors.length = 0;

        // rebuild the neighbor list 
        for (let i = 0; i < this.grains.length; i++) {
            for (let j = i + 1; j < this.grains.length; j++) {

                // branch vector from j to i
                let branch = p5.Vector.sub(this.grains[i].pos, this.grains[j].pos);
                let b = branch.mag();

                // distance between the grains
                let dn = b - (this.grains[i].radius + this.grains[j].radius);
                if (dn <= dmax) {
                    this.neighbors.push(new Neighbor(i, j));
                }
            }
        }

        // retrieve the previous pairs (that held fn and ft values)
        // This is important because ft is updated incrementally
        if (tmp.length > 0) {
            let kprev = 0;
            let nk = this.neighbors.length;
            let nkprev = tmp.length;
            for (let k = 0; k < nk; k++) {
                // FIXED: Proper lexicographic ordering - skip old neighbors that come before current new neighbor
                while (kprev < nkprev && (tmp[kprev].i < this.neighbors[k].i ||
                    (tmp[kprev].i == this.neighbors[k].i && tmp[kprev].j < this.neighbors[k].j))) {
                    ++kprev;
                }
                if (kprev == nkprev) break;

                // Check for exact match and restore force history
                if (kprev < nkprev && tmp[kprev].i == this.neighbors[k].i && tmp[kprev].j == this.neighbors[k].j) {
                    this.neighbors[k].fn = tmp[kprev].fn;
                    this.neighbors[k].ft = tmp[kprev].ft;
                    this.neighbors[k].touch = tmp[kprev].touch;
                }
            }
        }
    }

    computeNeighborForce(k) {
        let i = this.neighbors[k].i;
        let j = this.neighbors[k].j;


        kn = this.kn
        kt = this.kt;
        mu = this.mu;
        viscoRate = this.viscoRate;
        dt = this.dt;

        // branch vector from j to i
        let branch = p5.Vector.sub(this.grains[i].pos, this.grains[j].pos);
        let b = branch.mag();

        // distance between the grains
        let dn = b - (this.grains[i].radius + this.grains[j].radius);

        if (dn < 0.0) {
            this.neighbors[k].touch = true;

            // normal and tangent unit vectors
            let normal = p5.Vector.div(branch, b);  // normalize branch vector
            let tangent = createVector(-normal.y, normal.x);  // perpendicular to normal

            // relative velocity calculation
            let vrel = p5.Vector.sub(this.grains[i].vel, this.grains[j].vel);

            let vn = p5.Vector.dot(vrel, normal);  // normal component of relative velocity
            let vt = p5.Vector.dot(vrel, tangent); // tangential component of relative velocity

            let meff = (this.grains[i].mass * this.grains[j].mass) / (this.grains[i].mass + this.grains[j].mass);
            let visco = viscoRate * 2.0 * Math.sqrt(kn * meff);
            let fn = kn * (-dn) + visco * (-vn);  // Repulsive force (positive when in contact)
            this.neighbors[k].fn = fn;

            // Tangential force with both elastic and viscous components
            let ft = this.neighbors[k].ft - kt * vt * dt - visco * vt;
            let threshold = mu * abs(fn);  // fn is now positive, so abs() is still correct
            if (ft > threshold) {
                ft = threshold;
            } else if (ft < -threshold) {
                ft = -threshold;
            }
            this.neighbors[k].ft = ft;

            // force vector = projection of normal and tangential forces
            let forceVec = p5.Vector.mult(normal, fn);  // fn > 0, normal points from j to i
            forceVec.add(p5.Vector.mult(tangent, ft));

            this.grains[i].force.add(forceVec);    // grain i gets pushed away from grain j
            this.grains[j].force.sub(forceVec);    // grain j gets pushed away from grain i
            this.grains[i].p += fn;  // fn is now positive (compressive pressure)
            this.grains[j].p += fn;  // fn is now positive (compressive pressure)
        } else {
            this.neighbors[k].touch = false;
            this.neighbors[k].fn = 0.0;
            this.neighbors[k].ft = 0.0;
        }
    }

    computeGrainForces() {
        // Compute inter-particle forces for all grain pairs
        let fnmax = 0.0;

        for (let k = 0; k < this.neighbors.length; k++) {
            this.computeNeighborForce(k);
            if (this.neighbors[k].fn > fnmax) {
                fnmax = this.neighbors[k].fn;
            }
        }

        return fnmax;  // Return maximum force for visualization
    }

    computeWallForces() {
        // Unified wall force computation - replaces all 4 separate functions!
        viscoRate = this.viscoRate;
        // Compute forces for all grain-wall contacts
        for (let i = 0; i < this.grains.length; i++) {
            for (let wall of this.walls) {
                wall.computeForce(this.grains[i], viscoRate);
            }
        }
    }

    // Get wall forces for analysis/control (replaces rightf, topf variables)
    getRightWallForce() {
        return this.rightWall.force;
    }

    getTopWallForce() {
        return this.topWall.force;
    }

    getVolume() {
        // Calculate current volume of the system
        let width = this.rightWall.position.x - this.leftWall.position.x;
        let height = this.topWall.position.y - this.bottomWall.position.y;
        return width * height;
    }

    getGrainVolume() {
        // Calculate total volume of grains
        let volume = 0;
        for (let grain of this.grains) {
            volume += grain.getVolume();
        }
        return volume;
    }

    getKineticEnergy() {
        // Calculate total kinetic energy of grains
        let KE = 0;
        for (let grain of this.grains) {
            KE += 0.5 * grain.mass * grain.vel.magSq();
        }
        return KE;
    }

    // Draw functions for visualization

    drawContacts(scale, fnmax, rmin) {
        // Draw contact forces between grains
        stroke(255, 100, 0);  // Orange color for contacts
        let contactCount = 0;

        for (let k = 0; k < this.neighbors.length; k++) {
            if (this.neighbors[k].touch) {
                // Line thickness based on force magnitude
                let w = map(this.neighbors[k].fn, 0, fnmax, 0, rmin * scale);
                strokeWeight(w);

                // Get positions of the two grains in contact
                let xi = this.grains[this.neighbors[k].i].pos.x * scale;
                let yi = height - (this.grains[this.neighbors[k].i].pos.y * scale);
                let xj = this.grains[this.neighbors[k].j].pos.x * scale;
                let yj = height - (this.grains[this.neighbors[k].j].pos.y * scale);

                // Draw contact line
                line(xi, yi, xj, yj);
                contactCount++;
            }
        }

        strokeWeight(1);  // Reset stroke weight
        return contactCount;  // Return number of contacts for coordination number
    }
    drawParticles(scale, maxVelocity) {
        // Draw all grains
        for (let grain of this.grains) {
            grain.draw(scale, maxVelocity);
        }
    }

    drawWalls(scale) {
        // Draw the box created by the walls (walls are given by their normal and a position on the line)
        stroke(255);
        strokeWeight(2);
        // Get the corners of the box
        let x1 = this.leftWall.position.x;
        let x2 = this.rightWall.position.x;
        let y1 = this.bottomWall.position.y;
        let y2 = this.topWall.position.y;

        // Draw the box
        noFill();
        rect(x1 * scale, height - (y2 * scale), (x2 - x1) * scale, (y2 - y1) * scale);
    }
}