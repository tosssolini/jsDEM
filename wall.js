class Wall {
    constructor(position, normal, mass, kn) {
        this.position = position.copy();     // A point on the wall (p5.Vector)
        this.normal = normal.copy().normalize(); // Normal vector pointing inward
        this.mass = mass;             // Wall mass (Infinity for fixed walls)
        this.kn = kn;                 // Wall stiffness 
        this.velocity = createVector(0, 0); // Wall velocity
        this.acc = createVector(0, 0); // Wall acceleration
        this.force = 0;              // Total force on wall (scalar, along normal)
    }

    // Distance from grain center to wall (negative = penetration)
    distanceToGrain(grain) {
        let toGrain = p5.Vector.sub(grain.pos, this.position);
        return p5.Vector.dot(toGrain, this.normal) - grain.radius;
    }

    // Compute contact force on a grain
    computeForce(grain, viscoRate) {
        let dn = this.distanceToGrain(grain);
        if (dn < 0.0) {
            // Effective stiffness for grain-wall contact
            let k_eff;
            if (this.kn === Infinity) {
                k_eff = grain.kn;  // Rigid wall case
            } else {
                // Safety check and proper effective stiffness
                if (grain.kn > 0 && this.kn > 0) {
                    k_eff = (grain.kn * this.kn) / (grain.kn + this.kn);  // Springs in series
                } else {
                    k_eff = Math.max(grain.kn || 0, this.kn || 0);  // Fallback
                }
            }

            // Normal velocity component (relative to wall) - CORRECT
            let relVel = p5.Vector.sub(grain.vel, this.velocity);
            let vn = p5.Vector.dot(relVel, this.normal);

            // Proper damping coefficient calculation
            let meff = grain.mass;  // For grain-wall contact
            let visco = viscoRate * 2.0 * Math.sqrt(k_eff * meff);

            // Contact force magnitude - CORRECTED SIGNS
            let fn = k_eff * (-dn) + visco * (-vn);  // Repulsive spring + damping

            // Safety check for finite values
            if (!isFinite(fn)) {
                console.warn("Non-finite wall force:", fn, "k_eff:", k_eff, "dn:", dn, "vn:", vn, "visco:", visco);
                fn = 0;
            }

            // Force vector on grain (normal points inward to domain)
            let forceVec = p5.Vector.mult(this.normal, fn);
            grain.force.add(forceVec);
            grain.p += fn;  // Positive pressure during contact

            // Reaction force on wall (Newton's 3rd law)
            this.force += fn;

            return fn;
        }
        return 0;
    }

    updatePosition(dt) {
        // Euler integration 
        this.position.add(p5.Vector.mult(this.velocity, dt)); //update position
        this.velocity.add(p5.Vector.mult(this.acc, dt)); //update velocity        
        this.velocity.mult(0.9999); //artificial and brutal dissipation of energy
    }

    resetForce() {
        this.force = 0;
    }

}