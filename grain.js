class Grain {
    constructor(position, velocity, radius, density) {
        this.density = density;
        this.mass = Math.PI * radius * radius * density;
        this.radius = radius;
        this.pos = position;
        this.vel = velocity;
        this.acc = createVector(0.0, 0.0);
        this.force = createVector(0.0, 0.0);
        this.p = 0.0; // pressure
    }

    updatePosition(dt, damping) {
        // Semi-implicit integration with damping
        // m×dv/dt = F - γ×m×v
        // v_new = [1/(1 + γdt/2)] × [(1 - γdt/2)×v_old + dt×F/m]  // Always stable
        this.vel = p5.Vector.mult(this.vel, (1.0 - damping * dt / 2.0) / (1.0 + damping * dt / 2.0));
        this.vel.add(p5.Vector.mult(this.force, dt / this.mass / (1.0 + damping * dt / 2.0)));
        this.pos.add(p5.Vector.mult(this.vel, dt));
    }

    resetForce() {
        this.force.set(0.0, 0.0);
        this.p = 0.0;
    }

    getVolume() {
        // actually area in 2D
        return Math.PI * Math.pow(this.radius, 2);
    }

    draw(scale, maxVelocity) {

        let vel = this.vel.mag(); // get magnitude of velocity vector
        let colorVel = map(vel, 0.0, maxVelocity, 255, 0);
        fill(colorVel);
        push();
        translate(this.pos.x * scale, height - (this.pos.y * scale));
        stroke(0);
        let r = this.radius * scale;
        circle(0, 0, 2 * r);
        pop();
    }
}