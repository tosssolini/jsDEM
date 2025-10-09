class Neighbor {
    constructor(i, j) {
        this.i = i;
        this.j = j;
        this.touch = false;
        this.fn = 0.0;
        this.ft = 0.0;
    }

    draw(grains, scale, fnmax, rmin) {
        if (!this.touch) return; // Only draw if contact exists

        // Line thickness based on force magnitude
        let w = map(this.fn, 0, fnmax, 0, rmin * scale);
        strokeWeight(w);

        // Get positions of the two grains in contact
        let xi = grains[this.i].pos.x * scale;
        let yi = height - (grains[this.i].pos.y * scale);
        let xj = grains[this.j].pos.x * scale;
        let yj = height - (grains[this.j].pos.y * scale);

        // Draw contact line
        stroke(255, 100, 0);  // Orange color for contacts
        line(xi, yi, xj, yj);
    }
}