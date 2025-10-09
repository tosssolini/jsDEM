class PlotManager {
    constructor(dataManager) {
        this.data = dataManager;
    }

    timePlot(x, y, w, h, data, title) {
        stroke(255);
        line(x, y, x, y - h);
        line(x, y, x + w, y);

        if (data.length < 2) return;

        let ymin = Math.min(...data);
        let ymax = Math.max(...data);
        if (ymax === ymin) ymax = ymin + 1;

        let dx = w / (data.length - 1);
        let yscale = h / (ymax - ymin);

        stroke(0, 255, 0);
        for (let i = 0; i < data.length - 1; i++) {
            line(x + dx * i, y - (data[i] - ymin) * yscale,
                x + dx * (i + 1), y - (data[i + 1] - ymin) * yscale);
        }

        noStroke();
        fill(255);
        textSize(14);
        let current = data[data.length - 1] || 0;

        // Use scientific notation for Energy plots
        if (title.toLowerCase().includes('energy')) {
            text(`${title} min: ${ymin.toExponential(2)} max: ${ymax.toExponential(2)} current: ${current.toExponential(2)}`,
                x + 3, y - h + 10);
        } else {
            text(`${title} min: ${ymin.toFixed(3)} max: ${ymax.toFixed(3)} current: ${current.toFixed(3)}`,
                x + 3, y - h + 10);
        }
    }

    scatterPlot(x, y, w, h, xdata, ydata, title) {
        stroke(255);
        line(x, y, x, y - h);
        line(x, y, x + w, y);

        if (xdata.length < 2 || ydata.length < 2) {
            noStroke();
            fill(255);
            textSize(14);
            text(title, x + 3, y - h + 15);
            text('(No data yet)', x + 3, y - h + 35);
            return;
        }

        let xmin = Math.min(...xdata);
        let xmax = Math.max(...xdata);
        let ymin = Math.min(...ydata);
        let ymax = Math.max(...ydata);

        if (xmax === xmin) xmax = xmin + 1;
        if (ymax === ymin) ymax = ymin + 1;

        let xscale = w / (xmax - xmin);
        let yscale = h / (ymax - ymin);

        stroke(0, 255, 0);
        strokeWeight(2);
        for (let i = 0; i < Math.min(xdata.length, ydata.length) - 1; i++) {
            let x1 = x + (xdata[i] - xmin) * xscale;
            let y1 = y - (ydata[i] - ymin) * yscale;
            let x2 = x + (xdata[i + 1] - xmin) * xscale;
            let y2 = y - (ydata[i + 1] - ymin) * yscale;
            line(x1, y1, x2, y2);
        }

        strokeWeight(1);
        noStroke();
        fill(255);
        textSize(14);
        text(title, x + 3, y - h + 15);
        text(`X: ${xmin.toFixed(2)} to ${xmax.toFixed(2)}`, x + 3, y - h + 30);
        text(`Y: ${ymin.toFixed(2)} to ${ymax.toFixed(2)}`, x + 3, y - h + 45);
    }

    drawAll(currentPhase) {
        // Calculate plot dimensions for 2 columns × 3 rows grid on right half of canvas
        // Canvas is 1400x800, right half starts at x=700
        // Align bottom of plots with bottom wall (y=800)
        const startX = 700;
        const plotWidth = 300;   // (1400-700)/2 = 350, use 300 for margins
        const plotHeight = 240;  // Larger plots to fill more space
        const marginX = 25;
        const marginY = 20;

        // Column positions
        const col1X = startX + marginX;           // Left column
        const col2X = col1X + plotWidth + marginX; // Right column

        // Row positions - align bottom plot with bottom wall at y=800
        const bottomWallY = 800;  // Bottom wall canvas position
        const row3Y = bottomWallY; // Bottom row aligns with bottom wall
        const row2Y = row3Y - plotHeight - marginY;   // Middle row  
        const row1Y = row2Y - plotHeight - marginY;   // Top row

        // Left Column: Always visible plots
        this.timePlot(col1X, row1Y, plotWidth, plotHeight, this.data.kineticEnergy, 'Energy');
        this.timePlot(col1X, row2Y, plotWidth, plotHeight, this.data.solidFraction, 'Solid Fraction');
        this.timePlot(col1X, row3Y, plotWidth, plotHeight, this.data.p, 'Boundary Pressures');

        // Right Column: Only visible during biaxial loading (Phase 3)
        if (currentPhase >= 3) {
            const phase3Data = this.data.getPhase3Data();

            this.scatterPlot(col2X, row1Y, plotWidth, plotHeight, phase3Data.p, phase3Data.q, 'p vs q');
            this.scatterPlot(col2X, row2Y, plotWidth, plotHeight, phase3Data.strainY, phase3Data.eta, 'Axial Strain vs η');
            this.scatterPlot(col2X, row3Y, plotWidth, plotHeight, phase3Data.strainY, phase3Data.volumetricStrain, 'Axial vs Vol Strain');
        }
    }
}