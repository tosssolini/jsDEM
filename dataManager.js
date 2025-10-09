class DataManager {
    constructor() {
        // Time series data
        this.step = [];
        this.phase = [];

        //Scalar fabric terms
        this.coordNumber = [];
        this.solidFraction = [];
        this.voidRatio = [];

        // Energy terms
        this.kineticEnergy = [];
        this.strainEnergy = [];
        this.externalWork = [];
        this.totalEnergy = [];

        // Stress terms (stress from intergranular forces)
        this.stressX = [];
        this.stressY = [];
        this.stressXY = [];
        this.q = [];
        this.p = [];
        this.eta = [];

        // Boundary positions
        this.boundaryPosTop = [];
        this.boundaryPosBottom = [];
        this.boundaryPosLeft = [];
        this.boundaryPosRight = [];

        // Strain 
        this.strainX = [];
        this.strainY = [];
        this.volumetricStrain = [];

        // Boundary forces (not implemented yet)
        this.boundaryForceTop = [];
        this.boundaryForceBottom = [];
        this.boundaryForceLeft = [];
        this.boundaryForceRight = [];
    }

    collect(world, step, phase, referenceHeight, referenceWidth) {

        this.step.push(step);
        this.phase.push(phase);

        // Void ratio and solid fraction
        this.coordNumber.push(world.calculateCoordinationNumber());
        this.solidFraction.push(world.calculateSolidFraction());
        this.voidRatio.push(world.calculateVoidRatio());

        // Energy (only kinetic energy for now)
        this.kineticEnergy.push(world.calculateKineticEnergy());

        // Stress (from intergranular forces)
        let stress = world.calculateVoigtStress();
        this.stressX.push(stress[0]);
        this.stressY.push(stress[1]);
        this.stressXY.push(stress[2]);
        this.q.push((stress[1] - stress[0]) / 2);
        this.p.push((stress[0] + stress[1]) / 2);
        this.eta.push((stress[1] - stress[0]) / (stress[0] + stress[1]));

        // Boundary positions
        this.boundaryPosTop.push(world.getTopWallPosition());
        this.boundaryPosBottom.push(world.getBottomWallPosition());
        this.boundaryPosLeft.push(world.getLeftWallPosition());
        this.boundaryPosRight.push(world.getRightWallPosition());

        // Strain calculation using hybrid approach
        let strain = world.calculateStrainFromReference(referenceHeight, referenceWidth);
        this.strainX.push(strain.strainX);
        this.strainY.push(strain.strainY);
        this.volumetricStrain.push(strain.volumetricStrain);

    }

    reset() {
        this.step = [];
        this.phase = [];
        this.coordNumber = [];
        this.solidFraction = [];
        this.voidRatio = [];
        this.kineticEnergy = [];
        this.strainEnergy = [];
        this.externalWork = [];
        this.totalEnergy = [];
        this.stressX = [];
        this.stressY = [];
        this.stressXY = [];
        this.q = [];
        this.p = [];
        this.eta = [];
        this.strainX = [];
        this.strainY = [];
        this.volumetricStrain = [];
        this.boundaryPosTop = [];
        this.boundaryPosBottom = [];
        this.boundaryPosLeft = [];
        this.boundaryPosRight = [];
        this.boundaryForceTop = [];
        this.boundaryForceBottom = [];
        this.boundaryForceLeft = [];
        this.boundaryForceRight = [];
    }

    exportToCSV(filePath = 'output/simulation_data.csv') {
        const headers = ['step', 'phase', 'coordNumber', 'solidFraction', 'voidRatio',
            'kineticEnergy', 'strainEnergy', 'externalWork', 'totalEnergy',
            'stressX', 'stressY', 'stressXY', 'q', 'p', 'eta',
            'strainX', 'strainY', 'volumetricStrain',
            'boundaryPosTop', 'boundaryPosBottom', 'boundaryPosLeft', 'boundaryPosRight',
            'boundaryForceTop', 'boundaryForceBottom', 'boundaryForceLeft', 'boundaryForceRight'];

        let csv = headers.join(',') + '\n';

        const dataLength = this.step.length;

        for (let i = 0; i < dataLength; i++) {
            const row = [
                this.step[i] || '',
                this.phase[i] || '',
                this.coordNumber[i] || '',
                this.solidFraction[i] || '',
                this.voidRatio[i] || '',
                this.kineticEnergy[i] || '',
                this.strainEnergy[i] || '',
                this.externalWork[i] || '',
                this.totalEnergy[i] || '',
                this.stressX[i] || '',
                this.stressY[i] || '',
                this.stressXY[i] || '',
                this.q[i] || '',
                this.p[i] || '',
                this.eta[i] || '',
                this.strainX[i] || '',
                this.strainY[i] || '',
                this.volumetricStrain[i] || '',
                this.boundaryPosTop[i] || '',
                this.boundaryPosBottom[i] || '',
                this.boundaryPosLeft[i] || '',
                this.boundaryPosRight[i] || '',
                this.boundaryForceTop[i] || '',
                this.boundaryForceBottom[i] || '',
                this.boundaryForceLeft[i] || '',
                this.boundaryForceRight[i] || ''
            ];
            csv += row.join(',') + '\n';
        }

        // Download file
        const blob = new Blob([csv], { type: 'text/csv' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filePath;
        a.click();
        window.URL.revokeObjectURL(url);
    }

    // Get data arrays starting from Phase 3 (biaxial loading)
    getPhase3Data() {
        const phase3StartIndex = this.phase.findIndex(p => p === 3);
        if (phase3StartIndex === -1) {
            return {
                p: [],
                q: [],
                strainY: [],
                eta: [],
                volumetricStrain: []
            };
        }

        return {
            p: this.p.slice(phase3StartIndex),
            q: this.q.slice(phase3StartIndex),
            strainY: this.strainY.slice(phase3StartIndex),
            eta: this.eta.slice(phase3StartIndex),
            volumetricStrain: this.volumetricStrain.slice(phase3StartIndex)
        };
    }

}

