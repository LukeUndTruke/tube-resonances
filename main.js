createModule().then(Module => {
    const computeFieldCylinder = Module.cwrap("computeFieldCylinder", null, ["number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number"]);
    const computeFieldCone = Module.cwrap("computeFieldCone", null, ["number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number"]);
    const computeImpedanceAbsCylinder = Module.cwrap("computeImpedanceAbsCylinder", "number", ["number", "number", "number", "number", "number"])
    const computeImpedanceAbsCone = Module.cwrap("computeImpedanceAbsCone", "number", ["number", "number", "number", "number", "number", "number"])
    const computeImpedanceSpectrumCylinder = Module.cwrap("computeImpedanceSpectrumCylinder", null, ["number", "number", "number", "number", "number", "number", "number", "number"])
    const computeImpedanceSpectrumCone = Module.cwrap("computeImpedanceSpectrumCone", null, ["number", "number", "number", "number", "number", "number", "number", "number", "number"])

    function maxArray(arr) {
        let m = -Infinity;
        for (let i=0; i<arr.length; i++) {
            if (arr[i] > m) m = arr[i];
        }
        return m;
    }

    function safeRealloc(ptr, size) {
        if (ptr) Module._free(ptr);
        return Module._malloc(size);
    }

    //function for interpolating array to calculate displacement at arbitrary position from velocity array
    function sample(arr, x, L, N) {
        const idx = x/L * (N-1);
        const i = Math.floor(idx);
        const frac = idx-i;
        if (i>=N-1) {return arr[N-1]};
        return arr[i]*(1-frac) + arr[i+1]*frac;
    }

    function samplePhase(arr, x, L, N) {
    // unwrap phases to avoid jumps
        const idx = x/L * (N-1);
        const i = Math.floor(idx);
        const frac = idx - i;
        if (i >= N - 1) return arr[N - 1];

        let p0 = arr[i];
        let p1 = arr[i + 1];
        
        // unwrap
        while (p1 - p0 > Math.PI) p1 -= 2 * Math.PI;
        while (p1 - p0 < -Math.PI) p1 += 2 * Math.PI;
        
        return p0 * (1 - frac) + p1 * frac;
    }


    function displacement(g, x, t) {
        const A = g.type==="cylinder" ? Math.PI*g.R**2 : 2*Math.PI*(g.R1 + x*(g.R2-g.R1)/g.L)**2/(1+g.costh); //area
        const uAbs = sample(g.uAbs, x, g.L, g.N);
        const uArg = samplePhase(g.uArg, x, g.L, g.N);
        return uAbs/(2*Math.PI*f)/A * Math.sin(uArg + 2*Math.PI*f*t);
    }

    let pxPerUnit = 200;
    let timeScale = 0.001;
    let densityParticles = 5000;

    // let scaleFactor = 5; //particle displacement exaggeration factor


    // global frequency
    let f = 200;

    const fMin = 0;
    const fMax = 1500;


    // for array lengths
    const ptsPerUnit = 500;


    function createGeometry(type, params) {

        const g = {
            type, 
            
            L: params.L,
            R1: params.R1,
            R2: params.R2,
            R: params.R,
            Zr: params.Zr,
            Zi: params.Zi,
            
            N: params.N,
            NFreq: null,

            ptrs: {       
                xPtr: Module._malloc(params.N*8),
                pAbsPtr: Module._malloc(params.N*8),
                pArgPtr: Module._malloc(params.N*8),
                pPtr: Module._malloc(params.N*8),
                uAbsPtr: Module._malloc(params.N*8),
                uArgPtr: Module._malloc(params.N*8),
                uPtr: Module._malloc(params.N*8),
                freqPtr: null,
                ZinAbsPtr: null,
            },

            x: null,
            pAbs: null,
            pArg: null,
            p: null,
            uAbs: null,
            uArg: null,
            u: null,

            freqArray: null,
            ZinAbsArray: null,

            computeField: null,
            computeImpedanceAbs: null,
            updateDerived: null,
        }

        g.x = new Float64Array(Module.HEAPF64.buffer, g.ptrs.xPtr, g.N);
        g.pAbs = new Float64Array(Module.HEAPF64.buffer, g.ptrs.pAbsPtr, g.N);
        g.pArg = new Float64Array(Module.HEAPF64.buffer, g.ptrs.pArgPtr, g.N);
        g.p= new Float64Array(Module.HEAPF64.buffer, g.ptrs.pPtr, g.N);
        g.uAbs = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uAbsPtr, g.N);
        g.uArg = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uArgPtr, g.N);
        g.u = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uPtr, g.N);

        g.updateDerived = () => {
            for (let i = 0; i < g.N; i++) {
                g.x[i] = g.L * i/(g.N-1);
            }

            if (g.type === "cone") {
                g.dl = g.R1/(g.R2-g.R1)*g.L; //distance from left end to (virtual) apex
                g.costh = Math.sqrt(g.L**2 - (g.R2-g.R1)**2)/g.L;
                g.theta = Math.acos(g.costh);
            }

            if (g.type=="cylinder") {
                g.rMax = g.R;
            }
            else if (g.type=="cone") {
                g.rMax = Math.max( g.R1, g.R2);
            }

        };


        g.computeField = (f) => {
            if (g.type=="cylinder") {
                computeFieldCylinder(g.L, g.R, f, g.Zr, g.Zi, g.ptrs.xPtr, g.N, g.ptrs.pAbsPtr, g.ptrs.pArgPtr, g.ptrs.uAbsPtr, g.ptrs.uArgPtr);
            }
            else if (g.type=="cone") {
                computeFieldCone(g.L, g.R1, g.R2, f, g.Zr, g.Zi, g.ptrs.xPtr, g.N, g.ptrs.pAbsPtr, g.ptrs.pArgPtr, g.ptrs.uAbsPtr, g.ptrs.uArgPtr);
            }
        }


        g.computeImpedanceAbs = (f) => {
            if (g.type=="cylinder") {
                return computeImpedanceAbsCylinder(g.L, g.R, f, g.Zr, g.Zi);
            }
            else if (g.type=="cone") {
                return computeImpedanceAbsCone(g.L, g.R1, g.R2, f, g.Zr, g.Zi);
            }
        }

        g.updateDerived();
        return g;

    }
    
    const geometries = [
        createGeometry("cylinder", {L: 1.5, R: 0.5, R1: 0.1, R2: 0.5, freq: f, Zr: 1e-5, Zi: 1e-5, N: Math.round(ptsPerUnit*1.5)}),
        createGeometry("cone", {L: 1.5, R: 0.5, R1: 0.1, R2: 0.5, freq: f, Zr: 1e-5, Zi: 1e-5, N: Math.round(ptsPerUnit*1.5)})
    ];


    function adaptiveSpectrumSample(f1, f2, Z1, Z2, g, depth=0) {
        const fmid = 0.5 * (f1 + f2);
        const Zmid = g.computeImpedanceAbs(fmid);

        //linear interpolation in log space
        const Zlin = Math.sqrt(Z1*Z2);

        const error = Math.abs (Math.log(Zmid) - Math.log(Zlin) );

        const threshold = 0.0005;
        const maxDepth = 20;
        const dfMax = 343.0/(2*g.L)/10.0;

        if ((error<threshold && f2-f1<dfMax) || depth>maxDepth ) {
            return [
                {f: f1, Z: Z1},
                {f: f2, Z:Z2}
            ];
        }

        const left = adaptiveSpectrumSample(f1, fmid, Z1, Zmid, g, depth+1);
        const right = adaptiveSpectrumSample(fmid, f2, Zmid, Z2, g, depth+1);

        return left.slice(0,-1).concat(right);
    }

    function computeSpectrumAdaptive(g, fmin, fmax) {
        const Z1 = g.computeImpedanceAbs(fmin);
        const Z2 = g.computeImpedanceAbs(fmax);
        const data = adaptiveSpectrumSample(fmin, fmax, Z1, Z2, g);

        const freq = data.map(p => p.f);
        const Zin = data.map(p => p.Z);

        g.NFreq = freq.length;
        g.ptrs.freqPtr = safeRealloc(g.ptrs.freqPtr, g.NFreq*8); 
        g.ptrs.ZinAbsPtr = safeRealloc(g.ptrs.ZinAbsPtr, g.NFreq*8);

        g.freqArray = new Float64Array(Module.HEAPF64.buffer, g.ptrs.freqPtr, g.NFreq);
        g.ZinAbsArray = new Float64Array(Module.HEAPF64.buffer, g.ptrs.ZinAbsPtr, g.NFreq);

        g.freqArray.set(freq);
        g.ZinAbsArray.set(Zin); 
    }


    //initialize
    let globalpMax = 0;
    let globaluMax = 0;

    geometries.forEach( (g, i) => {
        computeSpectrumAdaptive(g, fMin, fMax);
        g.computeField(f);
        for (let i=0; i<g.N; i++) {
            g.p[i] = g.pAbs[i] * Math.cos(g.pArg[i]);
            g.u[i] = g.uAbs[i] * Math.cos(g.uArg[i]);
        }
        g.pMax = maxArray(g.pAbs)
        g.uMax = maxArray(g.uAbs);
        if (g.pMax > globalpMax) {
            globalpMax = g.pMax;
        }
        if (g.uMax > globaluMax) {
            globaluMax = g.uMax;
        }
    });



    function updateFrequency(f) {
        globalpMax = 0;
        globaluMax = 0;

        geometries.forEach( (g,i) => {
            g.computeField(f);
            g.pMax = maxArray(g.pAbs);
            g.uMax = maxArray(g.uAbs);
            if (g.pMax > globalpMax) {
                globalpMax = g.pMax;
            }
            if (g.uMax > globaluMax) {
                globaluMax = g.uMax;
            }
            Plotly.relayout(`spectrum${i}`, {
            "shapes[0].x0": f,
            "shapes[0].x1": f
            });
        })
    }

    function bindFrequencyControls() {

        //initialize
        const freqSlider = document.getElementById("freqSlider");
        const freqInput = document.getElementById("freqInput");
        freqSlider.value = f;
        freqSlider.min = fMin;
        freqSlider.max = fMax;
        freqInput.value = f;
        freqInput.min = fMin;
        freqInput.max = fMax;

        freqSlider.addEventListener("input", () => {
            f = +freqSlider.value;
            freqInput.value = +freqSlider.value;
            updateFrequency(f);
        })
        freqInput.addEventListener("change", () => {
            f = +freqInput.value;
            freqSlider.value = + freqInput.value;
            updateFrequency(f);
        })
    }
    bindFrequencyControls();
    


    const synchronize = document.getElementById("synchronizeToggle").checked;

    //  Longitudinal Plot
    function initParticles(g) {
        let area;
        if (g.type =="cylinder") {
            area = 2*g.R*g.L;
        }
        else if (g.type=="cone") {
            if (Math.abs(g.R1-g.R2) < 0.00001) area = 2*g.R1*g.L;
            else area = Math.abs(g.theta*((g.L+g.dl)**2-g.dl**2));
        }

        const numParticles = Math.round(densityParticles*area);
        g.particles = [];
        if (g.type=="cylinder") {
            for (let i=0; i<numParticles; i++) {
                const x0 = Math.random()*g.L;
                const y0 = 2.0*(Math.random()-0.5)*g.R;
                g.particles.push({x0:x0, y0:y0, x:x0, y:y0});
            }
        }
        else if (g.type=="cone") {

            if (Math.abs(g.R1-g.R2) < 0.00001) {
                for (let i=0; i<numParticles; i++) {
                    const x0 = Math.random()*g.L;
                    const y0 = 2.0*(Math.random()-0.5)*g.R1;
                    g.particles.push({x0:x0, y0:y0, x:x0, y:y0});
                }
            }
            else {
                for (let i=0; i<numParticles; i++) {
                    // via inversion sampling
                    const u = Math.random();
                    const r0 = g.L/(g.R2-g.R1) * (-g.R1 + Math.sqrt(g.R1*g.R1 + (g.R2*g.R2-g.R1*g.R1)*u));
                    const phi0 = 2*(Math.random()-0.5) * g.theta;
                    g.particles.push({r0:r0, phi0:phi0, r:r0});
                }
            } 
        }
    }
    
    geometries.forEach((g) => initParticles(g));


    function drawParticles(canvas, g) {
        const ctx = canvas.getContext("2d");
        ctx.fillStyle = "black";
        ctx.setTransform(1,0,0,1,0,0);
        ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
        const dpr = window.devicePixelRatio;
        const canvasHeightHalf = ctx.canvas.height / 2;
        
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);    
        const particles = g.particles;

        for (let k = 0; k < particles.length; k++) {
            const p = particles[k];

            const px = leftMargin + p.x * pxPerUnit;
            const py = canvasHeightHalf + (p.y * pxPerUnit);

            ctx.fillRect(px, py, 2, 2);
        }
    }

        
    function allocateFieldArrays(g) {
        // free old memory if it exists
        if (g.ptrs.xPtr) {
            Module._free(g.ptrs.xPtr);
            Module._free(g.ptrs.pAbsPtr);
            Module._free(g.ptrs.pArgPtr);
            Module._free(g.ptrs.pPtr);
            Module._free(g.ptrs.uAbsPtr);
            Module._free(g.ptrs.uArgPtr);
            Module._free(g.ptrs.uPtr);
        }

        if (g.ptrs.freqPtr) Module._free(g.ptrs.freqPtr);
        if (g.ptrs.ZinAbsPtr) Module._free(g.ptrs.ZinAbsPtr);

        // allocate new memory
        g.ptrs.xPtr = Module._malloc(g.N * 8);
        g.ptrs.pAbsPtr = Module._malloc(g.N * 8);
        g.ptrs.pArgPtr = Module._malloc(g.N * 8);
        g.ptrs.pPtr = Module._malloc(g.N * 8);
        g.ptrs.uAbsPtr = Module._malloc(g.N * 8);
        g.ptrs.uArgPtr = Module._malloc(g.N * 8);
        g.ptrs.uPtr = Module._malloc(g.N * 8);
        
        g.ptrs.freqPtr = null;
        g.ptrs.ZinAbsPtr = null;
        g.ptrs.ZinArgPtr = null;

        // recreate typed views
        g.x = new Float64Array(Module.HEAPF64.buffer, g.ptrs.xPtr, g.N);
        g.pAbs = new Float64Array(Module.HEAPF64.buffer, g.ptrs.pAbsPtr, g.N);
        g.pArg = new Float64Array(Module.HEAPF64.buffer, g.ptrs.pArgPtr, g.N);
        g.p = new Float64Array(Module.HEAPF64.buffer, g.ptrs.pPtr, g.N);
        g.uAbs = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uAbsPtr, g.N);
        g.uArg = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uArgPtr, g.N);
        g.u = new Float64Array(Module.HEAPF64.buffer, g.ptrs.uPtr, g.N);
    }

    
    function updateGeometry(g, i) {

        g.N = Math.round(ptsPerUnit * g.L);
        allocateFieldArrays(g);
        g.updateDerived();
        computeSpectrumAdaptive(g, fMin, fMax);
        g.computeField(f);
        g.pDynamic = new Float64Array(g.N);
        g.uDynamic = new Float64Array(g.N);
        initParticles(g);


        globalpMax = 0;
        globaluMax = 0;
        geometries.forEach( (g, i) => {
            g.pMax = maxArray(g.pAbs);
            g.uMax = maxArray(g.uAbs);
            if (g.pMax > globalpMax) {
                globalpMax = g.pMax;
            }
            if (g.uMax > globaluMax) {
                globaluMax = g.uMax;
            }
        }) 
        


        // update canvas geometry
        const dpr = window.devicePixelRatio;
        const canvas = animCanvases[i];
        canvas.style.width = `${g.L * pxPerUnit + leftMargin + rightMargin}px`;
        canvas.style.height = `${2 * g.rMax * pxPerUnit}px`;
        canvas.width = canvas.clientWidth * dpr;
        canvas.height = canvas.clientHeight * dpr;

        const fieldPlots = document.getElementById(`field${i}`);
        fieldPlots.style.width = `${g.L * pxPerUnit + leftMargin + rightMargin}px`;

        // update plots
        Plotly.update(`field${i}`, {
            x: [g.x, g.x],
            y: [g.p, g.u]
        });
        Plotly.relayout(`field${i}`, {
            width: g.L * pxPerUnit + leftMargin + rightMargin,
            'xaxis.range': [0, g.L]
        });
  
        Plotly.update(`spectrum${i}`, {
            x: g.freqArray,
            y: g.ZinAbsArray
        })
        
    }

    function updateGeometryControlState(g, i) {
        const Rslider = document.getElementById(`Rslider${i}`);
        const Rinput  = document.getElementById(`Rinput${i}`);
        const R1slider = document.getElementById(`R1slider${i}`);
        const R1input  = document.getElementById(`R1input${i}`);
        const R2slider = document.getElementById(`R2slider${i}`);
        const R2input  = document.getElementById(`R2input${i}`);

        if (g.type === "cylinder") {
            Rslider.disabled = false;
            Rinput.disabled = false;
            R1slider.disabled = true;
            R1input.disabled = true;
            R2slider.disabled = true;
            R2input.disabled = true;
        }

        else if (g.type === "cone") {
            Rslider.disabled = true;
            Rinput.disabled = true;
            R1slider.disabled = false;
            R1input.disabled = false;
            R2slider.disabled = false; 
            R2input.disabled = false;
        }
    }

    function bindGeometryControls(g, i) {
        const Lslider = document.getElementById(`Lslider${i}`);
        const Linput = document.getElementById(`Linput${i}`);
        const Rslider = document.getElementById(`Rslider${i}`);
        const Rinput = document.getElementById(`Rinput${i}`);
        const R1slider = document.getElementById(`R1slider${i}`);
        const R1input = document.getElementById(`R1input${i}`);
        const R2slider = document.getElementById(`R2slider${i}`);
        const R2input = document.getElementById(`R2input${i}`);
        const typeSelect = document.getElementById(`type${i}`);
        
        //initialize
        Lslider.value = g.L;
        Linput.value = g.L;
        Rslider.value = g.R;
        Rinput.value = g.R;
        R1slider.value = g.R1;
        R1input.value = g.R1;
        R2slider.value = g.R2;
        R2input.value = g.R2;
        typeSelect.value = g.type;

        updateGeometryControlState(g, i);

        typeSelect.addEventListener("change", () => {
            g.type = typeSelect.value;
            updateGeometry(g, i);
            updateGeometryControlState(g, i);
        });

        Lslider.addEventListener("input", () => {
            if (+Lslider.value<Math.abs(g.R2-g.R1) ) {
                g.L = Math.abs(g.R2-g.R1) + 1e-10;
                Lslider.value = Math.abs(g.R2-g.R1);
            }
            else {
                g.L = +Lslider.value;
            }
            Linput.value = +Lslider.value;
            updateGeometry(g, i);
        });
        Linput.addEventListener("change", () => {
            if (+Linput.value<Math.abs(g.R2-g.R1) ) {
                g.L = Math.abs(g.R2-g.R1) + 1e-10;
                Linput.value = Math.abs(g.R2-g.R1);
            }
            else {
                g.L = +Linput.value;
            }
            Lslider.value = +Linput.value;
            updateGeometry(g, i);
        })

        Rslider.addEventListener("input", () => {
            g.R = +Rslider.value;
            Rinput.value = +Rslider.value;
            updateGeometry(g, i);
        });
        Rinput.addEventListener("change", () => {
            g.R = +Rinput.value;
            Rslider.value = +Rinput.value;
            updateGeometry(g, i);
        })

        R1slider.addEventListener("input", () => {
            if (+R1slider.value<g.R2-g.L) {
                g.R1 = g.R2-g.L + 1e-10;
                R1slider.value = g.R2-g.L;
            }
            else if (+R1slider.value>g.R2+g.L) {
                g.R1 = g.R2+g.L - 1e-10;
                R1slider.value = g.R2+g.L;
            }
            else {
                g.R1 = +R1slider.value;
            }
            R1input.value = +R1slider.value;
            updateGeometry(g, i);
        });
        R1input.addEventListener("change", () => {
            if (+R1input.value<g.R2-g.L) {
                g.R1 = g.R2-g.L + 1e-10;
                R1input.value = g.R2-g.L;
            }
            else if (+R1input.value>g.R2+g.L) {
                g.R1 = g.R2+g.L - 1e-10;
                R1input.value = g.R2+g.L;
            }
            else {
                g.R1 = +R1input.value;
            }
            R1slider.value = +R1input.value;
            updateGeometry(g, i);
        })

        R2slider.addEventListener("input", () => {
            if (+R2slider.value<g.R1-g.L) {
                g.R2 = g.R1-g.L + 1e-10;
                R2slider.value = g.R1-g.L;
            }
            else if (+R2slider.value>g.R1+g.L) {
                g.R2 = g.R1+g.L - 1e-10;
                R2slider.value = g.R1+g.L;
            }
            else {
                g.R2 = +R2slider.value;
            }
            R2input.value = +R2slider.value;
            updateGeometry(g, i);
        });
        R2input.addEventListener("change", () => {
            if (+R2input.value<g.R1-g.L) {
                g.R2 = g.R1-g.L + 1e-10;
                R2input.value = g.R1-g.L;
            }
            else if (+R2input.value>g.R1+g.L) {
                g.R2 = g.R1+g.L - 1e-10;
                R2input.value = g.R1+g.L;
            }
            else {
                g.R2 = +R2input.value;
            }
            R2slider.value = +R2input.value;
            updateGeometry(g, i);
        })

    }

    geometries.forEach( (g, i) => {
        bindGeometryControls(g, i);
    })


    function updateVisualScale(g, i) {

        // update canvas geometry
        const dpr = window.devicePixelRatio;
        const canvas = animCanvases[i];
        const width = g.L * pxPerUnit + leftMargin + rightMargin;
        canvas.style.width = `${width}px`;
        canvas.style.height = `${2 * g.rMax * pxPerUnit}px`;
        canvas.width = canvas.clientWidth * dpr;
        canvas.height = canvas.clientHeight * dpr;

        const fieldPlots = document.getElementById(`field${i}`);
        fieldPlots.style.width = `${g.L * pxPerUnit + leftMargin + rightMargin}px`;

        Plotly.relayout(`field${i}`, {
            width: g.L * pxPerUnit + leftMargin + rightMargin,
            'xaxis.range': [0, g.L]
        });

        drawParticles(canvas, g);
    }

    function updateParticleDensity() {
        geometries.forEach( (g, i) => {
            initParticles(g);
            const canvas = animCanvases[i];
            drawParticles(canvas, g);
        })
    }

    function bindScaleControls() {
        const pxSlider = document.getElementById("pxPerUnitslider");
        const pxInput = document.getElementById("pxPerUnitinput");
        const timeSlider = document.getElementById("timeScaleslider");
        const timeInput = document.getElementById("timeScaleinput");
        const densitySlider = document.getElementById("particleDensityslider");
        const densityInput = document.getElementById("particleDensityinput");

        //initialize
        pxSlider.value = pxPerUnit;
        pxInput.value = pxPerUnit;
        timeSlider.value = timeScale;
        timeInput.value = timeScale;
        densitySlider.value = densityParticles;
        densityInput.value = densityParticles;

        pxSlider.addEventListener("input", () => {
            pxPerUnit = +pxSlider.value;
            pxInput.value = +pxSlider.value;
            geometries.forEach( (g, i) => {
                updateVisualScale(g, i);
            })  
        })
        pxInput.addEventListener("change", () => {
            pxPerUnit = +pxInput.value;
            pxSlider.value = +pxInput.value;
            geometries.forEach( (g, i) => {
                updateVisualScale(g, i);
            })  
        })

        timeSlider.addEventListener("input" , () => {
            timeScale = +timeSlider.value;
            timeInput.value = +timeSlider.value;
        })
        timeInput.addEventListener("change", () => {
            timeScale = +timeInput.value;
            timeSlider.value = +timeInput.value;
        })
        
        densitySlider.addEventListener("input" , () => {
            densityParticles = +densitySlider.value;
            densityInput.value = +densitySlider.value;
            updateParticleDensity();
        })
        densityInput.addEventListener("change", () => {
            densityParticles = +densityInput.value;
            densitySlider.value = +densityInput.value;
        })
    }
    bindScaleControls();


    function sliderToInput(t, min, max) {
        return Math.sign(t) * min*(max/min)**Math.abs(t);
    }
    function inputToSlider(input, min, max) {
        return Math.sign(input) * Math.log(Math.abs(input)/min)/Math.log(max/min);
    }


    function updateImpedance(g, i) {
        computeSpectrumAdaptive(g, fMin, fMax);
        g.computeField(f);
        globalpMax = 0;
        globaluMax = 0;

        geometries.forEach( (g, i) => {
            g.pMax = maxArray(g.pAbs);
            g.uMax = maxArray(g.uAbs);
            if (g.pMax > globalpMax) {
                globalpMax = g.pMax;
            }
            if (g.uMax > globaluMax) {
                globaluMax = g.uMax;
            }
        })
        

        const fieldPlots = document.getElementById(`field${i}`);
        fieldPlots.style.width = `${g.L * pxPerUnit + leftMargin + rightMargin}px`;

        // update plots
        Plotly.update(`field${i}`, {
            y: [g.p, g.u]
        });

        Plotly.update(`spectrum${i}`, {
            x: g.freqArray,
            y: g.ZinAbsArray
        })
    }

    function bindImpedanceControls() {
        geometries.forEach( (g, i) => {
            const ZrSlider = document.getElementById(`Zrslider${i}`);
            const ZrInput = document.getElementById(`Zrinput${i}`);
            const ZiSlider = document.getElementById(`Zislider${i}`);
            const ZiInput = document.getElementById(`Ziinput${i}`);

            Zrmin = 1e-5;
            Zrmax = 1e+8;
            ZrInput.min = Zrmin;
            ZrInput.max = Zrmax;
            Zimin = 1e-5;
            Zimax = 1e8;
            ZiInput.max = Zimax;
            ZiInput.min = -Zimax;
            

            //initialize
            ZrSlider.value = inputToSlider(g.Zr, Zrmin, Zrmax);
            ZrInput.value = g.Zr;
            ZiSlider.value = inputToSlider(g.Zi, Zimin, Zimax);
            ZiInput.value = g.Zi;


            ZrSlider.addEventListener("input", () => {
                g.Zr = sliderToInput(+ZrSlider.value, Zrmin, Zrmax);
                ZrInput.value = sliderToInput(+ZrSlider.value, Zrmin, Zrmax).toExponential(2);
                updateImpedance(g, i);
            })
            ZrInput.addEventListener("change", () => {
                g.Zr = +ZrInput.value;
                ZrSlider.value = inputToSlider(+ZrInput.value, Zrmin, Zrmax);
                updateImpedance(g, i);
            })
            
            ZiSlider.addEventListener("input", () => {
                g.Zi = sliderToInput(+ZiSlider.value, Zimin, Zimax);
                ZiInput.value = sliderToInput(+ZiSlider.value, Zimin, Zimax).toExponential(2);
                updateImpedance(g, i);
            })
            ZiInput.addEventListener("change", () => {
                g.Zi = +ZiInput.value;
                ZiSlider.value = inputToSlider(+ZiInput.value, Zimin, Zimax);
                updateImpedance(g, i);
            })
        })
        
    }
    bindImpedanceControls();

    function bindMovementScaleControls() {
        geometries.forEach( (g, i) => {
            const movementSlider = document.getElementById(`particleMovementSlider${i}`);
            const movementInput = document.getElementById(`particleMovementInput${i}`);

            //initialize
            g.scaleFactor = 1;
            movementMin = 1e-3;
            movementMax = 1e4;
            movementSlider.value = inputToSlider(g.scaleFactor, movementMin, movementMax);
            movementInput.value = g.scaleFactor;

            movementSlider.addEventListener("input", () => {
                g.scaleFactor = sliderToInput(+movementSlider.value, movementMin, movementMax);
                movementInput.value = sliderToInput(+movementSlider.value, movementMin, movementMax).toExponential(2);
            })
            movementInput.addEventListener("change", () => {
                g.scaleFactor = +movementInput.value;
                movementSlider.value = inputToSlider(+movementInput.value, movementMin, movementMax);
            })

        })
    }
    bindMovementScaleControls();




    const container = document.getElementById("plots");
    const fieldCharts = [];
    const animCharts = [];
    const spectrumDivs = [document.getElementById("spectrum0"), document.getElementById("spectrum1")];
    const fieldDivs = [document.getElementById("field0"), document.getElementById("field1")];
    const animCanvases = [document.getElementById("particles0"), document.getElementById("particles1")];

    const leftMargin = 75;
    const rightMargin = 75;
    
    geometries.forEach((g,i) => {

        // Field plot
        spectrumDivs[i].height = 200;
        const traceSpectrum = {
            x: Array.from(g.freqArray),
            y: Array.from(g.ZinAbsArray),
            mode: "lines",
            name: "Input Impedance",
        }
        const layoutSpectrum = {
            height: 200,
            margin: {
                t: 20,
                b: 20,
            },
            xaxis: {
                range: [0, fMax],
                title: {
                    text: 'f',
                },
                showline: true,
                mirror: true,
                automargin: true,
                ticks: "outside",
            },
            yaxis: {
                range: [-6, 12],
                type: "log",
                title: {
                    text: "Input Impedance",
                },
                showline: true,
                mirror: true,
            },
            shapes: [
                {
                    type: "line",
                    x0: f,
                    x1: f,
                    y0: 0,
                    y1: 1,
                    yref: "paper",
                    line: {
                        color: "red",
                        width: 2
                    }
                }
            ]
        }
        Plotly.newPlot(`spectrum${i}`, [traceSpectrum], layoutSpectrum);

        fieldDivs[i].height=250;
        const traceP = {
            x: Array.from(g.x),
            y: Array.from(g.p),
            mode: 'lines',
            name: 'Pressure'
        };

        const traceU = {
            x: Array.from(g.x),
            y: Array.from(g.u),
            mode: 'lines',
            name: 'Velocity',
            yaxis: 'y2'   // second axis
        };

        const layoutField = {
            width: g.L * pxPerUnit + leftMargin + rightMargin,
            height: 220,

            margin: {
                l: leftMargin,
                r: rightMargin,
                t: 20,
                b: 20,
            },

            legend: {
                yanchor: 'bottom',
                xanchor: 'center',
                orientation: "h",
                x: 0.5,
                y: 1.01,
            },

            xaxis: {
                range: [0, g.L],
                title: 'x',
                showline: true,
                mirror: true,
                automargin: true,
                ticks: "outside",
            },

            yaxis: {
                showline: false,
                title: {
                    text: 'Pressure [Pa]',
                },
                automargin: false,
                ticks: "outside",
                range: [synchronize ? -globalpMax : -g.pMax, synchronize ? globalpMax : g.pMax],
            },

            yaxis2: {
                showline: true,
                title: {
                    text: 'Volume Velocity [m<sup>3</sup>/s]'
                },
                overlaying: 'y',
                side: 'right',
                automargin: false,
                ticks: "outside",
                range: [synchronize ? -globaluMax : -g.uMax, synchronize ? globaluMax : g.uMax]
            },
        };

        Plotly.newPlot(`field${i}`, [traceP, traceU], layoutField);
        document.getElementById(`spectrum${i}`).on('plotly_click', (data) => {
            f = data.points[0].x;
            document.getElementById("freqSlider").value = f;
            updateFrequency(f);
        });

        const rMax = Math.max(g.R ?? 0, g.R1 ?? 0, g.R2 ?? 0);
        const width = g.L * pxPerUnit + leftMargin + rightMargin;
        
        const dpr = window.devicePixelRatio;
        animCanvases[i].style.width = `${width}px`;
        animCanvases[i].style.height = `${2*g.rMax * pxPerUnit}px`;
        animCanvases[i].width = animCanvases[i].clientWidth * dpr;
        animCanvases[i].height = animCanvases[i].clientHeight * dpr;
        drawParticles(animCanvases[i], g);
    })



   

   


    let start;
    geometries.forEach((g,i) => {
        g.pDynamic = new Float64Array(g.N);
        g.uDynamic = new Float64Array(g.N);
    })

    function step(timestamp) {

        if (start === undefined) {
            start = timestamp;
        }
        const t = (timestamp - start)/1000*timeScale; //elapsed time since start
        const synchronize = synchronizeToggle.checked;

        geometries.forEach((g, i) => {  
            const particles = g.particles;
            if (g.type=="cylinder") {
                particles.forEach( (p, j) => {
                    const dx = displacement(g, p.x0, t);
                    p.x = p.x0 + dx * g.scaleFactor; 
                })

            }

            else if (g.type=="cone") {
                if (Math.abs(g.R1-g.R2) < 0.00001) {
                    particles.forEach( (p, j) => {
                        const dx = displacement(g, p.x0, t);
                        p.x = p.x0 + dx * g.scaleFactor; 
                    })
                }
                else {
                    particles.forEach( (p, j) => {
                        const dr = displacement(g, p.r0, t);
                        p.r =  g.dl>0 ? Math.max(-g.dl, p.r0 + dr*g.scaleFactor) : Math.min(-g.dl, p.r0 + dr*g.scaleFactor);
                        p.x =  (p.r+g.dl)*Math.cos(p.phi0)- g.dl;
                        p.y =  (p.r+g.dl)*Math.sin(p.phi0);
                    })
                }
                
            }
            drawParticles(animCanvases[i], g);

        for (let j = 0; j < g.N; j++) {
            g.pDynamic[j] = Math.cos(g.pArg[j] + 2*Math.PI*f*t) * g.pAbs[j];
            g.uDynamic[j] = Math.cos(g.uArg[j] + 2*Math.PI*f*t) * g.uAbs[j];
        }

        Plotly.update(`field${i}`,
            {   
                x: [g.x, g.x],
                y: [g.pDynamic, g.uDynamic],                
            }, 
            {
                'yaxis.range': [synchronize ? -globalpMax : -g.pMax,
                                synchronize ?  globalpMax :  g.pMax],
                'yaxis2.range': [synchronize ? -globaluMax : -g.uMax,
                                synchronize ?  globaluMax :  g.uMax]
            }
            );

        })
       
    
    requestAnimationFrame(step);


    }
    requestAnimationFrame(step);




//     Module._free(xPtr);
//     Module._free(p1AbsPtr);
//     Module._free(p1ArgPtr);
//     Module._free(u1AbsPtr);
//     Module._free(u1ArgPtr);
//     Module._free(p1Ptr);
//     Module._free(u1Ptr);  
}
)

