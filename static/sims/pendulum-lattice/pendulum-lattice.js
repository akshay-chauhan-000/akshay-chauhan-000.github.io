/* ==========================================================================
   Double-Pendulum Lattice
   --------------------------------------------------------------------------
   An n x n lattice of identical double pendulums (m1 = m2 = l1 = l2 = 1),
   each coupled to its four nearest neighbours through torsional springs on
   the corresponding joint angles, with periodic boundaries.

   Per site, with D = th1 - th2, the equations of motion come from the
   standard double-pendulum Lagrangian written as M(q) q'' = F(q, q'):

       M = | 2      cosD |      F1 = -w2^2 sinD - 2 g sin(th1) + tau1
           | cosD   1    |      F2 =  w1^2 sinD -   g sin(th2) + tau2

   The coupling and damping enter as generalized torques,

       tau1 = K1 * sum_j sin(th1_j - th1_i) - gamma w1
       tau2 = K2 * sum_j sin(th2_j - th2_i) - gamma w2

   so that for gamma = 0 the full lattice is Hamiltonian: site energies plus
   the bond potential sum K (1 - cos(dth)) are conserved. Solving M q'' = F:

       det  = 2 - cos^2 D
       w1'  = (F1 - cosD * F2) / det
       w2'  = (2 F2 - cosD * F1) / det

   Integration is classical RK4 over the entire 4 n^2 - dimensional state,
   with the four trig values per site cached once per stage; every coupling
   term is then reconstructed from the cache via angle-difference identities,
   so the whole derivative costs exactly 4 trig calls per site.

   No dependencies, no build step. The engine (Lattice) is DOM-free and
   exported on window.PendulumLattice for headless use; the UI mounts into
   [data-sim="pendulum-lattice"].
   ========================================================================== */

(function () {
  'use strict';

  var TWO_PI = 2 * Math.PI;
  var INK = '#0b0d12';

  /* Energy ramp shared with the other sims: ink -> indigo -> copper -> sand */
  var RAMP = [
    [0.00, 11, 13, 18],
    [0.25, 36, 48, 94],
    [0.50, 86, 98, 142],
    [0.75, 189, 141, 87],
    [1.00, 246, 231, 178]
  ];

  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  function mulberry32(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function rampColor(e, out) {
    e = clamp(e, 0, 1);
    for (var s = 1; s < RAMP.length; s++) {
      if (e <= RAMP[s][0]) {
        var a = RAMP[s - 1], b = RAMP[s];
        var t = (e - a[0]) / (b[0] - a[0]);
        out[0] = (a[1] + (b[1] - a[1]) * t) | 0;
        out[1] = (a[2] + (b[2] - a[2]) * t) | 0;
        out[2] = (a[3] + (b[3] - a[3]) * t) | 0;
        return out;
      }
    }
    out[0] = 246; out[1] = 231; out[2] = 178;
    return out;
  }

  /* Cyclic phase colour: smooth hue wheel from three shifted cosines. */
  function phaseColor(theta, out) {
    out[0] = (127.5 + 110 * Math.cos(theta)) | 0;
    out[1] = (127.5 + 110 * Math.cos(theta - 2.0944)) | 0;
    out[2] = (127.5 + 110 * Math.cos(theta + 2.0944)) | 0;
    return out;
  }

  // -------------------------------------------------------------------------
  // Lattice engine
  // -------------------------------------------------------------------------

  function Lattice(n, seed) {
    this.g = 9.8;
    this.K1 = 1.2;
    this.K2 = 0.6;
    this.gamma = 0.04;
    this.t = 0;
    this.alloc(n);
    this.reset(seed || 1, 0.1);
  }

  Lattice.prototype.alloc = function (n) {
    var N = n * n;
    this.n = n;
    this.N = N;
    this.state = [];                        // [th1, w1, th2, w2]
    this.tmp = [];
    this.k = [[], [], [], []];              // RK4 stage derivatives
    for (var c = 0; c < 4; c++) {
      this.state.push(new Float64Array(N));
      this.tmp.push(new Float64Array(N));
      for (var s = 0; s < 4; s++) this.k[s].push(new Float64Array(N));
    }
    this.s1 = new Float64Array(N);
    this.c1 = new Float64Array(N);
    this.s2 = new Float64Array(N);
    this.c2 = new Float64Array(N);
    // Periodic von Neumann neighbourhood.
    this.nb = new Int32Array(4 * N);
    for (var y = 0; y < n; y++) {
      for (var x = 0; x < n; x++) {
        var i = y * n + x;
        this.nb[4 * i] = y * n + ((x + 1) % n);
        this.nb[4 * i + 1] = y * n + ((x + n - 1) % n);
        this.nb[4 * i + 2] = ((y + 1) % n) * n + x;
        this.nb[4 * i + 3] = ((y + n - 1) % n) * n + x;
      }
    }
  };

  Lattice.prototype.reset = function (seed, amp) {
    var rnd = mulberry32(seed);
    for (var i = 0; i < this.N; i++) {
      this.state[0][i] = amp * (2 * rnd() - 1);
      this.state[1][i] = 0;
      this.state[2][i] = amp * (2 * rnd() - 1);
      this.state[3][i] = 0;
    }
    this.t = 0;
  };

  Lattice.prototype.randomize = function (seed) {
    var rnd = mulberry32(seed);
    for (var i = 0; i < this.N; i++) {
      this.state[0][i] = Math.PI * (2 * rnd() - 1);
      this.state[1][i] = 0;
      this.state[2][i] = Math.PI * (2 * rnd() - 1);
      this.state[3][i] = 0;
    }
    this.t = 0;
  };

  Lattice.prototype.kick = function (x, y, dw) {
    var i = ((y % this.n) + this.n) % this.n * this.n + ((x % this.n) + this.n) % this.n;
    this.state[1][i] += dw;
    this.state[3][i] += 0.5 * dw;
  };

  /* Nearest-neighbour resampling: the pattern survives a live change of n. */
  Lattice.prototype.resample = function (newN) {
    var old = this.state, oldN = this.n;
    var keep = [old[0], old[1], old[2], old[3]];
    this.alloc(newN);
    for (var y = 0; y < newN; y++) {
      var oy = Math.min(oldN - 1, Math.floor(y * oldN / newN));
      for (var x = 0; x < newN; x++) {
        var ox = Math.min(oldN - 1, Math.floor(x * oldN / newN));
        var i = y * newN + x, o = oy * oldN + ox;
        for (var c = 0; c < 4; c++) this.state[c][i] = keep[c][o];
      }
    }
  };

  /* d/dt of (th1, w1, th2, w2) for the whole lattice. */
  Lattice.prototype.deriv = function (y, out) {
    var N = this.N, g = this.g, K1 = this.K1, K2 = this.K2, gam = this.gamma;
    var th1 = y[0], w1 = y[1], th2 = y[2], w2 = y[3];
    var s1 = this.s1, c1 = this.c1, s2 = this.s2, c2 = this.c2, nb = this.nb;
    var i, j;

    for (i = 0; i < N; i++) {
      s1[i] = Math.sin(th1[i]);
      c1[i] = Math.cos(th1[i]);
      s2[i] = Math.sin(th2[i]);
      c2[i] = Math.cos(th2[i]);
    }

    for (i = 0; i < N; i++) {
      // sin/cos of (th1 - th2) from the per-angle cache.
      var sinD = s1[i] * c2[i] - c1[i] * s2[i];
      var cosD = c1[i] * c2[i] + s1[i] * s2[i];

      // Coupling: sum over neighbours of sin(th_j - th_i), both joints.
      var C1 = 0, C2 = 0, b = 4 * i;
      for (var q = 0; q < 4; q++) {
        j = nb[b + q];
        C1 += s1[j] * c1[i] - c1[j] * s1[i];
        C2 += s2[j] * c2[i] - c2[j] * s2[i];
      }

      var F1 = -w2[i] * w2[i] * sinD - 2 * g * s1[i] + K1 * C1 - gam * w1[i];
      var F2 = w1[i] * w1[i] * sinD - g * s2[i] + K2 * C2 - gam * w2[i];
      var det = 2 - cosD * cosD;

      out[0][i] = w1[i];
      out[1][i] = (F1 - cosD * F2) / det;
      out[2][i] = w2[i];
      out[3][i] = (2 * F2 - cosD * F1) / det;
    }
  };

  Lattice.prototype.step = function (dt) {
    var y = this.state, tmp = this.tmp, k = this.k;
    var N = this.N, c, i;

    this.deriv(y, k[0]);
    for (c = 0; c < 4; c++) {
      for (i = 0; i < N; i++) tmp[c][i] = y[c][i] + 0.5 * dt * k[0][c][i];
    }
    this.deriv(tmp, k[1]);
    for (c = 0; c < 4; c++) {
      for (i = 0; i < N; i++) tmp[c][i] = y[c][i] + 0.5 * dt * k[1][c][i];
    }
    this.deriv(tmp, k[2]);
    for (c = 0; c < 4; c++) {
      for (i = 0; i < N; i++) tmp[c][i] = y[c][i] + dt * k[2][c][i];
    }
    this.deriv(tmp, k[3]);
    for (c = 0; c < 4; c++) {
      var yc = y[c], k1 = k[0][c], k2 = k[1][c], k3 = k[2][c], k4 = k[3][c];
      for (i = 0; i < N; i++) {
        yc[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      }
    }
    this.t += dt;
    return isFinite(y[1][0]) && isFinite(y[0][(N - 1) | 0]);
  };

  /* Site energy: kinetic + gravitational, offset so the hanging state is 0. */
  Lattice.prototype.energy = function (out) {
    var y = this.state, g = this.g;
    for (var i = 0; i < this.N; i++) {
      var cosD = Math.cos(y[0][i] - y[2][i]);
      out[i] = y[1][i] * y[1][i] + 0.5 * y[3][i] * y[3][i]
        + y[1][i] * y[3][i] * cosD
        + 2 * g * (1 - Math.cos(y[0][i]))
        + g * (1 - Math.cos(y[2][i]));
    }
  };

  /* Kuramoto-style order parameters for both joints. */
  Lattice.prototype.order = function () {
    var y = this.state, a1 = 0, b1 = 0, a2 = 0, b2 = 0;
    for (var i = 0; i < this.N; i++) {
      a1 += Math.cos(y[0][i]); b1 += Math.sin(y[0][i]);
      a2 += Math.cos(y[2][i]); b2 += Math.sin(y[2][i]);
    }
    return [
      Math.sqrt(a1 * a1 + b1 * b1) / this.N,
      Math.sqrt(a2 * a2 + b2 * b2) / this.N
    ];
  };

  // -------------------------------------------------------------------------
  // UI
  // -------------------------------------------------------------------------

  var FIELDS = [
    ['th1', 'theta-1 phase'],
    ['th2', 'theta-2 phase'],
    ['energy', 'site energy'],
    ['sync', 'neighbour sync']
  ];

  var SLIDERS = [
    { key: 'n', label: 'Lattice size n', min: 4, max: 32, step: 1, fmt: function (v) { return v + ' x ' + v; } },
    { key: 'K1', label: 'Coupling K1 (upper)', min: 0, max: 5, step: 0.05, fmt: function (v) { return v.toFixed(2); } },
    { key: 'K2', label: 'Coupling K2 (lower)', min: 0, max: 5, step: 0.05, fmt: function (v) { return v.toFixed(2); } },
    { key: 'g', label: 'Gravity g', min: 0, max: 20, step: 0.1, fmt: function (v) { return v.toFixed(1); } },
    { key: 'gamma', label: 'Damping', min: 0, max: 1, step: 0.01, fmt: function (v) { return v.toFixed(2); } },
    { key: 'dt', label: 'Time step', min: 0.002, max: 0.02, step: 0.001, fmt: function (v) { return v.toFixed(3); } },
    { key: 'spf', label: 'Steps / frame', min: 1, max: 12, step: 1, fmt: function (v) { return v.toFixed(0); } }
  ];

  var TEMPLATE = [
    '<div class="dp-stage">',
    '  <div class="dp-hud">',
    '    <span>t <b data-hud="t">0.0</b></span>',
    '    <span>&lt;E&gt; <b data-hud="e">0.00</b></span>',
    '    <span>R1 <b data-hud="r1">0.00</b></span>',
    '    <span>R2 <b data-hud="r2">0.00</b></span>',
    '    <span><b data-hud="fps">0</b> fps</span>',
    '  </div>',
    '  <div class="dp-panel">',
    '    <canvas class="dp-canvas dp-canvas--lattice" role="img" aria-label="Lattice of double pendulums"></canvas>',
    '    <span class="dp-caption">pendulums &middot; click to kick</span>',
    '  </div>',
    '  <div class="dp-panel">',
    '    <canvas class="dp-canvas dp-canvas--map" role="img" aria-label="Field heatmap"></canvas>',
    '    <span class="dp-caption" data-cap="map">field</span>',
    '  </div>',
    '</div>',
    '<div class="dp-bar">',
    '  <button type="button" data-act="run">Pause</button>',
    '  <button type="button" data-act="reset">Reset</button>',
    '  <button type="button" data-act="random">Randomize</button>',
    '  <button type="button" data-act="kick">Kick centre</button>',
    '  <span class="dp-spacer"></span>',
    '  <button type="button" data-act="link">Copy link</button>',
    '  <span class="dp-status" role="status"></span>',
    '</div>',
    '<div class="dp-body">',
    '  <div class="dp-controls">',
    '    <div class="dp-field"><label><span>Heatmap field</span></label><select data-sel="field"></select></div>',
    '  </div>',
    '</div>'
  ].join('\n');

  function mount(root) {
    if (root.dataset.dpMounted) return;
    root.dataset.dpMounted = '1';
    root.classList.add('dp');
    root.innerHTML = TEMPLATE;

    var latCanvas = root.querySelector('.dp-canvas--lattice');
    var mapCanvas = root.querySelector('.dp-canvas--map');
    var latCtx = latCanvas.getContext('2d', { alpha: false });
    var mapCtx = mapCanvas.getContext('2d', { alpha: false });
    if (!latCtx || !mapCtx) return;

    var conf = {
      n: 16, K1: 1.2, K2: 0.6, g: 9.8, gamma: 0.04,
      dt: 0.005, spf: 4, field: 'th1', seed: 20260809
    };

    var statusEl = root.querySelector('.dp-status');
    var runBtn = root.querySelector('[data-act="run"]');
    var mapCap = root.querySelector('[data-cap="map"]');
    var hud = {
      t: root.querySelector('[data-hud="t"]'),
      e: root.querySelector('[data-hud="e"]'),
      r1: root.querySelector('[data-hud="r1"]'),
      r2: root.querySelector('[data-hud="r2"]'),
      fps: root.querySelector('[data-hud="fps"]')
    };

    // ---- controls ---------------------------------------------------------

    var selField = root.querySelector('[data-sel="field"]');
    FIELDS.forEach(function (p) {
      var o = document.createElement('option');
      o.value = p[0];
      o.textContent = p[1];
      selField.appendChild(o);
    });
    selField.addEventListener('change', function () {
      conf.field = selField.value;
      syncCaption();
    });

    var controlsHost = root.querySelector('.dp-controls');
    var readouts = {};
    SLIDERS.forEach(function (ctl) {
      var wrap = document.createElement('div');
      wrap.className = 'dp-field';
      var label = document.createElement('label');
      var name = document.createElement('span');
      name.textContent = ctl.label;
      var value = document.createElement('span');
      label.appendChild(name);
      label.appendChild(value);
      var input = document.createElement('input');
      input.type = 'range';
      input.min = String(ctl.min);
      input.max = String(ctl.max);
      input.step = String(ctl.step);
      wrap.appendChild(label);
      wrap.appendChild(input);
      controlsHost.appendChild(wrap);
      readouts[ctl.key] = { input: input, value: value, ctl: ctl };

      input.addEventListener('input', function () {
        var v = parseFloat(input.value);
        if (!isFinite(v)) return;
        conf[ctl.key] = v;
        value.textContent = ctl.fmt(v);
        if (ctl.key === 'n') {
          lattice.resample(Math.round(v));   // live: pattern is resampled
          rebuildBuffers();
        } else if (ctl.key !== 'dt' && ctl.key !== 'spf') {
          applyParams();
        }
      });
    });

    function syncControls() {
      SLIDERS.forEach(function (ctl) {
        var r = readouts[ctl.key];
        r.input.value = String(conf[ctl.key]);
        r.value.textContent = ctl.fmt(conf[ctl.key]);
      });
      selField.value = conf.field;
      syncCaption();
    }

    function syncCaption() {
      var lbl = conf.field;
      for (var i = 0; i < FIELDS.length; i++) {
        if (FIELDS[i][0] === conf.field) lbl = FIELDS[i][1];
      }
      mapCap.textContent = 'field: ' + lbl;
    }

    // ---- engine + buffers -------------------------------------------------

    var lattice = new Lattice(conf.n, conf.seed);
    var running = true;
    var eField = null, sField = null;
    var mapImage = null;
    var mapBack = document.createElement('canvas');
    var eScale = 1;                          // smoothed energy normalisation
    var frames = 0, fpsAccum = 0, lastStamp = 0;

    function applyParams() {
      lattice.g = conf.g;
      lattice.K1 = conf.K1;
      lattice.K2 = conf.K2;
      lattice.gamma = conf.gamma;
    }

    function rebuildBuffers() {
      var n = lattice.n;
      eField = new Float64Array(n * n);
      sField = new Float64Array(n * n);
      mapBack.width = n;
      mapBack.height = n;
      mapImage = mapBack.getContext('2d').createImageData(n, n);
    }

    // ---- rendering --------------------------------------------------------

    function fitCanvas(canvas) {
      var dpr = Math.min(2, window.devicePixelRatio || 1);
      var rect = canvas.getBoundingClientRect();
      var w = Math.max(1, Math.round(rect.width * dpr));
      var h = Math.max(1, Math.round(rect.height * dpr));
      if (w !== canvas.width || h !== canvas.height) {
        canvas.width = w;
        canvas.height = h;
      }
    }

    var HUE_BINS = 24;
    var binPaths = [];                       // reused arrays of segment coords
    for (var hb = 0; hb < HUE_BINS; hb++) binPaths.push([]);
    var rgbTmp = [0, 0, 0];
    var binColor = [];
    for (var hb2 = 0; hb2 < HUE_BINS; hb2++) {
      phaseColor((hb2 + 0.5) / HUE_BINS * TWO_PI - Math.PI, rgbTmp);
      binColor.push('rgb(' + rgbTmp[0] + ',' + rgbTmp[1] + ',' + rgbTmp[2] + ')');
    }

    function drawLattice() {
      fitCanvas(latCanvas);
      var W = latCanvas.width, H = latCanvas.height;
      var n = lattice.n, cell = W / n;
      var y = lattice.state;
      latCtx.fillStyle = INK;
      latCtx.fillRect(0, 0, W, H);

      var arm = cell * 0.23;
      var lw = Math.max(1, cell * 0.07);
      var b;
      for (b = 0; b < HUE_BINS; b++) binPaths[b].length = 0;

      for (var gy = 0; gy < n; gy++) {
        for (var gx = 0; gx < n; gx++) {
          var i = gy * n + gx;
          var px = (gx + 0.5) * cell;
          var py = (gy + 0.5) * H / n;
          var x1 = px + arm * Math.sin(y[0][i]);
          var y1 = py + arm * Math.cos(y[0][i]);
          var x2 = x1 + arm * Math.sin(y[2][i]);
          var y2 = y1 + arm * Math.cos(y[2][i]);
          // Bin by th1 phase so strokes batch into HUE_BINS paths.
          var ang = y[0][i] % TWO_PI;
          if (ang < -Math.PI) ang += TWO_PI;
          else if (ang > Math.PI) ang -= TWO_PI;
          b = clamp(((ang + Math.PI) / TWO_PI * HUE_BINS) | 0, 0, HUE_BINS - 1);
          binPaths[b].push(px, py, x1, y1, x2, y2);
        }
      }

      latCtx.lineWidth = lw;
      latCtx.lineCap = 'round';
      for (b = 0; b < HUE_BINS; b++) {
        var seg = binPaths[b];
        if (!seg.length) continue;
        latCtx.strokeStyle = binColor[b];
        latCtx.beginPath();
        for (var s = 0; s < seg.length; s += 6) {
          latCtx.moveTo(seg[s], seg[s + 1]);
          latCtx.lineTo(seg[s + 2], seg[s + 3]);
          latCtx.lineTo(seg[s + 4], seg[s + 5]);
        }
        latCtx.stroke();
      }
    }

    function drawMap() {
      fitCanvas(mapCanvas);
      var n = lattice.n, N = n * n;
      var y = lattice.state, nb = lattice.nb;
      var px = mapImage.data;
      var rgb = rgbTmp;
      var i, v;

      if (conf.field === 'energy') {
        lattice.energy(eField);
        var mx = 1e-6;
        for (i = 0; i < N; i++) if (eField[i] > mx) mx = eField[i];
        eScale = Math.max(0.92 * eScale, mx);     // smooth the normalisation
        for (i = 0; i < N; i++) {
          rampColor(eField[i] / eScale, rgb);
          px[4 * i] = rgb[0]; px[4 * i + 1] = rgb[1]; px[4 * i + 2] = rgb[2]; px[4 * i + 3] = 255;
        }
      } else if (conf.field === 'sync') {
        for (i = 0; i < N; i++) {
          var acc = 0;
          for (var q = 0; q < 4; q++) {
            var j = nb[4 * i + q];
            acc += Math.cos(y[0][i] - y[0][j]);
          }
          sField[i] = acc * 0.25;                  // [-1, 1]
        }
        for (i = 0; i < N; i++) {
          rampColor(0.5 * (sField[i] + 1), rgb);
          px[4 * i] = rgb[0]; px[4 * i + 1] = rgb[1]; px[4 * i + 2] = rgb[2]; px[4 * i + 3] = 255;
        }
      } else {
        var src = conf.field === 'th2' ? y[2] : y[0];
        for (i = 0; i < N; i++) {
          v = src[i] % TWO_PI;
          phaseColor(v, rgb);
          px[4 * i] = rgb[0]; px[4 * i + 1] = rgb[1]; px[4 * i + 2] = rgb[2]; px[4 * i + 3] = 255;
        }
      }

      mapBack.getContext('2d').putImageData(mapImage, 0, 0);
      mapCtx.imageSmoothingEnabled = false;
      mapCtx.fillStyle = INK;
      mapCtx.fillRect(0, 0, mapCanvas.width, mapCanvas.height);
      mapCtx.drawImage(mapBack, 0, 0, mapCanvas.width, mapCanvas.height);
    }

    function updateHUD(stamp) {
      var ord = lattice.order();
      hud.t.textContent = lattice.t.toFixed(1);
      lattice.energy(eField);
      var mean = 0;
      for (var i = 0; i < lattice.N; i++) mean += eField[i];
      hud.e.textContent = (mean / lattice.N).toFixed(2);
      hud.r1.textContent = ord[0].toFixed(2);
      hud.r2.textContent = ord[1].toFixed(2);
      if (stamp) {
        fpsAccum++;
        if (stamp - lastStamp > 500) {
          hud.fps.textContent = String(Math.round(fpsAccum * 1000 / (stamp - lastStamp)));
          fpsAccum = 0;
          lastStamp = stamp;
        }
      }
    }

    // ---- main loop --------------------------------------------------------

    var statusTimer = 0;
    function say(msg) {
      statusEl.textContent = msg;
      clearTimeout(statusTimer);
      statusTimer = setTimeout(function () { statusEl.textContent = ''; }, 2600);
    }

    function frame(stamp) {
      if (running) {
        applyParams();
        var ok = true;
        for (var s = 0; s < conf.spf && ok; s++) ok = lattice.step(conf.dt);
        if (!ok) {
          lattice.reset(conf.seed, 0.1);
          say('numerical blow-up - state reset (lower dt or gravity)');
        }
      }
      drawLattice();
      drawMap();
      updateHUD(stamp);
      requestAnimationFrame(frame);
    }

    // ---- interaction ------------------------------------------------------

    latCanvas.addEventListener('pointerdown', function (ev) {
      ev.preventDefault();
      var rect = latCanvas.getBoundingClientRect();
      var gx = clamp(Math.floor((ev.clientX - rect.left) / rect.width * lattice.n), 0, lattice.n - 1);
      var gy = clamp(Math.floor((ev.clientY - rect.top) / rect.height * lattice.n), 0, lattice.n - 1);
      lattice.kick(gx, gy, ev.shiftKey ? -8 : 8);
    });

    root.addEventListener('click', function (ev) {
      var btn = ev.target.closest ? ev.target.closest('[data-act]') : null;
      if (!btn) return;
      var act = btn.dataset.act;
      if (act === 'run') {
        running = !running;
        runBtn.textContent = running ? 'Pause' : 'Run';
      } else if (act === 'reset') {
        conf.seed = (conf.seed + 1) >>> 0;
        lattice.reset(conf.seed, 0.1);
      } else if (act === 'random') {
        conf.seed = (Math.random() * 0xffffffff) >>> 0;
        lattice.randomize(conf.seed);
      } else if (act === 'kick') {
        lattice.kick(lattice.n >> 1, lattice.n >> 1, 8);
      } else if (act === 'link') {
        var url = buildURL();
        try { history.replaceState(null, '', url); } catch (e) { /* file:// */ }
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(url).then(
            function () { say('Link copied'); },
            function () { say('Link is in the address bar'); }
          );
        } else {
          say('Link is in the address bar');
        }
      }
    });

    // ---- state <-> URL ----------------------------------------------------

    function readHash() {
      var raw = (location.hash || '').replace(/^#/, '');
      if (raw.indexOf('dp=1') === -1) return;
      var p;
      try { p = new URLSearchParams(raw); } catch (e) { return; }
      var num = function (key, target, lo, hi) {
        if (!p.has(key)) return;
        var v = parseFloat(p.get(key));
        if (isFinite(v)) conf[target] = clamp(v, lo, hi);
      };
      num('n', 'n', 4, 32);
      num('k1', 'K1', 0, 5);
      num('k2', 'K2', 0, 5);
      num('g', 'g', 0, 20);
      num('gm', 'gamma', 0, 1);
      num('dt', 'dt', 0.002, 0.02);
      num('sf', 'spf', 1, 12);
      if (p.has('f')) {
        for (var i = 0; i < FIELDS.length; i++) {
          if (FIELDS[i][0] === p.get('f')) conf.field = p.get('f');
        }
      }
      if (p.has('s')) {
        var seed = parseInt(p.get('s'), 10);
        if (isFinite(seed)) conf.seed = seed >>> 0;
      }
      conf.n = Math.round(conf.n);
    }

    function buildURL() {
      var p = new URLSearchParams();
      p.set('dp', '1');
      p.set('n', String(lattice.n));
      p.set('k1', conf.K1.toFixed(2));
      p.set('k2', conf.K2.toFixed(2));
      p.set('g', conf.g.toFixed(1));
      p.set('gm', conf.gamma.toFixed(2));
      p.set('dt', conf.dt.toFixed(3));
      p.set('sf', String(Math.round(conf.spf)));
      p.set('f', conf.field);
      p.set('s', String(conf.seed));
      return location.origin + location.pathname + location.search + '#' + p.toString();
    }

    // ---- boot -------------------------------------------------------------

    readHash();
    lattice = new Lattice(conf.n, conf.seed);
    applyParams();
    rebuildBuffers();
    syncControls();

    if (typeof ResizeObserver !== 'undefined') {
      new ResizeObserver(function () { /* next frame redraws at new size */ }).observe(latCanvas);
    }

    requestAnimationFrame(frame);

    root._pendulumLattice = { getLattice: function () { return lattice; } };
  }

  // -------------------------------------------------------------------------

  if (typeof window !== 'undefined') {
    window.PendulumLattice = { Lattice: Lattice, mulberry32: mulberry32 };
  }
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { Lattice: Lattice, mulberry32: mulberry32 };
  }
  if (typeof document !== 'undefined') {
    var boot = function () {
      var nodes = document.querySelectorAll('[data-sim="pendulum-lattice"]');
      for (var i = 0; i < nodes.length; i++) mount(nodes[i]);
    };
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', boot);
    } else {
      boot();
    }
  }
})();
