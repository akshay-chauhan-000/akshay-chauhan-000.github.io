/* ==========================================================================
   Particle Life
   --------------------------------------------------------------------------
   N point particles of K species on a periodic (toroidal) domain. Every pair
   within a cutoff r_max interacts through a piecewise-linear radial force

       f(r; a) =  r/beta - 1                              ,  0 <= r < beta
                  a * (1 - |2r - 1 - beta| / (1 - beta))  ,  beta <= r < 1
                  0                                       ,  r >= 1

   with r = |r_ij| / r_max. The first branch is a species-independent core
   repulsion that keeps particles from collapsing; the second is a tent
   function scaled by the entry A[s_i][s_j] of the interaction matrix.

   A is *not* required to be symmetric, so the force law violates Newton's
   third law: species i can chase species j while j flees i. That single
   relaxation is what produces the self-propelled, non-equilibrium structures
   the model is known for. Dynamics are overdamped (velocity-Verlet without
   the acceleration half-step is pointless here), integrated as

       v <- v * 2^(-dt/tau) + F dt ,   x <- x + v dt   (mod L)

   where tau is the velocity half-life. Neighbour search is a uniform hash
   grid with cell size >= r_max, giving O(N) per step at fixed density.

   No dependencies, no build step. The engine (ParticleLife) is DOM-free and
   exported for headless use; the UI below it mounts into [data-sim].
   ========================================================================== */

(function () {
  'use strict';

  var TAU = Math.PI * 2;
  var INK = '#0b0d12';

  var PALETTE = [
    '#ff5d73', '#ffa64d', '#ffe066', '#5ddb8b',
    '#3fc8d8', '#6f9dff', '#b57bff', '#ff7ad1'
  ];

  var DEFAULTS = {
    seed: 20260807,
    species: 6,
    count: 1600,
    rMax: 0.075,      // cutoff, in units of the domain height
    beta: 0.30,       // core radius as a fraction of rMax
    forceScale: 10,   // overall force gain
    halfLife: 0.035,  // seconds for velocity to decay by half
    symmetry: 0,      // 0 = fully asymmetric matrix, 1 = symmetric
    dt: 1 / 120,
    maxSpeed: 4       // domain heights per second; blow-up guard only
  };

  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  /* Deterministic PRNG (mulberry32) — a seed reproduces a run exactly. */
  function mulberry32(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // -------------------------------------------------------------------------
  // Engine
  // -------------------------------------------------------------------------

  function ParticleLife(options) {
    this.config = {};
    for (var key in DEFAULTS) this.config[key] = DEFAULTS[key];
    if (options) for (var k in options) this.config[k] = options[k];

    this.w = 1;   // domain width  (tracks the canvas aspect ratio)
    this.h = 1;   // domain height (fixed: the unit of length)
    this.matrix = null;
    this.pointer = { active: false, x: 0, y: 0, radius: 0.14, strength: 3 };
    this._grid = { gx: 0, gy: 0, order: new Int32Array(0) };
    this.reset();
  }

  Object.defineProperty(ParticleLife.prototype, 'n', {
    get: function () { return this.config.count; }
  });
  Object.defineProperty(ParticleLife.prototype, 'k', {
    get: function () { return this.config.species; }
  });

  /* Draw A from the seed. Off-diagonal pairs are blended toward symmetry by
     config.symmetry, so the user can dial in how far from Newton's third law
     the system sits. */
  ParticleLife.prototype.randomizeMatrix = function (seed) {
    var k = this.k;
    var rng = mulberry32(seed === undefined ? (this.config.seed ^ 0x9e3779b9) : seed);
    var sym = clamp(this.config.symmetry, 0, 1);
    var m = new Float32Array(k * k);

    for (var i = 0; i < k; i++) {
      for (var j = 0; j <= i; j++) {
        var a = rng() * 2 - 1;
        var b = rng() * 2 - 1;
        m[i * k + j] = a;
        m[j * k + i] = (i === j) ? a : (a * sym + b * (1 - sym));
      }
    }
    this.matrix = m;
    return m;
  };

  /* Re-seed positions. The matrix is preserved unless its size no longer
     matches the species count. */
  ParticleLife.prototype.reset = function (seed) {
    var cfg = this.config;
    if (seed !== undefined) cfg.seed = seed >>> 0;

    var n = cfg.count;
    var k = cfg.species;
    var rng = mulberry32(cfg.seed);

    this.px = new Float32Array(n);
    this.py = new Float32Array(n);
    this.vx = new Float32Array(n);
    this.vy = new Float32Array(n);
    this.species = new Uint8Array(n);
    this.speciesStart = new Int32Array(k + 1);

    for (var i = 0; i < n; i++) {
      this.px[i] = rng() * this.w;
      this.py[i] = rng() * this.h;
    }
    // Species are laid out in contiguous blocks so rendering can batch one
    // canvas path per colour instead of one per particle.
    for (var s = 0; s <= k; s++) this.speciesStart[s] = Math.round(s * n / k);
    for (var t = 0; t < k; t++) {
      for (var j = this.speciesStart[t]; j < this.speciesStart[t + 1]; j++) {
        this.species[j] = t;
      }
    }

    if (!this.matrix || this.matrix.length !== k * k) this.randomizeMatrix();
    return this;
  };

  /* Keep the torus the same shape as the viewport. Positions are rescaled so
     a resize does not teleport anything. */
  ParticleLife.prototype.setAspect = function (aspect) {
    var w = clamp(aspect, 0.25, 8);
    if (Math.abs(w - this.w) < 1e-6) return;
    var s = w / this.w;
    for (var i = 0; i < this.px.length; i++) this.px[i] *= s;
    this.w = w;
  };

  ParticleLife.prototype._ensureGrid = function () {
    var g = this._grid;
    var r = Math.max(this.config.rMax, 1e-4);
    // >= 3 cells per axis, otherwise the 3x3 stencil wraps onto itself and
    // would count the same neighbour more than once.
    var gx = clamp(Math.floor(this.w / r), 3, 512) | 0;
    var gy = clamp(Math.floor(this.h / r), 3, 512) | 0;
    var n = this.n;

    if (g.gx === gx && g.gy === gy && g.order.length === n) return g;

    var cells = gx * gy;
    g.gx = gx;
    g.gy = gy;
    g.cells = cells;
    g.cellStart = new Int32Array(cells + 1);
    g.cursor = new Int32Array(cells);
    g.order = new Int32Array(n);
    g.cellOf = new Int32Array(n);
    return g;
  };

  ParticleLife.prototype.step = function (dtIn) {
    var cfg = this.config;
    var dt = dtIn === undefined ? cfg.dt : dtIn;
    var n = this.n, k = this.k;
    var px = this.px, py = this.py, vx = this.vx, vy = this.vy;
    var sp = this.species, m = this.matrix;
    var w = this.w, h = this.h;
    var halfW = w * 0.5, halfH = h * 0.5;

    var rMax = cfg.rMax, rMax2 = rMax * rMax;
    var beta = clamp(cfg.beta, 0.01, 0.95);
    var invBeta = 1 / beta;
    var invSpan = 1 / (1 - beta);

    // ---- bin particles into the uniform grid (counting sort, no allocation)
    var g = this._ensureGrid();
    var gx = g.gx, gy = g.gy, cells = g.cells;
    var cellStart = g.cellStart, cursor = g.cursor, order = g.order, cellOf = g.cellOf;
    var sx = gx / w, sy = gy / h;
    var i, c;

    cellStart.fill(0);
    for (i = 0; i < n; i++) {
      var cx = (px[i] * sx) | 0;
      var cy = (py[i] * sy) | 0;
      if (cx < 0) cx = 0; else if (cx >= gx) cx = gx - 1;
      if (cy < 0) cy = 0; else if (cy >= gy) cy = gy - 1;
      c = cy * gx + cx;
      cellOf[i] = c;
      cellStart[c + 1]++;
    }
    for (c = 0; c < cells; c++) cellStart[c + 1] += cellStart[c];
    cursor.set(cellStart.subarray(0, cells));
    for (i = 0; i < n; i++) order[cursor[cellOf[i]]++] = i;

    // ---- accumulate forces, update velocities
    var accel = rMax * cfg.forceScale;
    var damping = Math.pow(0.5, dt / Math.max(cfg.halfLife, 1e-4));
    var vMax = cfg.maxSpeed, vMax2 = vMax * vMax;
    var ptr = this.pointer;
    var ptrR2 = ptr.radius * ptr.radius;

    for (var ii = 0; ii < n; ii++) {
      i = order[ii];                       // grid order: better cache locality
      var xi = px[i], yi = py[i];
      var row = sp[i] * k;
      var ci = cellOf[i];
      var bx = ci % gx, by = (ci / gx) | 0;
      var fx = 0, fy = 0;

      for (var oy = -1; oy <= 1; oy++) {
        var ny = by + oy;
        if (ny < 0) ny += gy; else if (ny >= gy) ny -= gy;
        var base = ny * gx;

        for (var ox = -1; ox <= 1; ox++) {
          var nx = bx + ox;
          if (nx < 0) nx += gx; else if (nx >= gx) nx -= gx;
          var cc = base + nx;

          for (var jj = cellStart[cc], end = cellStart[cc + 1]; jj < end; jj++) {
            var j = order[jj];
            if (j === i) continue;

            var dx = px[j] - xi;
            var dy = py[j] - yi;
            if (dx > halfW) dx -= w; else if (dx < -halfW) dx += w;
            if (dy > halfH) dy -= h; else if (dy < -halfH) dy += h;

            var d2 = dx * dx + dy * dy;
            if (d2 >= rMax2 || d2 === 0) continue;

            var d = Math.sqrt(d2);
            var r = d / rMax;
            var f;
            if (r < beta) {
              f = r * invBeta - 1;                              // core repulsion
            } else {
              f = m[row + sp[j]] * (1 - Math.abs(2 * r - 1 - beta) * invSpan);
            }
            var inv = f / d;
            fx += dx * inv;
            fy += dy * inv;
          }
        }
      }

      if (ptr.active) {
        var qx = xi - ptr.x, qy = yi - ptr.y;
        if (qx > halfW) qx -= w; else if (qx < -halfW) qx += w;
        if (qy > halfH) qy -= h; else if (qy < -halfH) qy += h;
        var q2 = qx * qx + qy * qy;
        if (q2 > 1e-9 && q2 < ptrR2) {
          var qd = Math.sqrt(q2);
          var gq = ptr.strength * (1 - qd / ptr.radius) / qd;
          fx += qx * gq;
          fy += qy * gq;
        }
      }

      var nvx = vx[i] * damping + fx * accel * dt;
      var nvy = vy[i] * damping + fy * accel * dt;
      var s2 = nvx * nvx + nvy * nvy;
      if (!(s2 <= vMax2)) {                     // also catches NaN / Infinity
        if (s2 > 0 && s2 < Infinity) {
          var sc = vMax / Math.sqrt(s2);
          nvx *= sc; nvy *= sc;
        } else {
          nvx = 0; nvy = 0;
        }
      }
      vx[i] = nvx;
      vy[i] = nvy;
    }

    // ---- drift, wrapped onto the torus
    for (i = 0; i < n; i++) {
      var x = px[i] + vx[i] * dt;
      var y = py[i] + vy[i] * dt;
      x -= Math.floor(x / w) * w;
      y -= Math.floor(y / h) * h;
      if (!(x >= 0 && x <= w)) x = halfW;
      if (!(y >= 0 && y <= h)) y = halfH;
      px[i] = x;
      py[i] = y;
    }

    return this;
  };

  // -------------------------------------------------------------------------
  // Matrix <-> URL encoding (int8 quantisation, base64url)
  // -------------------------------------------------------------------------

  function encodeMatrix(m) {
    var s = '';
    for (var i = 0; i < m.length; i++) {
      s += String.fromCharCode((Math.round(clamp(m[i], -1, 1) * 127) + 128) & 255);
    }
    return btoa(s).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
  }

  function decodeMatrix(str, k) {
    try {
      var b = atob(str.replace(/-/g, '+').replace(/_/g, '/'));
      if (b.length !== k * k) return null;
      var m = new Float32Array(k * k);
      for (var i = 0; i < b.length; i++) m[i] = (b.charCodeAt(i) - 128) / 127;
      return m;
    } catch (e) {
      return null;
    }
  }

  // -------------------------------------------------------------------------
  // UI
  // -------------------------------------------------------------------------

  var CONTROLS = [
    { key: 'species', label: 'Species', min: 2, max: 8, step: 1, rebuild: 'matrix', fmt: function (v) { return v.toFixed(0); } },
    { key: 'count', label: 'Particles', min: 200, max: 4000, step: 100, rebuild: 'particles', fmt: function (v) { return v.toFixed(0); } },
    { key: 'rMax', label: 'Cutoff r_max', min: 0.02, max: 0.16, step: 0.005, fmt: function (v) { return v.toFixed(3); } },
    { key: 'beta', label: 'Core radius beta', min: 0.05, max: 0.80, step: 0.01, fmt: function (v) { return v.toFixed(2); } },
    { key: 'forceScale', label: 'Force gain', min: 0, max: 30, step: 0.5, fmt: function (v) { return v.toFixed(1); } },
    { key: 'halfLife', label: 'Velocity half-life', min: 0.005, max: 0.200, step: 0.005, fmt: function (v) { return v.toFixed(3) + ' s'; } },
    { key: 'symmetry', label: 'Matrix symmetry', min: 0, max: 1, step: 0.05, rebuild: 'matrix', fmt: function (v) { return v.toFixed(2); } },
    { key: 'trails', label: 'Trails', min: 0, max: 0.90, step: 0.05, render: true, fmt: function (v) { return v.toFixed(2); } }
  ];

  var TEMPLATE = [
    '<div class="pl-stage">',
    '  <canvas class="pl-canvas" role="img" aria-label="Particle life simulation viewport"></canvas>',
    '  <div class="pl-hud">',
    '    <span><b data-hud="n">0</b> particles</span>',
    '    <span><b data-hud="k">0</b> species</span>',
    '    <span><b data-hud="fps">--</b> fps</span>',
    '  </div>',
    '  <p class="pl-hint">drag to stir &middot; shift-drag to pull</p>',
    '  <div class="pl-overlay" hidden><button type="button" data-act="start">Start simulation</button></div>',
    '</div>',
    '<div class="pl-bar">',
    '  <button type="button" data-act="toggle">Pause</button>',
    '  <button type="button" data-act="step">Step</button>',
    '  <button type="button" data-act="reseed">Reseed</button>',
    '  <button type="button" data-act="random">Randomize matrix</button>',
    '  <span class="pl-spacer"></span>',
    '  <button type="button" data-act="link">Copy link</button>',
    '  <span class="pl-status" role="status"></span>',
    '</div>',
    '<div class="pl-body">',
    '  <div class="pl-sliders"></div>',
    '  <div class="pl-matrix">',
    '    <h4>Interaction matrix A</h4>',
    '    <div class="pl-table"></div>',
    '    <p class="pl-legend">row feels column &middot; drag a cell to edit<br>green attracts &middot; red repels</p>',
    '  </div>',
    '</div>'
  ].join('\n');

  function mount(root) {
    if (root.dataset.plMounted) return;
    root.dataset.plMounted = '1';
    root.classList.add('pl');
    root.innerHTML = TEMPLATE;

    var canvas = root.querySelector('.pl-canvas');
    var ctx = canvas.getContext('2d', { alpha: false });
    if (!ctx) return;

    var engine = new ParticleLife();
    var view = { trails: 0.25 };
    var running = true;
    var needsClear = true;
    var fps = 0;
    var last = 0;
    var acc = 0;
    var raf = 0;
    var frames = 0;
    var onScreen = true;
    var statusTimer = 0;

    var hud = {
      n: root.querySelector('[data-hud="n"]'),
      k: root.querySelector('[data-hud="k"]'),
      fps: root.querySelector('[data-hud="fps"]')
    };
    var statusEl = root.querySelector('.pl-status');
    var overlay = root.querySelector('.pl-overlay');
    var table = root.querySelector('.pl-table');
    var toggleBtn = root.querySelector('[data-act="toggle"]');

    // ---- state <-> URL -----------------------------------------------------

    function readHash() {
      var raw = (location.hash || '').replace(/^#/, '');
      if (raw.indexOf('pl=1') === -1) return;
      var p;
      try { p = new URLSearchParams(raw); } catch (e) { return; }

      var num = function (key, target, obj, lo, hi) {
        if (!p.has(key)) return;
        var v = parseFloat(p.get(key));
        if (isFinite(v)) obj[target] = clamp(v, lo, hi);
      };
      num('k', 'species', engine.config, 2, 8);
      num('n', 'count', engine.config, 200, 4000);
      num('r', 'rMax', engine.config, 0.02, 0.16);
      num('b', 'beta', engine.config, 0.05, 0.80);
      num('f', 'forceScale', engine.config, 0, 30);
      num('hl', 'halfLife', engine.config, 0.005, 0.200);
      num('sy', 'symmetry', engine.config, 0, 1);
      num('t', 'trails', view, 0, 0.90);
      engine.config.species = Math.round(engine.config.species);
      engine.config.count = Math.round(engine.config.count);

      if (p.has('s')) {
        var seed = parseInt(p.get('s'), 10);
        if (isFinite(seed)) engine.config.seed = seed >>> 0;
      }
      engine.reset();
      var m = p.has('m') ? decodeMatrix(p.get('m'), engine.k) : null;
      if (m) engine.matrix = m; else engine.randomizeMatrix();
    }

    function buildURL() {
      var c = engine.config;
      var p = new URLSearchParams();
      p.set('pl', '1');
      p.set('k', String(c.species));
      p.set('n', String(c.count));
      p.set('r', c.rMax.toFixed(3));
      p.set('b', c.beta.toFixed(2));
      p.set('f', c.forceScale.toFixed(1));
      p.set('hl', c.halfLife.toFixed(3));
      p.set('sy', c.symmetry.toFixed(2));
      p.set('t', view.trails.toFixed(2));
      p.set('s', String(c.seed));
      p.set('m', encodeMatrix(engine.matrix));
      return location.origin + location.pathname + location.search + '#' + p.toString();
    }

    function say(msg) {
      statusEl.textContent = msg;
      clearTimeout(statusTimer);
      statusTimer = setTimeout(function () { statusEl.textContent = ''; }, 2400);
    }

    // ---- sliders -----------------------------------------------------------

    var sliderHost = root.querySelector('.pl-sliders');
    var readouts = {};

    CONTROLS.forEach(function (ctl) {
      var wrap = document.createElement('div');
      wrap.className = 'pl-field';

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
      var id = 'pl-' + ctl.key + '-' + Math.random().toString(36).slice(2, 7);
      input.id = id;
      label.htmlFor = id;

      wrap.appendChild(label);
      wrap.appendChild(input);
      sliderHost.appendChild(wrap);

      readouts[ctl.key] = { input: input, value: value, ctl: ctl };

      input.addEventListener('input', function () {
        var v = parseFloat(input.value);
        if (!isFinite(v)) return;
        if (ctl.key === 'trails') {
          view.trails = v;
          needsClear = true;
        } else {
          engine.config[ctl.key] = (ctl.step === 1) ? Math.round(v) : v;
        }
        value.textContent = ctl.fmt(v);

        if (ctl.rebuild === 'particles') {
          engine.reset();
          needsClear = true;
        } else if (ctl.rebuild === 'matrix') {
          engine.randomizeMatrix();
          if (ctl.key === 'species') engine.reset();
          buildMatrix();
          needsClear = true;
        }
        syncHUD();
      });
    });

    function syncSliders() {
      CONTROLS.forEach(function (ctl) {
        var r = readouts[ctl.key];
        var v = ctl.key === 'trails' ? view.trails : engine.config[ctl.key];
        r.input.value = String(v);
        r.value.textContent = ctl.fmt(v);
      });
    }

    // ---- interaction matrix editor ----------------------------------------

    var cells = [];

    function cellColor(v) {
      var a = 0.10 + 0.78 * Math.min(1, Math.abs(v));
      return v >= 0
        ? 'rgba(64, 200, 160, ' + a.toFixed(3) + ')'
        : 'rgba(230, 88, 110, ' + a.toFixed(3) + ')';
    }

    function paintCell(cell) {
      var v = engine.matrix[cell._i * engine.k + cell._j];
      cell.style.background = cellColor(v);
      cell.setAttribute('aria-label',
        'Species ' + (cell._i + 1) + ' feels species ' + (cell._j + 1) +
        ': ' + v.toFixed(2) + '. Use arrow keys to adjust.');
    }

    function setCell(cell, v) {
      engine.matrix[cell._i * engine.k + cell._j] = clamp(v, -1, 1);
      paintCell(cell);
    }

    function buildMatrix() {
      var k = engine.k;
      table.innerHTML = '';
      table.style.gridTemplateColumns = 'repeat(' + (k + 1) + ', 26px)';
      cells = [];

      var corner = document.createElement('div');
      corner.className = 'pl-key';
      table.appendChild(corner);

      for (var j = 0; j < k; j++) {
        var head = document.createElement('div');
        head.className = 'pl-key';
        head.style.color = PALETTE[j % PALETTE.length];
        head.innerHTML = '<i></i>';
        table.appendChild(head);
      }

      for (var i = 0; i < k; i++) {
        var side = document.createElement('div');
        side.className = 'pl-key';
        side.style.color = PALETTE[i % PALETTE.length];
        side.innerHTML = '<i></i>';
        table.appendChild(side);

        for (var c = 0; c < k; c++) {
          var cell = document.createElement('button');
          cell.type = 'button';
          cell.className = 'pl-cell';
          cell._i = i;
          cell._j = c;
          paintCell(cell);
          table.appendChild(cell);
          cells.push(cell);
        }
      }
    }

    var drag = null;

    table.addEventListener('pointerdown', function (ev) {
      var cell = ev.target.closest ? ev.target.closest('.pl-cell') : null;
      if (!cell) return;
      ev.preventDefault();                    // suppress text selection on drag
      cells.forEach(function (c) { c.removeAttribute('data-selected'); });
      cell.dataset.selected = 'true';
      cell.focus();                           // preventDefault ate the focus

      drag = { cell: cell, y: ev.clientY, start: engine.matrix[cell._i * engine.k + cell._j] };
      try { table.setPointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    });

    table.addEventListener('pointermove', function (ev) {
      if (!drag) return;
      setCell(drag.cell, drag.start + (drag.y - ev.clientY) / 70);
    });

    function endDrag(ev) {
      if (!drag) return;
      drag = null;
      try { table.releasePointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    }
    table.addEventListener('pointerup', endDrag);
    table.addEventListener('pointercancel', endDrag);

    table.addEventListener('keydown', function (ev) {
      var cell = ev.target.closest ? ev.target.closest('.pl-cell') : null;
      if (!cell) return;
      var v = engine.matrix[cell._i * engine.k + cell._j];
      var d = 0;
      if (ev.key === 'ArrowUp' || ev.key === 'ArrowRight') d = 0.05;
      else if (ev.key === 'ArrowDown' || ev.key === 'ArrowLeft') d = -0.05;
      else if (ev.key === 'Home') { setCell(cell, 0); ev.preventDefault(); return; }
      else return;
      ev.preventDefault();
      setCell(cell, v + d);
    });

    // ---- pointer stirring --------------------------------------------------

    function toWorld(ev) {
      var rect = canvas.getBoundingClientRect();
      var s = engine.h / rect.height;                 // world units per CSS px
      return {
        x: clamp((ev.clientX - rect.left) * s, 0, engine.w),
        y: clamp((ev.clientY - rect.top) * s, 0, engine.h)
      };
    }

    canvas.addEventListener('pointerdown', function (ev) {
      if (ev.pointerType === 'touch') return;         // leave touch to scrolling
      var p = toWorld(ev);
      engine.pointer.active = true;
      engine.pointer.x = p.x;
      engine.pointer.y = p.y;
      engine.pointer.strength = ev.shiftKey ? -3 : 3;
      try { canvas.setPointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    });

    canvas.addEventListener('pointermove', function (ev) {
      if (!engine.pointer.active) return;
      var p = toWorld(ev);
      engine.pointer.x = p.x;
      engine.pointer.y = p.y;
      engine.pointer.strength = ev.shiftKey ? -3 : 3;
    });

    function releasePointer(ev) {
      engine.pointer.active = false;
      try { canvas.releasePointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    }
    canvas.addEventListener('pointerup', releasePointer);
    canvas.addEventListener('pointercancel', releasePointer);
    canvas.addEventListener('pointerleave', releasePointer);

    // ---- toolbar -----------------------------------------------------------

    root.addEventListener('click', function (ev) {
      var btn = ev.target.closest ? ev.target.closest('[data-act]') : null;
      if (!btn) return;
      var act = btn.dataset.act;

      if (act === 'toggle') {
        setRunning(!running);
      } else if (act === 'start') {
        overlay.hidden = true;
        setRunning(true);
      } else if (act === 'step') {
        setRunning(false);
        engine.step();
        render();
      } else if (act === 'reseed') {
        engine.reset((Math.random() * 0xffffffff) >>> 0);
        needsClear = true;
        render();
      } else if (act === 'random') {
        engine.randomizeMatrix((Math.random() * 0xffffffff) >>> 0);
        buildMatrix();
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

    function setRunning(next) {
      running = next;
      toggleBtn.textContent = running ? 'Pause' : 'Play';
      toggleBtn.setAttribute('aria-pressed', String(!running));
      if (running) { last = 0; acc = 0; }
    }

    // ---- rendering ---------------------------------------------------------

    function resize() {
      var dpr = Math.min(2, window.devicePixelRatio || 1);
      var rect = canvas.getBoundingClientRect();
      var w = Math.max(1, Math.round(rect.width * dpr));
      var h = Math.max(1, Math.round(rect.height * dpr));
      if (w === canvas.width && h === canvas.height) return;
      canvas.width = w;
      canvas.height = h;
      engine.setAspect(w / h);
      needsClear = true;
    }

    function render() {
      var cw = canvas.width, ch = canvas.height;
      ctx.fillStyle = INK;
      if (needsClear || view.trails < 0.001) {
        ctx.fillRect(0, 0, cw, ch);
        needsClear = false;
      } else {
        ctx.globalAlpha = 1 - view.trails;
        ctx.fillRect(0, 0, cw, ch);
        ctx.globalAlpha = 1;
      }

      var scale = ch;                       // world height (= 1) fills the canvas
      var n = engine.n;
      var radius = Math.max(1.1, ch * 0.0035 * Math.sqrt(2000 / n));
      var px = engine.px, py = engine.py, start = engine.speciesStart;

      for (var s = 0; s < engine.k; s++) {
        var a = start[s], b = start[s + 1];
        if (b <= a) continue;
        ctx.fillStyle = PALETTE[s % PALETTE.length];
        ctx.beginPath();
        for (var i = a; i < b; i++) {
          var x = px[i] * scale, y = py[i] * scale;
          ctx.moveTo(x + radius, y);
          ctx.arc(x, y, radius, 0, TAU);
        }
        ctx.fill();
      }
    }

    function syncHUD() {
      hud.n.textContent = String(engine.n);
      hud.k.textContent = String(engine.k);
      hud.fps.textContent = fps > 0 ? fps.toFixed(0) : '--';
    }

    // ---- main loop ---------------------------------------------------------

    function frame(now) {
      raf = requestAnimationFrame(frame);
      if (!running || !onScreen) { last = 0; return; }

      resize();

      if (last) {
        var elapsed = (now - last) / 1000;
        if (elapsed > 0) fps += (1 / elapsed - fps) * 0.08;
        acc += Math.min(0.05, elapsed);
      }
      last = now;

      var dt = engine.config.dt;
      var steps = 0;
      while (acc >= dt && steps < 4) { engine.step(dt); acc -= dt; steps++; }
      if (acc > dt * 4) acc = 0;            // never spiral on a slow device

      render();
      if (++frames % 15 === 0) syncHUD();     // ~4 Hz readout
    }

    // ---- boot --------------------------------------------------------------

    readHash();
    buildMatrix();
    syncSliders();
    syncHUD();
    resize();
    render();

    if (typeof ResizeObserver !== 'undefined') {
      new ResizeObserver(function () { resize(); render(); }).observe(canvas);
    } else {
      window.addEventListener('resize', function () { resize(); render(); });
    }

    if (typeof IntersectionObserver !== 'undefined') {
      new IntersectionObserver(function (entries) {
        onScreen = entries[0].isIntersecting;
      }, { threshold: 0.01 }).observe(root);
    }

    document.addEventListener('visibilitychange', function () {
      if (document.hidden) { last = 0; acc = 0; }
    });

    var reduced = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    setRunning(!reduced);
    if (reduced) overlay.hidden = false;

    raf = requestAnimationFrame(frame);
    root._particleLife = { engine: engine, stop: function () { cancelAnimationFrame(raf); } };
  }

  // -------------------------------------------------------------------------

  if (typeof window !== 'undefined') {
    window.ParticleLife = ParticleLife;
  }
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { ParticleLife: ParticleLife, mulberry32: mulberry32 };
  }
  if (typeof document !== 'undefined') {
    var boot = function () {
      var nodes = document.querySelectorAll('[data-sim="particle-life"]');
      for (var i = 0; i < nodes.length; i++) mount(nodes[i]);
    };
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', boot);
    } else {
      boot();
    }
  }
})();
