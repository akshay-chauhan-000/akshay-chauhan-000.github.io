/* ==========================================================================
   Pathfinding on Energy Landscapes
   --------------------------------------------------------------------------
   Shortest paths between two minima of a scalar field E(x, y) sampled on a
   grid. The cost of stepping into a cell is Arrhenius-weighted,

       c(a -> b) = ell_ab * exp(E_b / T)

   with ell_ab the geometric step length (1 or sqrt(2); 8-connected moves)
   and T a temperature. As T -> infinity the optimal route tends to the
   geometrically shortest path; as T -> 0 the exponential makes high-energy
   cells prohibitively expensive and the route hugs valleys and crosses the
   landscape at its lowest saddle — a discrete stand-in for the minimum
   energy path of chemical kinetics. A maze is the degenerate case where
   barriers are infinite: walls are simply impassable and every open cell
   costs the same.

   Search algorithms: Dijkstra, A* (admissible octile heuristic scaled by a
   weight w; w <= 1 preserves optimality), greedy best-first, and BFS. The
   open set is a binary heap with lazy deletion; all per-node state lives in
   preallocated typed arrays.

   No dependencies, no build step. The engine (Field/Search/TERRAINS) is
   DOM-free and exported on window.PathfindingSim for headless use; the UI
   below mounts into [data-sim="pathfinding"].
   ========================================================================== */

(function () {
  'use strict';

  var SQRT2 = Math.SQRT2;
  var INK = '#0b0d12';

  /* 8-connected neighbour offsets and step lengths. */
  var OX = [1, -1, 0, 0, 1, 1, -1, -1];
  var OY = [0, 0, 1, -1, 1, -1, 1, -1];
  var OL = [1, 1, 1, 1, SQRT2, SQRT2, SQRT2, SQRT2];

  /* Energy colour ramp: ink valleys -> indigo -> slate -> copper -> sand. */
  var RAMP = [
    [0.00, 11, 13, 18],
    [0.25, 36, 48, 94],
    [0.50, 86, 98, 142],
    [0.75, 189, 141, 87],
    [1.00, 246, 231, 178]
  ];

  var COL_VISITED = [64, 190, 178, 66];
  var COL_FRONTIER = [80, 216, 255, 235];
  var COL_START = '#5ddb8b';
  var COL_GOAL = '#ff5d73';
  var COL_PATH = '#ffffff';

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

  // -------------------------------------------------------------------------
  // Field: the landscape
  // -------------------------------------------------------------------------

  function Field(gw, gh) {
    this.gw = gw;
    this.gh = gh;
    this.energy = new Float32Array(gw * gh);   // normalised to [0, 1]
    this.wall = new Uint8Array(gw * gh);
    this.maze = false;                          // walls + uniform cost mode
  }

  Field.prototype.idx = function (x, y) { return y * this.gw + x; };

  function normalizeEnergy(f) {
    var e = f.energy;
    var lo = Infinity, hi = -Infinity;
    for (var i = 0; i < e.length; i++) {
      if (e[i] < lo) lo = e[i];
      if (e[i] > hi) hi = e[i];
    }
    var span = hi - lo;
    if (span < 1e-12) { e.fill(0); return; }
    for (var j = 0; j < e.length; j++) e[j] = (e[j] - lo) / span;
  }

  /* Deepest cell, then the deepest cell at least ~half the domain away. */
  function deepestPair(f) {
    var e = f.energy, gw = f.gw, gh = f.gh;
    var s = 0, i;
    for (i = 1; i < e.length; i++) if (e[i] < e[s]) s = i;
    var sx = s % gw, sy = (s / gw) | 0;
    var minD = 0.45 * Math.sqrt(gw * gw + gh * gh);
    var g = -1, bestE = Infinity, bestD = -1, far = s;
    for (i = 0; i < e.length; i++) {
      var dx = (i % gw) - sx, dy = ((i / gw) | 0) - sy;
      var d = Math.sqrt(dx * dx + dy * dy);
      if (d > bestD) { bestD = d; far = i; }
      if (d >= minD && e[i] < bestE) { bestE = e[i]; g = i; }
    }
    if (g < 0) g = far;
    return { start: [sx, sy], goal: [g % gw, (g / gw) | 0] };
  }

  var TERRAINS = {
    maze: {
      label: 'Maze',
      make: function (f, rng) {
        var gw = f.gw, gh = f.gh;
        f.maze = true;
        f.wall.fill(1);
        // Recursive backtracker on the odd sublattice.
        var stack = [[1, 1]];
        f.wall[f.idx(1, 1)] = 0;
        var DIRS = [[2, 0], [-2, 0], [0, 2], [0, -2]];
        while (stack.length) {
          var top = stack[stack.length - 1];
          var x = top[0], y = top[1];
          var options = [];
          for (var d = 0; d < 4; d++) {
            var nx = x + DIRS[d][0], ny = y + DIRS[d][1];
            if (nx > 0 && ny > 0 && nx < gw - 1 && ny < gh - 1 &&
                f.wall[f.idx(nx, ny)]) {
              options.push(d);
            }
          }
          if (!options.length) { stack.pop(); continue; }
          var pick = DIRS[options[(rng() * options.length) | 0]];
          var mx = x + pick[0] / 2, my = y + pick[1] / 2;
          var tx = x + pick[0], ty = y + pick[1];
          f.wall[f.idx(mx, my)] = 0;
          f.wall[f.idx(tx, ty)] = 0;
          stack.push([tx, ty]);
        }
        for (var i = 0; i < f.energy.length; i++) f.energy[i] = f.wall[i];
        return { start: [1, 1], goal: [gw - 2, gh - 2] };
      }
    },

    hills: {
      label: 'Rolling hills',
      make: function (f, rng) {
        var gw = f.gw, gh = f.gh;
        f.maze = false;
        f.wall.fill(0);
        var K = 14, cx = [], cy = [], amp = [], sx = [], sy = [];
        for (var k = 0; k < K; k++) {
          cx.push(0.05 + 0.9 * rng());
          cy.push(0.05 + 0.9 * rng());
          amp.push((rng() < 0.5 ? -1 : 1) * (0.35 + 0.65 * rng()));
          sx.push(0.05 + 0.13 * rng());
          sy.push(0.05 + 0.13 * rng());
        }
        for (var y = 0; y < gh; y++) {
          var v = y / (gh - 1);
          for (var x = 0; x < gw; x++) {
            var u = x / (gw - 1), e = 0;
            for (var q = 0; q < K; q++) {
              var du = (u - cx[q]) / sx[q], dv = (v - cy[q]) / sy[q];
              e += amp[q] * Math.exp(-(du * du + dv * dv));
            }
            f.energy[f.idx(x, y)] = e;
          }
        }
        normalizeEnergy(f);
        return deepestPair(f);
      }
    },

    doublewell: {
      label: 'Double well',
      make: function (f) {
        // Two basins separated by a ridge pierced by two saddles of
        // different heights — the temperature slider picks between them.
        var gw = f.gw, gh = f.gh;
        f.maze = false;
        f.wall.fill(0);
        function sq(z) { return z * z; }
        for (var y = 0; y < gh; y++) {
          var v = y / (gh - 1);
          for (var x = 0; x < gw; x++) {
            var u = x / (gw - 1);
            var ridge = Math.exp(-sq((u - 0.5) / 0.06));
            var slot = Math.exp(-sq((u - 0.5) / 0.10));
            var e = ridge
              - 0.75 * slot * Math.exp(-sq((v - 0.28) / 0.07))
              - 0.45 * slot * Math.exp(-sq((v - 0.72) / 0.07))
              - 0.55 * Math.exp(-(sq((u - 0.20) / 0.10) + sq((v - 0.50) / 0.16)))
              - 0.55 * Math.exp(-(sq((u - 0.80) / 0.10) + sq((v - 0.50) / 0.16)));
            f.energy[f.idx(x, y)] = e;
          }
        }
        normalizeEnergy(f);
        return {
          start: [Math.round(0.20 * (gw - 1)), Math.round(0.5 * (gh - 1))],
          goal: [Math.round(0.80 * (gw - 1)), Math.round(0.5 * (gh - 1))]
        };
      }
    },

    muller: {
      label: 'Muller-Brown',
      make: function (f) {
        // The standard test surface of reaction-path methods
        // (Muller & Brown, Theor. Chim. Acta 53, 75 (1979)).
        var A = [-200, -100, -170, 15];
        var a = [-1, -1, -6.5, 0.7];
        var b = [0, 0, 11, 0.6];
        var c = [-10, -10, -6.5, 0.7];
        var X0 = [1, 0, -0.5, -1];
        var Y0 = [0, 0.5, 1.5, 1];
        var xmin = -1.7, xmax = 1.3, ymin = -0.4, ymax = 2.1;
        var gw = f.gw, gh = f.gh;
        f.maze = false;
        f.wall.fill(0);
        for (var gy = 0; gy < gh; gy++) {
          var y = ymax - (gy / (gh - 1)) * (ymax - ymin);   // canvas y is down
          for (var gx = 0; gx < gw; gx++) {
            var x = xmin + (gx / (gw - 1)) * (xmax - xmin);
            var e = 0;
            for (var k = 0; k < 4; k++) {
              var dx = x - X0[k], dy = y - Y0[k];
              e += A[k] * Math.exp(a[k] * dx * dx + b[k] * dx * dy + c[k] * dy * dy);
            }
            f.energy[f.idx(gx, gy)] = clamp(e, -150, 150);
          }
        }
        normalizeEnergy(f);
        function map(px, py) {
          return [
            Math.round((px - xmin) / (xmax - xmin) * (gw - 1)),
            Math.round((ymax - py) / (ymax - ymin) * (gh - 1))
          ];
        }
        return { start: map(-0.558, 1.442), goal: map(0.623, 0.028) };
      }
    },

    blank: {
      label: 'Blank canvas',
      make: function (f) {
        f.maze = false;
        f.wall.fill(0);
        f.energy.fill(0);
        return {
          start: [Math.round(0.12 * (f.gw - 1)), (f.gh / 2) | 0],
          goal: [Math.round(0.88 * (f.gw - 1)), (f.gh / 2) | 0]
        };
      }
    }
  };

  // -------------------------------------------------------------------------
  // Search
  // -------------------------------------------------------------------------

  function Search(field, opts) {
    var n = field.gw * field.gh;
    this.field = field;
    this.algo = opts.algo;                     // astar | dijkstra | greedy | bfs
    this.w = opts.w === undefined ? 1 : opts.w;
    this.start = opts.start;
    this.goal = opts.goal;
    this.expanded = 0;
    this.status = 'running';
    this.path = null;
    this.barrier = 0;

    this.g = new Float64Array(n);
    this.g.fill(Infinity);
    this.state = new Uint8Array(n);            // 0 unseen, 1 open, 2 closed
    this.prev = new Int32Array(n);
    this.prev.fill(-1);
    this.key = new Float64Array(n);

    // Arrhenius cost multiplier per cell; 1 everywhere in maze mode.
    this.mult = new Float64Array(n);
    if (field.maze) {
      this.mult.fill(1);
    } else {
      var beta = 1 / Math.max(opts.T || 0.25, 1e-3);
      for (var i = 0; i < n; i++) this.mult[i] = Math.exp(beta * field.energy[i]);
    }

    this.heap = [];
    this.queue = new Int32Array(n);            // BFS visits each node once
    this.qh = 0;
    this.qt = 0;

    var s = this.start;
    this.g[s] = 0;
    this.state[s] = 1;
    if (this.algo === 'bfs') {
      this.queue[this.qt++] = s;
    } else {
      this.key[s] = this.algo === 'dijkstra' ? 0 : this.h(s);
      this.push(s);
    }
  }

  /* Octile distance — admissible for 8-connected moves with unit multiplier. */
  Search.prototype.h = function (id) {
    var gw = this.field.gw;
    var dx = Math.abs((id % gw) - (this.goal % gw));
    var dy = Math.abs(((id / gw) | 0) - ((this.goal / gw) | 0));
    return dx > dy ? dx + (SQRT2 - 1) * dy : dy + (SQRT2 - 1) * dx;
  };

  Search.prototype.push = function (id) {
    var hp = this.heap, k = this.key;
    hp.push(id);
    var ci = hp.length - 1;
    while (ci > 0) {
      var p = (ci - 1) >> 1;
      if (k[hp[p]] <= k[hp[ci]]) break;
      var t = hp[p]; hp[p] = hp[ci]; hp[ci] = t;
      ci = p;
    }
  };

  Search.prototype.pop = function () {
    var hp = this.heap, k = this.key;
    var top = hp[0];
    var last = hp.pop();
    if (hp.length) {
      hp[0] = last;
      var ci = 0, len = hp.length;
      for (;;) {
        var l = 2 * ci + 1, r = l + 1, m = ci;
        if (l < len && k[hp[l]] < k[hp[m]]) m = l;
        if (r < len && k[hp[r]] < k[hp[m]]) m = r;
        if (m === ci) break;
        var t = hp[m]; hp[m] = hp[ci]; hp[ci] = t;
        ci = m;
      }
    }
    return top;
  };

  /* Expand one node. Returns 'running' | 'done' | 'blocked'. */
  Search.prototype.expand = function () {
    if (this.status !== 'running') return this.status;

    var id = -1;
    if (this.algo === 'bfs') {
      if (this.qh >= this.qt) { this.status = 'blocked'; return this.status; }
      id = this.queue[this.qh++];
    } else {
      while (this.heap.length) {              // lazy deletion of stale entries
        var t = this.pop();
        if (this.state[t] !== 2) { id = t; break; }
      }
      if (id < 0) { this.status = 'blocked'; return this.status; }
    }

    this.state[id] = 2;
    this.expanded++;
    if (id === this.goal) {
      this.status = 'done';
      this.trace();
      return this.status;
    }

    var f = this.field, gw = f.gw, gh = f.gh, wall = f.wall;
    var x = id % gw, y = (id / gw) | 0;

    for (var o = 0; o < 8; o++) {
      var nx = x + OX[o], ny = y + OY[o];
      if (nx < 0 || ny < 0 || nx >= gw || ny >= gh) continue;
      var nid = ny * gw + nx;
      if (wall[nid]) continue;
      // No corner cutting: a diagonal move needs both orthogonals open.
      if (OX[o] !== 0 && OY[o] !== 0 &&
          (wall[y * gw + nx] || wall[ny * gw + x])) continue;

      if (this.algo === 'bfs') {
        if (this.state[nid] === 0) {
          this.state[nid] = 1;
          this.prev[nid] = id;
          this.g[nid] = this.g[id] + OL[o] * this.mult[nid];   // stats only
          this.queue[this.qt++] = nid;
        }
        continue;
      }

      if (this.state[nid] === 2) continue;
      var ng = this.g[id] + OL[o] * this.mult[nid];
      if (ng < this.g[nid] - 1e-12) {
        this.g[nid] = ng;
        this.prev[nid] = id;
        this.key[nid] =
          this.algo === 'greedy' ? this.h(nid) :
          this.algo === 'astar' ? ng + this.w * this.h(nid) : ng;
        this.state[nid] = 1;
        this.push(nid);
      }
    }
    return this.status;
  };

  Search.prototype.run = function () {
    while (this.status === 'running') this.expand();
    return this.status;
  };

  Search.prototype.trace = function () {
    var path = [];
    var barrier = 0;
    var id = this.goal;
    while (id >= 0) {
      path.push(id);
      if (this.field.energy[id] > barrier) barrier = this.field.energy[id];
      id = this.prev[id];
    }
    path.reverse();
    this.path = path;
    this.barrier = barrier;
  };

  // -------------------------------------------------------------------------
  // UI
  // -------------------------------------------------------------------------

  var ALGOS = [
    ['astar', 'A*'],
    ['dijkstra', 'Dijkstra'],
    ['greedy', 'Greedy best-first'],
    ['bfs', 'Breadth-first (BFS)']
  ];

  var SLIDERS = [
    { key: 'T', label: 'Temperature T', min: 0.05, max: 1.5, step: 0.05, fmt: function (v) { return v.toFixed(2); } },
    { key: 'w', label: 'Heuristic weight w', min: 0, max: 1.5, step: 0.05, fmt: function (v) { return v.toFixed(2); } },
    { key: 'speed', label: 'Speed (cells/frame)', min: 1, max: 400, step: 1, fmt: function (v) { return v.toFixed(0); } },
    { key: 'res', label: 'Grid resolution', min: 41, max: 201, step: 20, fmt: function (v) { return v.toFixed(0); } }
  ];

  var TEMPLATE = [
    '<div class="pf-stage">',
    '  <canvas class="pf-canvas" role="img" aria-label="Pathfinding on an energy landscape"></canvas>',
    '  <div class="pf-hud">',
    '    <span><b data-hud="cells">0</b> cells</span>',
    '    <span><b data-hud="expanded">0</b> expanded</span>',
    '    <span>cost <b data-hud="cost">--</b></span>',
    '    <span>barrier <b data-hud="barrier">--</b></span>',
    '  </div>',
    '  <p class="pf-hint"></p>',
    '</div>',
    '<div class="pf-bar">',
    '  <button type="button" data-act="run">Run</button>',
    '  <button type="button" data-act="step">Step</button>',
    '  <button type="button" data-act="reset">Reset</button>',
    '  <button type="button" data-act="reseed">New landscape</button>',
    '  <span class="pf-spacer"></span>',
    '  <button type="button" data-act="link">Copy link</button>',
    '  <span class="pf-status" role="status"></span>',
    '</div>',
    '<div class="pf-body">',
    '  <div class="pf-controls">',
    '    <div class="pf-field"><label><span>Landscape</span></label><select data-sel="terrain"></select></div>',
    '    <div class="pf-field"><label><span>Algorithm</span></label><select data-sel="algo"></select></div>',
    '  </div>',
    '  <div class="pf-legend">',
    '    <h4>Energy</h4>',
    '    <div class="pf-ramp"></div>',
    '    <div class="pf-ramp-labels"><span>valley</span><span>ridge</span></div>',
    '    <div class="pf-keys">',
    '      <span><i style="background:' + COL_START + '"></i> start</span>',
    '      <span><i style="background:' + COL_GOAL + '"></i> goal</span>',
    '      <span><i style="background:rgb(' + COL_FRONTIER[0] + ',' + COL_FRONTIER[1] + ',' + COL_FRONTIER[2] + ')"></i> frontier</span>',
    '      <span><i style="background:rgba(' + COL_VISITED[0] + ',' + COL_VISITED[1] + ',' + COL_VISITED[2] + ',0.6)"></i> visited</span>',
    '      <span><i style="background:#fff"></i> path</span>',
    '    </div>',
    '  </div>',
    '</div>'
  ].join('\n');

  function mount(root) {
    if (root.dataset.pfMounted) return;
    root.dataset.pfMounted = '1';
    root.classList.add('pf');
    root.innerHTML = TEMPLATE;

    var canvas = root.querySelector('.pf-canvas');
    var ctx = canvas.getContext('2d', { alpha: false });
    if (!ctx) return;

    var conf = {
      terrain: 'maze', algo: 'astar',
      T: 0.25, w: 1, speed: 80, res: 101,
      seed: 20260807
    };

    var field = null;
    var start = [1, 1], goal = [3, 3];
    var search = null;
    var animating = false;
    var terrainCanvas = document.createElement('canvas');
    var overlayCanvas = document.createElement('canvas');
    var overlayData = null;
    var statusTimer = 0;

    var hud = {
      cells: root.querySelector('[data-hud="cells"]'),
      expanded: root.querySelector('[data-hud="expanded"]'),
      cost: root.querySelector('[data-hud="cost"]'),
      barrier: root.querySelector('[data-hud="barrier"]')
    };
    var statusEl = root.querySelector('.pf-status');
    var hintEl = root.querySelector('.pf-hint');
    var runBtn = root.querySelector('[data-act="run"]');

    // ---- selects and sliders ----------------------------------------------

    var selTerrain = root.querySelector('[data-sel="terrain"]');
    var selAlgo = root.querySelector('[data-sel="algo"]');
    Object.keys(TERRAINS).forEach(function (k) {
      var o = document.createElement('option');
      o.value = k;
      o.textContent = TERRAINS[k].label;
      selTerrain.appendChild(o);
    });
    ALGOS.forEach(function (p) {
      var o = document.createElement('option');
      o.value = p[0];
      o.textContent = p[1];
      selAlgo.appendChild(o);
    });

    var controlsHost = root.querySelector('.pf-controls');
    var readouts = {};
    SLIDERS.forEach(function (ctl) {
      var wrap = document.createElement('div');
      wrap.className = 'pf-field';
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
        if (ctl.key === 'res') {
          regen();
        } else if (ctl.key !== 'speed') {
          stopAnim();
          instant();
          draw();
        }
      });
    });

    selTerrain.addEventListener('change', function () {
      conf.terrain = selTerrain.value;
      regen();
    });
    selAlgo.addEventListener('change', function () {
      conf.algo = selAlgo.value;
      stopAnim();
      instant();
      draw();
    });

    function syncControls() {
      selTerrain.value = conf.terrain;
      selAlgo.value = conf.algo;
      SLIDERS.forEach(function (ctl) {
        var r = readouts[ctl.key];
        r.input.value = String(conf[ctl.key]);
        r.value.textContent = ctl.fmt(conf[ctl.key]);
      });
      readouts.T.input.disabled = conf.terrain === 'maze';
      readouts.w.input.disabled = conf.algo !== 'astar';
      hintEl.textContent = conf.terrain === 'maze'
        ? 'drag markers - paint walls (shift erases)'
        : 'drag markers - paint hills (shift carves)';
    }

    // ---- landscape lifecycle ----------------------------------------------

    function oddify(v) { v = Math.round(v); return v % 2 ? v : v + 1; }

    function regen() {
      stopAnim();
      var rect = canvas.getBoundingClientRect();
      var aspect = rect.width > 0 && rect.height > 0 ? rect.width / rect.height : 1.9;
      var gw = oddify(clamp(conf.res, 21, 301));
      var gh = oddify(clamp(Math.round(gw / aspect), 15, 301));
      field = new Field(gw, gh);
      var def = TERRAINS[conf.terrain].make(field, mulberry32(conf.seed));
      start = def.start.slice();
      goal = def.goal.slice();
      terrainCanvas.width = gw;
      terrainCanvas.height = gh;
      overlayCanvas.width = gw;
      overlayCanvas.height = gh;
      overlayData = overlayCanvas.getContext('2d').createImageData(gw, gh);
      bake();
      search = null;
      instant();
      syncControls();
      draw();
    }

    function bake() {
      var gw = field.gw, gh = field.gh;
      var tctx = terrainCanvas.getContext('2d');
      var img = tctx.createImageData(gw, gh);
      var px = img.data, rgb = [0, 0, 0];
      for (var i = 0; i < gw * gh; i++) {
        rampColor(field.energy[i], rgb);
        px[4 * i] = rgb[0];
        px[4 * i + 1] = rgb[1];
        px[4 * i + 2] = rgb[2];
        px[4 * i + 3] = 255;
      }
      tctx.putImageData(img, 0, 0);
    }

    function makeSearch() {
      return new Search(field, {
        algo: conf.algo, w: conf.w, T: conf.T,
        start: field.idx(start[0], start[1]),
        goal: field.idx(goal[0], goal[1])
      });
    }

    /* Solve synchronously with no frontier movie — used while dragging. */
    function instant() {
      search = makeSearch();
      search.run();
      overlayData.data.fill(0);              // path only, no visited wash
      overlayCanvas.getContext('2d').putImageData(overlayData, 0, 0);
      updateHUD();
    }

    // ---- animation --------------------------------------------------------

    function stopAnim() { animating = false; }

    function rebuildOverlay() {
      var st = search.state, d = overlayData.data;
      for (var i = 0; i < st.length; i++) {
        var o = 4 * i, c = null;
        if (st[i] === 2) c = COL_VISITED;
        else if (st[i] === 1) c = COL_FRONTIER;
        if (c) { d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = c[3]; }
        else d[o + 3] = 0;
      }
      overlayCanvas.getContext('2d').putImageData(overlayData, 0, 0);
    }

    function tick() {
      if (!animating) return;
      var k = Math.max(1, Math.round(conf.speed));
      while (k-- > 0 && search.status === 'running') search.expand();
      rebuildOverlay();
      updateHUD();
      draw();
      if (search.status === 'running') {
        requestAnimationFrame(tick);
      } else {
        animating = false;
        runBtn.textContent = 'Run';
      }
    }

    function startAnim(fresh) {
      if (fresh || !search || search.status !== 'running' || search.expanded === 0) {
        search = makeSearch();
      }
      animating = true;
      runBtn.textContent = 'Pause';
      requestAnimationFrame(tick);
    }

    // ---- rendering --------------------------------------------------------

    function resize() {
      var dpr = Math.min(2, window.devicePixelRatio || 1);
      var rect = canvas.getBoundingClientRect();
      var w = Math.max(1, Math.round(rect.width * dpr));
      var h = Math.max(1, Math.round(rect.height * dpr));
      if (w !== canvas.width || h !== canvas.height) {
        canvas.width = w;
        canvas.height = h;
      }
    }

    function cellCenter(id) {
      var gw = field.gw;
      return [
        ((id % gw) + 0.5) * canvas.width / gw,
        (((id / gw) | 0) + 0.5) * canvas.height / field.gh
      ];
    }

    function drawMarker(cellXY, color) {
      var cw = canvas.width / field.gw;
      var x = (cellXY[0] + 0.5) * cw;
      var y = (cellXY[1] + 0.5) * canvas.height / field.gh;
      var r = Math.max(5, cw * 0.8);
      ctx.beginPath();
      ctx.arc(x, y, r, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.fill();
      ctx.lineWidth = Math.max(1.5, r * 0.22);
      ctx.strokeStyle = INK;
      ctx.stroke();
    }

    function draw() {
      resize();
      var W = canvas.width, H = canvas.height;
      ctx.imageSmoothingEnabled = !field.maze;
      ctx.fillStyle = INK;
      ctx.fillRect(0, 0, W, H);
      ctx.drawImage(terrainCanvas, 0, 0, W, H);
      ctx.drawImage(overlayCanvas, 0, 0, W, H);

      if (search && search.path && search.status === 'done') {
        ctx.beginPath();
        var p = search.path;
        for (var i = 0; i < p.length; i++) {
          var c = cellCenter(p[i]);
          if (i === 0) ctx.moveTo(c[0], c[1]);
          else ctx.lineTo(c[0], c[1]);
        }
        ctx.strokeStyle = COL_PATH;
        ctx.lineWidth = Math.max(2, canvas.width / field.gw * 0.35);
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        ctx.stroke();
      }

      drawMarker(start, COL_START);
      drawMarker(goal, COL_GOAL);
    }

    function updateHUD() {
      hud.cells.textContent = field.gw + 'x' + field.gh;
      hud.expanded.textContent = search ? String(search.expanded) : '0';
      if (search && search.status === 'done') {
        hud.cost.textContent = search.g[search.goal].toFixed(1);
        hud.barrier.textContent = field.maze ? '--' : search.barrier.toFixed(2);
      } else {
        hud.cost.textContent = '--';
        hud.barrier.textContent = '--';
      }
      statusEl.textContent = statusText();
    }

    function statusText() {
      if (!search) return '';
      if (search.status === 'running') {
        return animating ? 'expanding...' : 'paused';
      }
      if (search.status === 'blocked') return 'no path - goal sealed off';
      if (conf.algo === 'dijkstra' || (conf.algo === 'astar' && conf.w <= 1)) {
        return 'optimal path';
      }
      if (conf.algo === 'astar') return 'w > 1: fast, may miss the optimum';
      if (conf.algo === 'greedy') return 'greedy: fast, not optimal';
      return 'fewest steps - BFS ignores the terrain';
    }

    function say(msg) {
      statusEl.textContent = msg;
      clearTimeout(statusTimer);
      statusTimer = setTimeout(updateHUD, 2400);
    }

    // ---- pointer: markers and painting ------------------------------------

    function toCell(ev) {
      var rect = canvas.getBoundingClientRect();
      return [
        clamp(Math.floor((ev.clientX - rect.left) / rect.width * field.gw), 0, field.gw - 1),
        clamp(Math.floor((ev.clientY - rect.top) / rect.height * field.gh), 0, field.gh - 1)
      ];
    }

    var drag = null;   // {mode: 'start'|'goal'|'paint', wallValue, last}

    function nearMarker(ev, cellXY) {
      var rect = canvas.getBoundingClientRect();
      var cw = rect.width / field.gw;
      var mx = (cellXY[0] + 0.5) * cw + rect.left;
      var my = (cellXY[1] + 0.5) * rect.height / field.gh + rect.top;
      return Math.hypot(ev.clientX - mx, ev.clientY - my) < Math.max(14, cw);
    }

    function paintCell(cell, erase) {
      var gw = field.gw, gh = field.gh;
      if (field.maze) {
        var id = field.idx(cell[0], cell[1]);
        if ((cell[0] === start[0] && cell[1] === start[1]) ||
            (cell[0] === goal[0] && cell[1] === goal[1])) return;
        field.wall[id] = drag.wallValue;
        field.energy[id] = drag.wallValue;
      } else {
        var r = Math.max(2, gw * 0.045);
        var amp = (erase ? -1 : 1) * 0.16;
        var x0 = Math.max(0, Math.round(cell[0] - 3 * r));
        var x1 = Math.min(gw - 1, Math.round(cell[0] + 3 * r));
        var y0 = Math.max(0, Math.round(cell[1] - 3 * r));
        var y1 = Math.min(gh - 1, Math.round(cell[1] + 3 * r));
        for (var y = y0; y <= y1; y++) {
          for (var x = x0; x <= x1; x++) {
            var dx = (x - cell[0]) / r, dy = (y - cell[1]) / r;
            var i = field.idx(x, y);
            field.energy[i] = clamp(field.energy[i] + amp * Math.exp(-(dx * dx + dy * dy)), 0, 1);
          }
        }
      }
    }

    canvas.addEventListener('pointerdown', function (ev) {
      ev.preventDefault();
      var cell = toCell(ev);
      stopAnim();
      if (nearMarker(ev, start)) {
        drag = { mode: 'start' };
      } else if (nearMarker(ev, goal)) {
        drag = { mode: 'goal' };
      } else {
        var erase = ev.shiftKey;
        drag = {
          mode: 'paint',
          erase: erase,
          wallValue: field.maze ? (erase ? 0 : 1) : 0,
          last: cell
        };
        paintCell(cell, erase);
        bake();
        instant();
        draw();
      }
      try { canvas.setPointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    });

    canvas.addEventListener('pointermove', function (ev) {
      if (!drag) return;
      var cell = toCell(ev);
      if (drag.mode === 'start' || drag.mode === 'goal') {
        if (field.wall[field.idx(cell[0], cell[1])]) return;
        if (drag.mode === 'start') start = cell; else goal = cell;
        instant();
        draw();
      } else {
        // Interpolate so fast strokes leave no gaps.
        var lx = drag.last[0], ly = drag.last[1];
        var steps = Math.max(Math.abs(cell[0] - lx), Math.abs(cell[1] - ly), 1);
        for (var s = 1; s <= steps; s++) {
          paintCell([
            Math.round(lx + (cell[0] - lx) * s / steps),
            Math.round(ly + (cell[1] - ly) * s / steps)
          ], drag.erase);
        }
        drag.last = cell;
        bake();
        instant();
        draw();
      }
    });

    function endDrag(ev) {
      drag = null;
      try { canvas.releasePointerCapture(ev.pointerId); } catch (e) { /* not fatal */ }
    }
    canvas.addEventListener('pointerup', endDrag);
    canvas.addEventListener('pointercancel', endDrag);

    // ---- toolbar ----------------------------------------------------------

    root.addEventListener('click', function (ev) {
      var btn = ev.target.closest ? ev.target.closest('[data-act]') : null;
      if (!btn) return;
      var act = btn.dataset.act;
      if (act === 'run') {
        if (animating) {
          stopAnim();
          runBtn.textContent = 'Run';
          updateHUD();
        } else {
          startAnim(search === null || search.status !== 'running');
        }
      } else if (act === 'step') {
        stopAnim();
        runBtn.textContent = 'Run';
        if (!search || search.status !== 'running') search = makeSearch();
        search.expand();
        rebuildOverlay();
        updateHUD();
        draw();
      } else if (act === 'reset') {
        stopAnim();
        runBtn.textContent = 'Run';
        instant();
        draw();
      } else if (act === 'reseed') {
        conf.seed = (Math.random() * 0xffffffff) >>> 0;
        regen();
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
      if (raw.indexOf('pf=1') === -1) return null;
      var p;
      try { p = new URLSearchParams(raw); } catch (e) { return null; }
      if (p.has('t') && TERRAINS[p.get('t')]) conf.terrain = p.get('t');
      if (p.has('a')) {
        for (var i = 0; i < ALGOS.length; i++) {
          if (ALGOS[i][0] === p.get('a')) conf.algo = p.get('a');
        }
      }
      var num = function (key, target, lo, hi) {
        if (!p.has(key)) return;
        var v = parseFloat(p.get(key));
        if (isFinite(v)) conf[target] = clamp(v, lo, hi);
      };
      num('T', 'T', 0.05, 1.5);
      num('w', 'w', 0, 1.5);
      num('sp', 'speed', 1, 400);
      num('res', 'res', 41, 201);
      if (p.has('s')) {
        var seed = parseInt(p.get('s'), 10);
        if (isFinite(seed)) conf.seed = seed >>> 0;
      }
      return p;
    }

    function applyMarkersFromHash(p) {
      if (!p) return;
      var pair = function (key, into) {
        if (!p.has(key)) return;
        var m = p.get(key).split(',');
        var x = parseInt(m[0], 10), y = parseInt(m[1], 10);
        if (isFinite(x) && isFinite(y)) {
          x = clamp(x, 0, field.gw - 1);
          y = clamp(y, 0, field.gh - 1);
          if (!field.wall[field.idx(x, y)]) { into[0] = x; into[1] = y; }
        }
      };
      pair('st', start);
      pair('go', goal);
    }

    function buildURL() {
      var p = new URLSearchParams();
      p.set('pf', '1');
      p.set('t', conf.terrain);
      p.set('a', conf.algo);
      p.set('T', conf.T.toFixed(2));
      p.set('w', conf.w.toFixed(2));
      p.set('sp', String(Math.round(conf.speed)));
      p.set('res', String(Math.round(conf.res)));
      p.set('s', String(conf.seed));
      p.set('st', start[0] + ',' + start[1]);
      p.set('go', goal[0] + ',' + goal[1]);
      return location.origin + location.pathname + location.search + '#' + p.toString();
    }

    // ---- boot -------------------------------------------------------------

    var hashParams = readHash();
    regen();
    applyMarkersFromHash(hashParams);
    if (hashParams) { instant(); draw(); }

    if (typeof ResizeObserver !== 'undefined') {
      new ResizeObserver(function () { draw(); }).observe(canvas);
    } else {
      window.addEventListener('resize', function () { draw(); });
    }

    var reduced = window.matchMedia &&
      window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    if (!reduced) startAnim(true);       // opening movie: watch the frontier

    root._pathfinding = {
      getField: function () { return field; },
      getSearch: function () { return search; }
    };
  }

  // -------------------------------------------------------------------------

  if (typeof window !== 'undefined') {
    window.PathfindingSim = {
      Field: Field,
      Search: Search,
      TERRAINS: TERRAINS,
      mulberry32: mulberry32
    };
  }
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
      Field: Field, Search: Search,
      TERRAINS: TERRAINS, mulberry32: mulberry32
    };
  }
  if (typeof document !== 'undefined') {
    var boot = function () {
      var nodes = document.querySelectorAll('[data-sim="pathfinding"]');
      for (var i = 0; i < nodes.length; i++) mount(nodes[i]);
    };
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', boot);
    } else {
      boot();
    }
  }
})();
