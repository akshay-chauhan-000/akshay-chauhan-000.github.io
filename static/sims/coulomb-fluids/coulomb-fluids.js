/* ==========================================================================
   Coulomb Fluids
   --------------------------------------------------------------------------
   Charged discs (diameter sigma = 1, mass 1) confined to a periodic plane,
   interacting through the full three-dimensional Coulomb law plus a WCA
   core, in an implicit solvent at temperature kT = 1:

       U_ij = lB q_i q_j / r  +  4 [ r^-12 - r^-6 ] + 1   (r < 2^(1/6))

   lB is the Bjerrum length in units of sigma; q are valences. Dynamics are
   Langevin, integrated with the BAOAB splitting (one force call per step).

   Electrostatics: smooth particle-mesh Ewald for point charges in a plane
   that is periodic in x and y but not in z. With the split
   1/r = erfc(a r)/r + erf(a r)/r,

       real space   lB q_i q_j erfc(a r)/r            cell list, r < rc
       k space      (lB / 2A) sum_{k != 0} G(k) |Q(k)|^2
                    G(k) = (2 pi / k) erfc(k / 2a)    [2-D FT of erf(ar)/r]
       self         -lB (a / sqrt(pi)) sum_i q_i^2
       background   -sqrt(pi) lB Q_net^2 / (a A)      (zero if neutral)

   Charges are spread onto the mesh with cubic cardinal B-splines and
   transformed with a radix-2 FFT. The field comes from ik-differentiation,
   with E_x and E_y packed into ONE complex inverse transform (both are real,
   so IFFT[E_x(k) + i E_y(k)] = E_x + i E_y). It is interpolated back with
   the same splines. The SPME Euler-spline factor
   |b(k)|^2 = 36 / (4 + 2 cos(2 pi m / M))^2 per axis corrects the
   interpolation. The splitting parameter adapts to the mesh spacing h,
   a = pi / (7 h), with a rc = 3.5, so k-space and real-space truncation
   errors stay balanced; the mean mesh self-energy (aliasing) error is
   computed analytically at setup and subtracted.

   No dependencies, no build step. The engine (System) is DOM-free and
   exported on window.CoulombFluids for headless use; the UI mounts into
   [data-sim="coulomb-fluids"].
   ========================================================================== */

(function () {
  'use strict';

  var INK = '#0b0d12';
  var TWO_PI = 2 * Math.PI;
  var SQRT_PI = Math.sqrt(Math.PI);
  var WCA_R = Math.pow(2, 1 / 6);
  var WCA_R2 = WCA_R * WCA_R;
  var FENE_K = 30;
  var FENE_R0SQ = 1.5 * 1.5;
  var MADELUNG_2D = 1.6155426267;              // square NaCl-type net, 1/r
  var LB_WATER_NM = 0.7135;                    // water, 298 K
  var EPS_WATER = 78.4;
  var SIGMA_NM = 0.5;                          // display mapping only
  var KT_MV = 25.69;                           // kT/e at 298 K
  var RDF_MAX = 4;
  var RDF_NB = 80;
  var THETA_MAX = 4;
  var THETA_NB = 40;
  var ZEFF_SPAN = 12;
  var ZEFF_NB = 60;
  var PROF_NB = 64;
  var SK_KMAX = 6;
  var TRACE_LEN = 360;

  var CAT_COL = ['#ff5d73', '#ffa64d', '#ffe066'];   // +1 +2 +3
  var ANI_COL = ['#6f9dff', '#3fc8d8', '#5ddb8b'];   // -1 -2 -3
  var MONO_COL = '#b57bff';
  var POS_RGB = [255, 93, 115];
  var NEG_RGB = [111, 157, 255];
  var INK_RGB = [11, 13, 18];

  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  function gcd(a, b) { while (b) { var t = a % b; a = b; b = t; } return a; }

  function pow2ceil(v) { var p = 1; while (p < v) p <<= 1; return p; }

  function mulberry32(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* erfc(x), x >= 0, given ex2 = exp(-x^2). Abramowitz & Stegun 7.1.26,
     |error| < 1.5e-7. Reusing the exponential saves one exp per pair. */
  function erfcExp(x, ex2) {
    var t = 1 / (1 + 0.3275911 * x);
    return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
      t * (-1.453152027 + t * 1.061405429)))) * ex2;
  }

  function erfc(x) {
    return x >= 0 ? erfcExp(x, Math.exp(-x * x)) : 2 - erfcExp(-x, Math.exp(-x * x));
  }

  // -------------------------------------------------------------------------
  // Radix-2 complex FFT, in place, strided (so rows and columns of a 2-D
  // row-major mesh can be transformed without copying).
  // -------------------------------------------------------------------------

  function FFT(n) {
    this.n = n;
    var levels = Math.round(Math.log(n) / Math.LN2);
    this.rev = new Uint32Array(n);
    for (var i = 0; i < n; i++) {
      var r = 0, v = i;
      for (var b = 0; b < levels; b++) { r = (r << 1) | (v & 1); v >>= 1; }
      this.rev[i] = r;
    }
    this.cos = new Float64Array(n >> 1);
    this.sin = new Float64Array(n >> 1);
    for (i = 0; i < (n >> 1); i++) {
      this.cos[i] = Math.cos(TWO_PI * i / n);
      this.sin[i] = Math.sin(TWO_PI * i / n);
    }
  }

  /* sign = -1: forward (e^{-ikx}); sign = +1: inverse, unnormalised. */
  FFT.prototype.transform = function (re, im, off, stride, sign) {
    var n = this.n, rev = this.rev, cs = this.cos, sn = this.sin;
    var i, j, k, a, b, t;
    for (i = 0; i < n; i++) {
      j = rev[i];
      if (j > i) {
        a = off + i * stride; b = off + j * stride;
        t = re[a]; re[a] = re[b]; re[b] = t;
        t = im[a]; im[a] = im[b]; im[b] = t;
      }
    }
    for (var size = 2; size <= n; size <<= 1) {
      var half = size >> 1, step = n / size;
      for (i = 0; i < n; i += size) {
        for (j = i, k = 0; j < i + half; j++, k += step) {
          a = off + j * stride; b = a + half * stride;
          var wr = cs[k], wi = sign * sn[k];
          var xr = re[b] * wr - im[b] * wi;
          var xi = re[b] * wi + im[b] * wr;
          re[b] = re[a] - xr; im[b] = im[a] - xi;
          re[a] += xr; im[a] += xi;
        }
      }
    }
  };

  // -------------------------------------------------------------------------
  // Particle-mesh Ewald for a 2-D periodic plane of 3-D Coulomb charges
  // -------------------------------------------------------------------------

  function PME() {
    this.Mx = 0;
    this.My = 0;
    this.cap = 0;
  }

  /* Splitting: a = pi / (KFAC h) puts the mesh Nyquist at k/2a = KFAC/2 = 3.5
     (erfc ~ 7e-7), and a rc = ARC = 3.5 (erfc ~ 7e-7). Measured against a
     direct Ewald sum: rms force error ~ 6e-4, energy error ~ 1e-5 kT per q^2.
     rc is capped at L/3 so the cell list stays exact; the mesh is refined to
     compensate. */
  var KFAC = 7.0, ARC = 3.5;

  PME.prototype.setup = function (Lx, Ly, hmax) {
    var rcMax = Math.min(Lx, Ly) / 3;
    var Mx = clamp(pow2ceil(Math.ceil(Lx / hmax)), 16, 256);
    var My = clamp(pow2ceil(Math.ceil(Ly / hmax)), 16, 256);
    var alpha = Math.PI / (KFAC * Math.max(Lx / Mx, Ly / My));
    var rc = Math.max(ARC / alpha, RDF_MAX);
    if (rc > rcMax) {
      rc = Math.max(rcMax, RDF_MAX);
      alpha = ARC / rc;
      var hNeed = Math.PI / (KFAC * alpha);
      Mx = clamp(pow2ceil(Math.ceil(Lx / hNeed)), 16, 256);
      My = clamp(pow2ceil(Math.ceil(Ly / hNeed)), 16, 256);
    }
    this.Lx = Lx;
    this.Ly = Ly;
    this.A = Lx * Ly;
    this.hx = Lx / Mx;
    this.hy = Ly / My;
    this.alpha = alpha;
    this.rc = rc;
    if (Mx !== this.Mx || My !== this.My) {
      this.Mx = Mx;
      this.My = My;
      var S = Mx * My;
      this.re = new Float64Array(S);
      this.im = new Float64Array(S);
      this.qr = new Float64Array(S);
      this.qi = new Float64Array(S);
      this.G = new Float64Array(S);
      this.ex = new Float64Array(S);
      this.ey = new Float64Array(S);
      this.phi = new Float64Array(S);
      this.kxd = new Float64Array(Mx);
      this.kyd = new Float64Array(My);
      this.fftX = new FFT(Mx);
      this.fftY = My === Mx ? this.fftX : new FFT(My);
    }
    var a = this.alpha, ix, iy;
    for (ix = 0; ix < Mx; ix++) {
      var nx = ix <= Mx / 2 ? ix : ix - Mx;
      this.kxd[ix] = 2 * ix === Mx ? 0 : TWO_PI * nx / Lx;   // Nyquist -> 0
    }
    for (iy = 0; iy < My; iy++) {
      var ny = iy <= My / 2 ? iy : iy - My;
      var ky = TWO_PI * ny / Ly;
      var cy = 4 + 2 * Math.cos(TWO_PI * iy / My);
      var by = 36 / (cy * cy);
      this.kyd[iy] = 2 * iy === My ? 0 : ky;
      for (ix = 0; ix < Mx; ix++) {
        var nx2 = ix <= Mx / 2 ? ix : ix - Mx;
        var kx = TWO_PI * nx2 / Lx;
        var k = Math.sqrt(kx * kx + ky * ky);
        var cx = 4 + 2 * Math.cos(TWO_PI * ix / Mx);
        this.G[iy * Mx + ix] = k === 0 ? 0 :
          TWO_PI * erfc(k / (2 * a)) / k * by * 36 / (cx * cx);
      }
    }

    // Mean mesh self-energy error per unit q^2 (per unit lB). Averaged over
    // where a particle sits in its cell, |Q(k)|^2 of a lone unit charge is
    // the alias sum  sum_j |W(k + 2 pi j / h)|^2,  W = sinc^4 per axis, so
    //   err = (1/2A) [ sum_mesh G(k) Sx Sy  -  sum_all (2 pi / k) erfc(k/2a) ]
    // and the exact Ewald energy is E_mesh - lB err sum q^2.
    var aliasX = new Float64Array(Mx), aliasY = new Float64Array(My);
    var aliasSum = function (theta) {
      var sum = 0;
      for (var j = -3; j <= 3; j++) {
        var z = 0.5 * theta + Math.PI * j;
        var sc = Math.abs(z) < 1e-12 ? 1 : Math.sin(z) / z;
        var s2 = sc * sc, s4 = s2 * s2;
        sum += s4 * s4;
      }
      return sum;
    };
    for (ix = 0; ix < Mx; ix++) aliasX[ix] = aliasSum(TWO_PI * (ix <= Mx / 2 ? ix : ix - Mx) / Mx);
    for (iy = 0; iy < My; iy++) aliasY[iy] = aliasSum(TWO_PI * (iy <= My / 2 ? iy : iy - My) / My);
    var meshSelf = 0, exactSelf = 0;
    for (iy = 0; iy < My; iy++) {
      var nyy = iy <= My / 2 ? iy : iy - My;
      for (ix = 0; ix < Mx; ix++) {
        var nxx = ix <= Mx / 2 ? ix : ix - Mx;
        meshSelf += this.G[iy * Mx + ix] * aliasX[ix] * aliasY[iy];
        for (var jy = -2; jy <= 2; jy++) {
          var kyy = TWO_PI * (nyy + jy * My) / Ly;
          for (var jx = -2; jx <= 2; jx++) {
            var kxx = TWO_PI * (nxx + jx * Mx) / Lx;
            var kk = Math.sqrt(kxx * kxx + kyy * kyy);
            if (kk > 0) exactSelf += TWO_PI * erfc(kk / (2 * a)) / kk;
          }
        }
      }
    }
    this.selfErr = (meshSelf - exactSelf) / (2 * this.A);
    this.phiValid = false;
  };

  PME.prototype.fft2 = function (re, im, sign) {
    var Mx = this.Mx, My = this.My, i;
    for (i = 0; i < My; i++) this.fftX.transform(re, im, i * Mx, 1, sign);
    for (i = 0; i < Mx; i++) this.fftY.transform(re, im, i, Mx, sign);
  };

  /* Cubic cardinal B-spline weights for fractional offset t in [0, 1),
     assigned to mesh points floor(u) - 1 ... floor(u) + 2. */
  function splineW(t, w, o) {
    var t2 = t * t, t3 = t2 * t, u = 1 - t;
    w[o] = u * u * u / 6;
    w[o + 1] = (3 * t3 - 6 * t2 + 4) / 6;
    w[o + 2] = (-3 * t3 + 3 * t2 + 3 * t + 1) / 6;
    w[o + 3] = t3 / 6;
  }

  /* Adds the long-range forces to (fx, fy); returns the k-space energy. */
  PME.prototype.compute = function (n, x, y, q, fx, fy, lB) {
    var Mx = this.Mx, My = this.My, re = this.re, im = this.im;
    if (this.cap < n) {
      this.cap = Math.max(64, 2 * n);
      this.w = new Float64Array(8 * this.cap);
      this.base = new Int32Array(2 * this.cap);
    }
    var w = this.w, base = this.base;
    var ihx = 1 / this.hx, ihy = 1 / this.hy;
    var i, a, b, o, bx, by, yy, xx, row, wy, id;
    re.fill(0);
    im.fill(0);

    for (i = 0; i < n; i++) {
      var ux = x[i] * ihx, uy = y[i] * ihy;
      var fx0 = Math.floor(ux), fy0 = Math.floor(uy);
      o = 8 * i;
      splineW(ux - fx0, w, o);
      splineW(uy - fy0, w, o + 4);
      bx = ((fx0 - 1) % Mx + Mx) % Mx;
      by = ((fy0 - 1) % My + My) % My;
      base[2 * i] = bx;
      base[2 * i + 1] = by;
      var qq = q[i];
      if (qq === 0) continue;
      for (a = 0; a < 4; a++) {
        yy = by + a; if (yy >= My) yy -= My;
        row = yy * Mx;
        wy = qq * w[o + 4 + a];
        for (b = 0; b < 4; b++) {
          xx = bx + b; if (xx >= Mx) xx -= Mx;
          re[row + xx] += wy * w[o + b];
        }
      }
    }

    this.fft2(re, im, -1);

    var G = this.G, kxd = this.kxd, kyd = this.kyd, qr = this.qr, qi = this.qi;
    var E = 0;
    for (var iy = 0; iy < My; iy++) {
      var ky = kyd[iy];
      row = iy * Mx;
      for (var ix = 0; ix < Mx; ix++) {
        id = row + ix;
        var g = G[id], ar = re[id], ai = im[id];
        qr[id] = ar;
        qi[id] = ai;
        E += g * (ar * ar + ai * ai);
        var gr = g * ar, gi = g * ai, kx = kxd[ix];
        // (-i kx + ky) (gr + i gi): real part -> E_x, imaginary -> E_y
        re[id] = ky * gr + kx * gi;
        im[id] = ky * gi - kx * gr;
      }
    }

    this.fft2(re, im, 1);

    var sc = lB / this.A, ex = this.ex, ey = this.ey, S = Mx * My;
    for (id = 0; id < S; id++) {
      ex[id] = re[id] * sc;
      ey[id] = im[id] * sc;
    }

    for (i = 0; i < n; i++) {
      if (q[i] === 0) continue;
      o = 8 * i;
      bx = base[2 * i];
      by = base[2 * i + 1];
      var Ex = 0, Ey = 0;
      for (a = 0; a < 4; a++) {
        yy = by + a; if (yy >= My) yy -= My;
        row = yy * Mx;
        wy = w[o + 4 + a];
        for (b = 0; b < 4; b++) {
          xx = bx + b; if (xx >= Mx) xx -= Mx;
          var ww = wy * w[o + b];
          Ex += ww * ex[row + xx];
          Ey += ww * ey[row + xx];
        }
      }
      fx[i] += q[i] * Ex;
      fy[i] += q[i] * Ey;
    }

    this.phiValid = false;
    this.lB = lB;
    return E * lB / (2 * this.A);
  };

  /* Smooth (Ewald-scale) potential on the mesh, from the last compute(). */
  PME.prototype.potential = function () {
    if (this.phiValid) return this.phi;
    var S = this.Mx * this.My, re = this.re, im = this.im, G = this.G, id;
    for (id = 0; id < S; id++) {
      re[id] = G[id] * this.qr[id];
      im[id] = G[id] * this.qi[id];
    }
    this.fft2(re, im, 1);
    var sc = (this.lB || 0) / this.A;
    for (id = 0; id < S; id++) this.phi[id] = re[id] * sc;
    this.phiValid = true;
    return this.phi;
  };

  /* Bilinear sample of a mesh field at (px, py). */
  PME.prototype.sample = function (field, px, py) {
    var Mx = this.Mx, My = this.My;
    var u = px / this.hx, v = py / this.hy;
    var i0 = Math.floor(u), j0 = Math.floor(v);
    var tu = u - i0, tv = v - j0;
    i0 = ((i0 % Mx) + Mx) % Mx;
    j0 = ((j0 % My) + My) % My;
    var i1 = (i0 + 1) % Mx, j1 = (j0 + 1) % My;
    return (1 - tv) * ((1 - tu) * field[j0 * Mx + i0] + tu * field[j0 * Mx + i1]) +
      tv * ((1 - tu) * field[j1 * Mx + i0] + tu * field[j1 * Mx + i1]);
  };

  // -------------------------------------------------------------------------
  // Particle system
  // -------------------------------------------------------------------------

  /* kind: 0 = small ion, 1 = chain monomer, 2 = macroion (pinned).
     Storage order is always [monomers][macroion][ions], so ions can be
     appended and swap-removed without disturbing bond indices. */
  function System(seed) {
    this.pme = new PME();
    this.lB = 1.43;
    this.gamma = 1;
    this.dt = 0.004;
    this.E0 = 0;
    this.wall = false;
    this.geom = 'bulk';
    this.n = 0;
    this.cap = 0;
    this.nMono = 0;
    this.macro = -1;
    this.macroR = 4;
    this.latticeA = 0;
    this.bonds = new Int32Array(0);
    this.nBonds = 0;
    this.wand = { on: false, x: 0, y: 0, Q: 20, a: 1.0 };
    this.seed = (seed || 1) >>> 0;
    this.rnd = mulberry32(this.seed);
    this.t = 0;
    this.collect = false;
    this.capHits = 0;
    this.Jacc = 0;
    this.Jshort = 0;
    this.Jlong = 0;
    this.ncx = 0;
    this.ensure(256);
    this.setBox(40);
    this.resetStats();
  }

  System.prototype.ensure = function (cap) {
    if (cap <= this.cap) return;
    cap = Math.max(cap, 2 * this.cap);
    var self = this;
    ['x', 'y', 'vx', 'vy', 'fx', 'fy', 'q'].forEach(function (k) {
      var nu = new Float64Array(cap);
      if (self[k]) nu.set(self[k].subarray(0, self.n));
      self[k] = nu;
    });
    ['kind', 'fixed', 'paired'].forEach(function (k) {
      var nu = new Uint8Array(cap);
      if (self[k]) nu.set(self[k].subarray(0, self.n));
      self[k] = nu;
    });
    this.order = new Int32Array(cap);
    this.cellOf = new Int32Array(cap);
    this.pairs = new Int32Array(2 * cap);
    this.psiRe = new Float64Array(cap);
    this.psiIm = new Float64Array(cap);
    this.psiN = new Uint16Array(cap);
    this.ux = new Float64Array(cap);
    this.uy = new Float64Array(cap);
    this.cap = cap;
  };

  System.prototype.setBox = function (L) {
    this.L = L;
    this.pme.setup(L, L, 1.0);
    this.rc = this.pme.rc;
    this.alpha = this.pme.alpha;
    var nc = Math.max(3, Math.floor(L / this.rc));   // L >= 3 rc by the cap (L >= 12)
    if (nc !== this.ncx) {
      this.ncx = nc;
      this.cellCount = new Int32Array(nc * nc);
      this.cellStart = new Int32Array(nc * nc + 1);
    }
    this.cellSize = L / nc;
  };

  System.prototype.add = function (x, y, q, kind) {
    this.ensure(this.n + 1);
    var i = this.n++;
    this.x[i] = x; this.y[i] = y;
    this.vx[i] = 0; this.vy[i] = 0;
    this.fx[i] = 0; this.fy[i] = 0;
    this.q[i] = q;
    this.kind[i] = kind;
    this.fixed[i] = kind === 2 ? 1 : 0;
    return i;
  };

  System.prototype.wrap = function (v) {
    var L = this.L;
    if (v >= L) v -= L; else if (v < 0) v += L;
    if (v >= L || v < 0) v -= L * Math.floor(v / L);
    return v;
  };

  /* Scenario construction. p: { zp, zm, salt, phi, Z, nm }. */
  System.prototype.build = function (geom, p) {
    var zp = clamp(Math.round(p.zp) || 1, 1, 3);
    var zm = clamp(Math.round(p.zm) || 1, 1, 3);
    var salt = Math.max(0, Math.round(p.salt || 0));
    var phi = clamp(p.phi || 0.1, 0.01, 0.5);
    var rnd = this.rnd = mulberry32(this.seed);
    var i, k;
    this.n = 0;
    this.nMono = 0;
    this.macro = -1;
    this.nBonds = 0;
    this.t = 0;
    this.geom = geom;
    this.wall = geom === 'electrode';
    this.latticeA = 0;
    this.Jacc = this.Jshort = this.Jlong = 0;

    if (geom === 'crystal') {
      var side = Math.max(12, 2 * Math.round(Math.sqrt(2 * Math.max(salt, 72)) / 2));
      var a = 1.0;                                    // contact spacing
      this.ensure(side * side + 8);
      this.setBox(side * a);
      for (var gy = 0; gy < side; gy++) {
        for (var gx = 0; gx < side; gx++) {
          this.add((gx + 0.5) * a, (gy + 0.5) * a, ((gx + gy) & 1) ? -1 : 1, 0);
        }
      }
      this.latticeA = a;
      this.thermalize(rnd);
      this.resetStats();
      return;
    }

    var charges = [];
    var nm = 0, Z = 0, R = 0;
    if (geom === 'bulk' || geom === 'electrode') {
      var g = gcd(zp, zm), units = Math.max(1, salt);
      for (k = 0; k < units * zm / g; k++) charges.push(zp);
      for (k = 0; k < units * zp / g; k++) charges.push(-zm);
    } else {
      if (geom === 'polymer') nm = clamp(6 * Math.round((p.nm || 60) / 6), 12, 240);
      if (geom === 'colloid') {
        Z = clamp(6 * Math.round((p.Z || 48) / 6), 6, 240);
        R = this.macroR = 4;
      }
      var counter = Math.round((nm || Z) / zp);
      for (k = 0; k < counter; k++) charges.push(zp);
      for (k = 0; k < salt; k++) charges.push(zp);
      for (k = 0; k < salt * zp; k++) charges.push(-1);
    }

    var area = (charges.length + nm) * Math.PI / 4;
    if (geom === 'colloid') area += Math.PI * (R + 0.5) * (R + 0.5);
    var L = Math.max(20, Math.sqrt(area / phi));
    var ringR = 0;
    if (geom === 'polymer') {
      ringR = nm * 0.97 / (1.8 * Math.PI);
      L = Math.max(L, 2 * ringR + 12);
    }
    if (geom === 'colloid') L = Math.max(L, 2 * R + 18);

    this.ensure(charges.length + nm + 8);
    this.setBox(L);

    if (geom === 'polymer') {
      var cx = L / 2, cy = L / 2, dth = 0.97 / ringR;
      for (k = 0; k < nm; k++) {
        this.add(cx + ringR * Math.cos(k * dth), cy + ringR * Math.sin(k * dth), -1, 1);
      }
      this.nMono = nm;
      this.bonds = new Int32Array(2 * (nm - 1));
      for (k = 0; k < nm - 1; k++) { this.bonds[2 * k] = k; this.bonds[2 * k + 1] = k + 1; }
      this.nBonds = nm - 1;
    }
    if (geom === 'colloid') {
      this.macro = this.add(L / 2, L / 2, -Z, 2);
    }
    this.placeIons(charges, rnd);
    this.thermalize(rnd);
    this.resetStats();
  };

  /* Jittered-lattice placement avoiding walls, the chain and the colloid. */
  System.prototype.placeIons = function (charges, rnd) {
    var need = charges.length;
    if (!need) return;
    var L = this.L, self = this, i;
    var ok = function (px, py) {
      if (self.wall && (px < 1.25 || px > L - 1.25)) return false;
      if (self.macro >= 0) {
        var dx = px - self.x[self.macro], dy = py - self.y[self.macro];
        if (dx * dx + dy * dy < (self.macroR + 1.1) * (self.macroR + 1.1)) return false;
      }
      for (var m = 0; m < self.nMono; m++) {
        var ex = px - self.x[m], ey = py - self.y[m];
        if (ex * ex + ey * ey < 1.1) return false;
      }
      return true;
    };
    var s = L / Math.ceil(Math.sqrt(need * 1.2 + this.nMono + 2 * this.macroR * this.macroR * (this.macro >= 0 ? 4 : 0)));
    var sites = [];
    for (var attempt = 0; attempt < 30; attempt++) {
      sites.length = 0;
      var m = Math.floor(L / s);
      for (var gy = 0; gy < m; gy++) {
        for (var gx = 0; gx < m; gx++) {
          var px = (gx + 0.5) * L / m, py = (gy + 0.5) * L / m;
          if (ok(px, py)) sites.push(px, py);
        }
      }
      if (sites.length / 2 >= need) break;
      s *= 0.95;
    }
    var nSites = sites.length / 2;
    for (i = nSites - 1; i > 0; i--) {                // Fisher-Yates
      var j = Math.floor(rnd() * (i + 1));
      var tx = sites[2 * i], ty = sites[2 * i + 1];
      sites[2 * i] = sites[2 * j]; sites[2 * i + 1] = sites[2 * j + 1];
      sites[2 * j] = tx; sites[2 * j + 1] = ty;
    }
    var spacing = L / Math.max(1, Math.floor(L / s));
    var jit = Math.max(0, (spacing - 1.05) * 0.5);
    for (i = 0; i < need; i++) {
      var sx, sy;
      if (i < nSites) {
        sx = sites[2 * i] + (2 * rnd() - 1) * jit;
        sy = sites[2 * i + 1] + (2 * rnd() - 1) * jit;
        if (!ok(sx, sy)) { sx = sites[2 * i]; sy = sites[2 * i + 1]; }   // jitter landed badly
      } else {
        sx = rnd() * L; sy = rnd() * L;               // overfull: let WCA sort it out
      }
      if (this.wall) sx = clamp(sx, 1.2, L - 1.2);
      this.add(this.wrap(sx), this.wrap(sy), charges[i], 0);
    }
  };

  System.prototype.thermalize = function (rnd) {
    for (var i = 0; i < this.n; i++) {
      if (this.fixed[i]) { this.vx[i] = this.vy[i] = 0; continue; }
      var u1 = rnd() || 1e-12, u2 = rnd();
      var r = Math.sqrt(-2 * Math.log(u1));
      this.vx[i] = r * Math.cos(TWO_PI * u2);
      this.vy[i] = r * Math.sin(TWO_PI * u2);
    }
    this.forces();
  };

  System.prototype.rescale = function (Lnew) {
    var f = Lnew / this.L;
    for (var i = 0; i < this.n; i++) { this.x[i] *= f; this.y[i] *= f; }
    this.setBox(Lnew);
    for (i = 0; i < this.n; i++) { this.x[i] = this.wrap(this.x[i]); this.y[i] = this.wrap(this.y[i]); }
    this.forces();
  };

  System.prototype.addIon = function (px, py, q) {
    var i = this.add(this.wrap(px), this.wrap(py), q, 0);
    var u1 = this.rnd() || 1e-12, u2 = this.rnd(), r = Math.sqrt(-2 * Math.log(u1));
    this.vx[i] = r * Math.cos(TWO_PI * u2);
    this.vy[i] = r * Math.sin(TWO_PI * u2);
    this.forces();
    return i;
  };

  System.prototype.removeNear = function (px, py, rmax) {
    var best = -1, bd = rmax * rmax, L = this.L, h = L / 2;
    for (var i = 0; i < this.n; i++) {
      if (this.kind[i] !== 0) continue;
      var dx = this.x[i] - px, dy = this.y[i] - py;
      if (dx > h) dx -= L; else if (dx < -h) dx += L;
      if (dy > h) dy -= L; else if (dy < -h) dy += L;
      var d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0) return false;
    var last = this.n - 1;
    if (best !== last) {
      var keys = ['x', 'y', 'vx', 'vy', 'fx', 'fy', 'q', 'kind', 'fixed'];
      for (var k = 0; k < keys.length; k++) this[keys[k]][best] = this[keys[k]][last];
    }
    this.n--;
    this.forces();
    return true;
  };

  System.prototype.buildCells = function () {
    var n = this.n, nc = this.ncx, inv = 1 / this.cellSize;
    var count = this.cellCount, start = this.cellStart, order = this.order, cellOf = this.cellOf;
    var x = this.x, y = this.y, i, c;
    count.fill(0);
    for (i = 0; i < n; i++) {
      var cx = (x[i] * inv) | 0, cy = (y[i] * inv) | 0;
      if (cx >= nc) cx = nc - 1; else if (cx < 0) cx = 0;
      if (cy >= nc) cy = nc - 1; else if (cy < 0) cy = 0;
      c = cy * nc + cx;
      cellOf[i] = c;
      count[c]++;
    }
    start[0] = 0;
    for (c = 0; c < nc * nc; c++) start[c + 1] = start[c] + count[c];
    for (c = 0; c < nc * nc; c++) count[c] = start[c];
    for (i = 0; i < n; i++) order[count[cellOf[i]]++] = i;
  };

  var STENCIL = [0, 0, 1, 0, -1, 1, 0, 1, 1, 1];   // half shell (dx, dy)

  System.prototype.forces = function () {
    var n = this.n, x = this.x, y = this.y, q = this.q, kind = this.kind;
    var fx = this.fx, fy = this.fy, L = this.L, h = L / 2;
    var lB = this.lB, alpha = this.alpha, rc2 = this.rc * this.rc;
    var c2a = 2 * alpha / SQRT_PI;
    var i, j, k;
    fx.fill(0, 0, n);
    fy.fill(0, 0, n);

    var sumQ = 0, sumQ2 = 0;
    for (i = 0; i < n; i++) { sumQ += q[i]; sumQ2 += q[i] * q[i]; }
    this.netQ = sumQ;
    this.eLong = lB > 0 ? this.pme.compute(n, x, y, q, fx, fy, lB) : 0;
    // analytic self term, plus removal of the mean mesh self-energy error
    this.eSelf = -lB * (alpha / SQRT_PI + this.pme.selfErr) * sumQ2;
    this.eBg = -SQRT_PI * lB * sumQ * sumQ / (alpha * L * L);

    this.buildCells();

    var collect = this.collect;
    var rdf = this.rdfFrame, rdfInv = RDF_NB / RDF_MAX, rdf2 = RDF_MAX * RDF_MAX;
    var paired = this.paired, pairs = this.pairs, nPairs = 0;
    var psi = collect && this.geom === 'crystal';
    var psiR2 = psi ? 1.25 * 1.25 : 0;                // between contact and sqrt2 x contact
    var psiRe = this.psiRe, psiIm = this.psiIm, psiN = this.psiN;
    if (collect) {
      rdf.fill(0);
      paired.fill(0, 0, n);
      if (psi) { psiRe.fill(0, 0, n); psiIm.fill(0, 0, n); psiN.fill(0, 0, n); }
    }

    var eReal = 0, eWca = 0;
    var nc = this.ncx, start = this.cellStart, order = this.order;
    for (var cy = 0; cy < nc; cy++) {
      for (var cx = 0; cx < nc; cx++) {
        var c = cy * nc + cx;
        for (var s = 0; s < 10; s += 2) {
          var ox = cx + STENCIL[s], oy = cy + STENCIL[s + 1];
          if (ox >= nc) ox -= nc; else if (ox < 0) ox += nc;
          if (oy >= nc) oy -= nc;
          var c2 = oy * nc + ox;
          for (var ia = start[c]; ia < start[c + 1]; ia++) {
            i = order[ia];
            var xi = x[i], yi = y[i], qi = q[i], ki = kind[i];
            var jb0 = s === 0 ? ia + 1 : start[c2];
            for (var jb = jb0; jb < start[c2 + 1]; jb++) {
              j = order[jb];
              var dx = x[j] - xi, dy = y[j] - yi;
              if (dx > h) dx -= L; else if (dx < -h) dx += L;
              if (dy > h) dy -= L; else if (dy < -h) dy += L;
              var r2 = dx * dx + dy * dy;
              if (r2 >= rc2) continue;
              var fr = 0, qq = qi * q[j], kj = kind[j];
              var r = Math.sqrt(r2);
              if (qq !== 0 && lB > 0) {
                var ar = alpha * r, ex2 = Math.exp(-ar * ar);
                var ec = erfcExp(ar, ex2);
                fr += lB * qq * (ec / r + c2a * ex2) / r2;
                eReal += lB * qq * ec / r;
              }
              if (r2 < WCA_R2 && ki !== 2 && kj !== 2) {
                var ir2 = 1 / r2, s6 = ir2 * ir2 * ir2;
                fr += 24 * s6 * (2 * s6 - 1) * ir2;
                eWca += 4 * s6 * (s6 - 1) + 1;
              }
              if (fr !== 0) {
                fx[i] -= fr * dx; fy[i] -= fr * dy;
                fx[j] += fr * dx; fy[j] += fr * dy;
              }
              if (collect) {
                if (ki === 0 && kj === 0) {
                  if (r2 < rdf2 && qq !== 0) {
                    var cls = qq < 0 ? 0 : (qi > 0 ? 1 : 2);
                    rdf[cls * RDF_NB + ((r * rdfInv) | 0)]++;
                  }
                  if (qq < 0) {
                    // Bjerrum radius for a planar (r dr) measure, capped at 2 sigma
                    var rB = lB * -qq;
                    if (rB > 1 && r < (rB < 2 ? rB : 2)) {
                      paired[i] = 1; paired[j] = 1;
                      if (nPairs < this.cap) { pairs[2 * nPairs] = i; pairs[2 * nPairs + 1] = j; nPairs++; }
                    }
                  }
                }
                if (psi && r2 < psiR2) {
                  var cs = dx / r, sn = dy / r;
                  var c2t = cs * cs - sn * sn, s2t = 2 * cs * sn;
                  var c4 = c2t * c2t - s2t * s2t, s4 = 2 * c2t * s2t;
                  psiRe[i] += c4; psiIm[i] += s4; psiN[i]++;
                  psiRe[j] += c4; psiIm[j] += s4; psiN[j]++;
                }
              }
            }
          }
        }
      }
    }
    this.nPairs = nPairs;

    // Colloid: WCA shifted by its radius (contact at r = R + 1/2).
    if (this.macro >= 0) {
      var M = this.macro, D = this.macroR - 0.5, rcut = D + WCA_R;
      var xm = x[M], ym = y[M];
      for (i = 0; i < n; i++) {
        if (i === M) continue;
        var mx = x[i] - xm, my = y[i] - ym;
        if (mx > h) mx -= L; else if (mx < -h) mx += L;
        if (my > h) my -= L; else if (my < -h) my += L;
        var mr2 = mx * mx + my * my;
        if (mr2 >= rcut * rcut) continue;
        var mr = Math.sqrt(mr2), rs = Math.max(mr - D, 0.6);
        var is2 = 1 / (rs * rs), m6 = is2 * is2 * is2;
        var fm = 24 * m6 * (2 * m6 - 1) / rs / mr;
        eWca += 4 * m6 * (m6 - 1) + 1;
        fx[i] += fm * mx; fy[i] += fm * my;
      }
    }

    // FENE bonds (Kremer-Grest), clamped short of the divergence.
    var eBond = 0, bonds = this.bonds;
    for (k = 0; k < this.nBonds; k++) {
      i = bonds[2 * k]; j = bonds[2 * k + 1];
      var bx = x[j] - x[i], by = y[j] - y[i];
      if (bx > h) bx -= L; else if (bx < -h) bx += L;
      if (by > h) by -= L; else if (by < -h) by += L;
      var u = Math.min((bx * bx + by * by) / FENE_R0SQ, 0.98);
      var fb = -FENE_K / (1 - u);
      eBond += -0.5 * FENE_K * FENE_R0SQ * Math.log(1 - u);
      fx[i] -= fb * bx; fy[i] -= fb * by;
      fx[j] += fb * bx; fy[j] += fb * by;
    }

    // Electrode walls at x = 0 (= L): WCA in the distance to the wall.
    var eWall = 0;
    if (this.wall) {
      for (i = 0; i < n; i++) {
        var d1 = Math.max(x[i], 0.6), d2 = Math.max(L - x[i], 0.6);
        if (d1 < WCA_R) {
          var w1 = 1 / (d1 * d1), w6 = w1 * w1 * w1;
          fx[i] += 24 * w6 * (2 * w6 - 1) / d1;
          eWall += 4 * w6 * (w6 - 1) + 1;
        }
        if (d2 < WCA_R) {
          var v1 = 1 / (d2 * d2), v6 = v1 * v1 * v1;
          fx[i] -= 24 * v6 * (2 * v6 - 1) / d2;
          eWall += 4 * v6 * (v6 - 1) + 1;
        }
      }
    }

    // Uniform applied field.
    if (this.E0 !== 0) {
      for (i = 0; i < n; i++) fx[i] += q[i] * this.E0;
    }

    // The wand: a soft-cored point charge held by the pointer, tapered to
    // zero between 0.3 L and 0.5 L so the minimum image never jumps.
    var wd = this.wand;
    if (wd.on && wd.Q !== 0 && lB > 0) {
      var a2 = wd.a * wd.a, r1 = 0.3 * L, r2t = 0.5 * L;
      for (i = 0; i < n; i++) {
        if (q[i] === 0 || this.fixed[i]) continue;
        var wx = x[i] - wd.x, wy = y[i] - wd.y;
        if (wx > h) wx -= L; else if (wx < -h) wx += L;
        if (wy > h) wy -= L; else if (wy < -h) wy += L;
        var wr2 = wx * wx + wy * wy, wr = Math.sqrt(wr2);
        if (wr >= r2t) continue;
        var sw = wr <= r1 ? 1 : 0.5 * (1 + Math.cos(Math.PI * (wr - r1) / (r2t - r1)));
        var d3 = wr2 + a2, fw = sw * lB * wd.Q * q[i] / (d3 * Math.sqrt(d3));
        fx[i] += fw * wx; fy[i] += fw * wy;
      }
    }

    this.eReal = eReal;
    this.eWca = eWca;
    this.eBond = eBond;
    this.eWall = eWall;
  };

  /* One BAOAB Langevin step. Returns false on a non-finite state. */
  System.prototype.step = function () {
    var n = this.n, dt = this.dt, hdt = 0.5 * dt;
    var x = this.x, y = this.y, vx = this.vx, vy = this.vy;
    var fx = this.fx, fy = this.fy, q = this.q, fixed = this.fixed;
    var L = this.L, wall = this.wall, rnd = this.rnd;
    var c1 = Math.exp(-this.gamma * dt), c2 = Math.sqrt(Math.max(0, 1 - c1 * c1));
    var vmax = Math.min(0.3 / dt, 100), vmax2 = vmax * vmax;
    var jq = 0, i;

    for (i = 0; i < n; i++) {
      if (fixed[i]) continue;
      vx[i] += hdt * fx[i];
      vy[i] += hdt * fy[i];
      var v2 = vx[i] * vx[i] + vy[i] * vy[i];
      if (v2 > vmax2) {
        var sc = vmax / Math.sqrt(v2);
        vx[i] *= sc; vy[i] *= sc;
        this.capHits++;
      }
      var dx = hdt * vx[i];
      x[i] += dx;
      y[i] += hdt * vy[i];
      jq += q[i] * dx;
    }
    if (c2 > 0) {
      for (i = 0; i < n; i++) {
        if (fixed[i]) continue;
        var u1 = rnd() || 1e-12, u2 = rnd();
        var rr = c2 * Math.sqrt(-2 * Math.log(u1)), th = TWO_PI * u2;
        vx[i] = c1 * vx[i] + rr * Math.cos(th);
        vy[i] = c1 * vy[i] + rr * Math.sin(th);
      }
    }
    for (i = 0; i < n; i++) {
      if (fixed[i]) continue;
      var dx2 = hdt * vx[i];
      x[i] += dx2;
      y[i] += hdt * vy[i];
      jq += q[i] * dx2;
      if (wall) {
        if (x[i] < 0.5) { x[i] = 0.5; if (vx[i] < 0) vx[i] = -vx[i]; }
        else if (x[i] > L - 0.5) { x[i] = L - 0.5; if (vx[i] > 0) vx[i] = -vx[i]; }
      } else if (x[i] >= L || x[i] < 0) {
        x[i] = this.wrap(x[i]);
      }
      if (y[i] >= L || y[i] < 0) y[i] = this.wrap(y[i]);
    }

    this.forces();

    for (i = 0; i < n; i++) {
      if (fixed[i]) continue;
      vx[i] += hdt * fx[i];
      vy[i] += hdt * fy[i];
    }

    // Charge current density J_x = sum q dx / (A dt): running means over ~5 tau
    // and ~60 tau. Weight 1/n until n reaches the window, so neither average is
    // biased towards its zero start.
    var Jinst = jq / (L * L * dt);
    this.jn++;
    this.Jshort += (Jinst - this.Jshort) * Math.max(dt / 5, 1 / this.jn);
    this.Jlong += (Jinst - this.Jlong) * Math.max(dt / 60, 1 / this.jn);
    this.t += dt;
    return isFinite(this.eLong) && isFinite(this.eReal) && isFinite(x[0]) && isFinite(vx[n - 1 > 0 ? n - 1 : 0]);
  };

  // -------------------------------------------------------------------------
  // Measurements
  // -------------------------------------------------------------------------

  System.prototype.resetStats = function () {
    this.rdfFrame = new Float64Array(3 * RDF_NB);
    this.rdfAvg = new Float64Array(3 * RDF_NB);
    this.nkMax = Math.max(1, Math.min(96, Math.floor(SK_KMAX * this.L / TWO_PI)));
    this.skAxis = new Float64Array(this.nkMax + 1);
    this.skDiag = new Float64Array(this.nkMax + 1);
    this.skTmp = new Float64Array(8 * (this.nkMax + 1));
    this.profP = new Float64Array(PROF_NB);
    this.profM = new Float64Array(PROF_NB);
    this.profPhi = new Float64Array(this.pme.Mx);
    this.theta = new Float64Array(THETA_NB);
    this.zeff = new Float64Array(ZEFF_NB);
    this.samples = 0;
    this.obs = this.obs || {};
    this.Jshort = this.Jlong = 0;
    this.jn = 0;
  };

  function ema(avg, frame, w) {
    for (var i = 0; i < avg.length; i++) avg[i] += (frame[i] - avg[i]) * w;
  }

  /* Call once per rendered frame, right after a step taken with collect on. */
  System.prototype.measure = function () {
    var n = this.n, L = this.L, A = L * L, q = this.q, kind = this.kind;
    var o = this.obs, i;
    this.samples++;
    var w = 1 / Math.min(this.samples, 240);

    // Kinetic temperature from fluctuations about each class's mean velocity
    // (cations, anions, monomers), so field-driven drift is not counted as heat.
    var cls = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];   // sum vx, vy, v^2, count
    var nIon = 0, nPos = 0, nNeg = 0, sq2Ion = 0, nCharged = 0, sumQ2 = 0, nPaired = 0;
    for (i = 0; i < n; i++) {
      if (!this.fixed[i]) {
        var cc = cls[kind[i] === 1 ? 2 : (q[i] > 0 ? 0 : 1)];
        cc[0] += this.vx[i]; cc[1] += this.vy[i];
        cc[2] += this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i]; cc[3]++;
      }
      if (q[i] !== 0) { nCharged++; sumQ2 += q[i] * q[i]; }
      if (kind[i] === 0 && q[i] !== 0) {
        nIon++; sq2Ion += q[i] * q[i];
        if (q[i] > 0) nPos++; else nNeg++;
        if (this.paired[i]) nPaired++;
      }
    }
    var keF = 0, dof = 0;
    for (var ci = 0; ci < 3; ci++) {
      var C = cls[ci];
      if (C[3] < 2) continue;
      keF += C[2] - (C[0] * C[0] + C[1] * C[1]) / C[3];
      dof += 2 * (C[3] - 1);
    }
    o.n = n; o.nIon = nIon; o.nPos = nPos; o.nNeg = nNeg; o.geom = this.geom;
    o.T = dof ? keF / dof : 0;
    o.eC = nCharged ? (this.eReal + this.eLong + this.eSelf + this.eBg) / nCharged : 0;
    o.eCtot = this.eReal + this.eLong + this.eSelf + this.eBg;
    o.pairFrac = nIon ? nPaired / nIon : 0;
    o.kappa = TWO_PI * this.lB * sq2Ion / A;
    o.q2 = nIon ? sq2Ion / nIon : 0;
    o.aws = nIon ? Math.sqrt(A / (Math.PI * nIon)) : 0;
    o.Gamma = o.aws ? this.lB * o.q2 / o.aws : 0;
    o.density = nIon / A;
    o.netQ = this.netQ;
    o.Jshort = this.Jshort;
    o.Jlong = this.Jlong;
    var mobQ2 = 0;
    for (i = 0; i < n; i++) if (!this.fixed[i]) mobQ2 += q[i] * q[i];
    o.lamNE = mobQ2 / (this.gamma * A);              // Nernst-Einstein, D = 1/gamma
    o.Jideal = o.lamNE * this.E0;
    o.lamRatio = Math.abs(this.E0) > 1e-3 ? this.Jlong / o.Jideal : NaN;
    o.capHits = this.capHits;
    this.capHits = 0;

    // g(r) by charge class among small ions.
    var nP = nPos, nN = nNeg, dr = RDF_MAX / RDF_NB, frame = this.rdfFrame;
    var norm = [nP * nN / A, nP * (nP - 1) / (2 * A), nN * (nN - 1) / (2 * A)];
    for (var cl = 0; cl < 3; cl++) {
      for (var b = 0; b < RDF_NB; b++) {
        var r1 = b * dr, r2 = r1 + dr;
        var ideal = norm[cl] * Math.PI * (r2 * r2 - r1 * r1);
        frame[cl * RDF_NB + b] = ideal > 0 ? frame[cl * RDF_NB + b] / ideal : 0;
      }
    }
    ema(this.rdfAvg, frame, w);

    if (this.geom === 'bulk' || this.geom === 'crystal') this.chargeStructure(w);
    if (this.geom === 'crystal') this.measureCrystal();
    if (this.geom === 'electrode') this.measureElectrode(w);
    if (this.geom === 'polymer') this.measurePolymer(w);
    if (this.geom === 'colloid') this.measureColloid(w);
    return o;
  };

  /* S_ZZ(k) / <q^2> along the axes and the diagonals, by direct sum. */
  System.prototype.chargeStructure = function (w) {
    var n = this.n, x = this.x, y = this.y, q = this.q, kind = this.kind;
    var nk = this.nkMax, k1 = TWO_PI / this.L, acc = this.skTmp, i, m, d;
    acc.fill(0);
    var sq2 = 0;
    for (i = 0; i < n; i++) {
      var qi = q[i];
      if (kind[i] !== 0 || qi === 0) continue;
      sq2 += qi * qi;
      var ph = [k1 * x[i], k1 * y[i], k1 * (x[i] + y[i]), k1 * (x[i] - y[i])];
      for (d = 0; d < 4; d++) {
        var c1 = Math.cos(ph[d]), s1 = Math.sin(ph[d]), cr = c1, ci = s1;
        var off = d * 2 * (nk + 1);
        for (m = 1; m <= nk; m++) {
          acc[off + 2 * m] += qi * cr;
          acc[off + 2 * m + 1] += qi * ci;
          var t = cr * c1 - ci * s1;
          ci = cr * s1 + ci * c1;
          cr = t;
        }
      }
    }
    if (!sq2) return;
    for (m = 1; m <= nk; m++) {
      var s = [0, 0, 0, 0];
      for (d = 0; d < 4; d++) {
        var o2 = d * 2 * (nk + 1) + 2 * m;
        s[d] = (acc[o2] * acc[o2] + acc[o2 + 1] * acc[o2 + 1]) / sq2;
      }
      this.skAxis[m] += (0.5 * (s[0] + s[1]) - this.skAxis[m]) * w;
      this.skDiag[m] += (0.5 * (s[2] + s[3]) - this.skDiag[m]) * w;
    }
  };

  System.prototype.measureCrystal = function () {
    var re = 0, im = 0, cnt = 0;
    for (var i = 0; i < this.n; i++) {
      if (this.psiN[i]) {
        re += this.psiRe[i] / this.psiN[i];
        im += this.psiIm[i] / this.psiN[i];
      }
      cnt++;
    }
    this.obs.psi4 = cnt ? Math.sqrt(re * re + im * im) / cnt : 0;
    this.obs.madelung = this.latticeA ? this.obs.eC / (this.lB / this.latticeA) : NaN;
  };

  System.prototype.measureElectrode = function (w) {
    var n = this.n, L = this.L, x = this.x, q = this.q, pme = this.pme, i, b;
    var fp = new Float64Array(PROF_NB), fm = new Float64Array(PROF_NB);
    var binA = (L / PROF_NB) * L, stored = 0;
    for (i = 0; i < n; i++) {
      if (q[i] === 0) continue;
      b = clamp((x[i] / L * PROF_NB) | 0, 0, PROF_NB - 1);
      if (q[i] > 0) fp[b] += 1 / binA; else fm[b] += 1 / binA;
      if (x[i] < L / 2) stored += q[i];
    }
    ema(this.profP, fp, w);
    ema(this.profM, fm, w);
    var Mx = pme.Mx, My = pme.My, phi = pme.potential(), ex = pme.ex;
    var col = new Float64Array(Mx), Eb = 0, nb = 0;
    for (var ix = 0; ix < Mx; ix++) {
      var acc = 0, eacc = 0;
      for (var iy = 0; iy < My; iy++) { acc += phi[iy * Mx + ix]; eacc += ex[iy * Mx + ix]; }
      col[ix] = acc / My - this.E0 * ix * pme.hx;
      if (ix >= Mx / 3 && ix < 2 * Mx / 3) { Eb += eacc / My; nb++; }
    }
    if (this.profPhi.length !== Mx) this.profPhi = new Float64Array(Mx);
    ema(this.profPhi, col, w);
    var o = this.obs;
    o.stored = stored / L;
    o.V = this.E0 * L;
    o.Ebulk = this.E0 + (nb ? Eb / nb : 0);
    var scr = Math.abs(this.E0) > 1e-6 ? o.Ebulk / this.E0 : NaN;
    // time-averaged like the profiles (the instantaneous value is noisy)
    o.screen = (isFinite(o.screen) && isFinite(scr) && this.samples > 1) ? o.screen + (scr - o.screen) * Math.max(w, 0.02) : scr;
    o.cap = Math.abs(o.V) > 1e-6 ? -o.stored / o.V : NaN;
  };

  System.prototype.measurePolymer = function (w) {
    var nm = this.nMono, x = this.x, y = this.y, L = this.L, h = L / 2, i;
    if (nm < 2) return;
    var ux = this.ux, uy = this.uy, bsum = 0;
    ux[0] = x[0]; uy[0] = y[0];
    for (i = 1; i < nm; i++) {
      var dx = x[i] - x[i - 1], dy = y[i] - y[i - 1];
      if (dx > h) dx -= L; else if (dx < -h) dx += L;
      if (dy > h) dy -= L; else if (dy < -h) dy += L;
      ux[i] = ux[i - 1] + dx; uy[i] = uy[i - 1] + dy;
      bsum += Math.sqrt(dx * dx + dy * dy);
    }
    var cx = 0, cy = 0;
    for (i = 0; i < nm; i++) { cx += ux[i]; cy += uy[i]; }
    cx /= nm; cy /= nm;
    var rg2 = 0;
    for (i = 0; i < nm; i++) rg2 += (ux[i] - cx) * (ux[i] - cx) + (uy[i] - cy) * (uy[i] - cy);
    var o = this.obs;
    o.rg = Math.sqrt(rg2 / nm);
    o.re = Math.sqrt(Math.pow(ux[nm - 1] - ux[0], 2) + Math.pow(uy[nm - 1] - uy[0], 2));
    o.bond = bsum / (nm - 1);
    o.chainQ = 0;
    for (i = 0; i < nm; i++) o.chainQ += this.q[i];
    o.xi = this.lB / o.bond;

    // Cumulative neutralisation theta(d): net small-ion charge whose nearest
    // monomer lies within d, over |chain charge|.
    var frame = new Float64Array(THETA_NB), dBin = THETA_MAX / THETA_NB;
    var nc = this.ncx, start = this.cellStart, order = this.order, cellOf = this.cellOf;
    for (i = 0; i < this.n; i++) {
      if (this.kind[i] !== 0 || this.q[i] === 0) continue;
      var c = cellOf[i], ccx = c % nc, ccy = (c / nc) | 0, best = THETA_MAX * THETA_MAX;
      for (var oy = -1; oy <= 1; oy++) {
        for (var ox = -1; ox <= 1; ox++) {
          var nx = (ccx + ox + nc) % nc, ny = (ccy + oy + nc) % nc, c2 = ny * nc + nx;
          for (var a = start[c2]; a < start[c2 + 1]; a++) {
            var j = order[a];
            if (this.kind[j] !== 1) continue;
            var ex = x[j] - x[i], ey = y[j] - y[i];
            if (ex > h) ex -= L; else if (ex < -h) ex += L;
            if (ey > h) ey -= L; else if (ey < -h) ey += L;
            var d2 = ex * ex + ey * ey;
            if (d2 < best) best = d2;
          }
        }
      }
      if (best < THETA_MAX * THETA_MAX) {
        var bb = clamp((Math.sqrt(best) / dBin) | 0, 0, THETA_NB - 1);
        frame[bb] += this.q[i];
      }
    }
    var Qc = Math.abs(o.chainQ) || 1, run = 0;
    for (var b = 0; b < THETA_NB; b++) { run += frame[b]; frame[b] = run / Qc; }
    ema(this.theta, frame, w);
    o.theta2 = this.theta[Math.round(2 / dBin) - 1];
    var z = 1;
    for (i = 0; i < this.n; i++) if (this.kind[i] === 0 && this.q[i] > 0) { z = this.q[i]; break; }
    o.zc = z;
    o.manning = z * o.xi > 1 ? 1 - 1 / (z * o.xi) : 0;
  };

  System.prototype.measureColloid = function (w) {
    var M = this.macro;
    if (M < 0) return;
    var x = this.x, y = this.y, q = this.q, L = this.L, h = L / 2, R = this.macroR;
    var frame = new Float64Array(ZEFF_NB), dBin = ZEFF_SPAN / ZEFF_NB;
    for (var i = 0; i < this.n; i++) {
      if (i === M || q[i] === 0) continue;
      var dx = x[i] - x[M], dy = y[i] - y[M];
      if (dx > h) dx -= L; else if (dx < -h) dx += L;
      if (dy > h) dy -= L; else if (dy < -h) dy += L;
      var r = Math.sqrt(dx * dx + dy * dy) - R;
      if (r < ZEFF_SPAN) frame[clamp((r / dBin) | 0, 0, ZEFF_NB - 1)] += q[i];
    }
    var Z = q[M], run = Z;
    for (var b = 0; b < ZEFF_NB; b++) { run += frame[b]; frame[b] = run / Math.abs(Z); }
    ema(this.zeff, frame, w);
    var o = this.obs, mx = -Infinity;
    for (b = 0; b < ZEFF_NB; b++) if (this.zeff[b] > mx) mx = this.zeff[b];
    o.Z = Z;
    o.zeffMax = mx;
    o.zeff1 = this.zeff[Math.round(1.5 / dBin) - 1];
    o.zeff2 = this.zeff[Math.round(2.5 / dBin) - 1];
  };

  // -------------------------------------------------------------------------
  // Scenarios
  // -------------------------------------------------------------------------

  var SCENARIOS = [
    { id: 'salt', geom: 'bulk', label: 'Salt solution (1:1, water)' },
    { id: 'strong', geom: 'bulk', label: 'Strong coupling (3:1 salt)' },
    { id: 'crystal', geom: 'crystal', label: 'Ionic crystal (melting)' },
    { id: 'electrode', geom: 'electrode', label: 'Electrode double layers' },
    { id: 'polymer', geom: 'polymer', label: 'Polyelectrolyte + counterions' },
    { id: 'colloid', geom: 'colloid', label: 'Charged colloid (counterion cloud)' }
  ];

  /* Presets were chosen from headless runs (see the write-up): e.g. the
     polyelectrolyte at lB = 5 sits right on Manning for z = 1 and collapses
     (Rg 10.4 -> 5.3) for z = 3; the crystal at lB = 40 keeps psi4 ~ 0.87. */
  var PRESETS = {
    salt: { lB: 1.43, zp: 1, zm: 1, salt: 300, phi: 0.12, E0: 0, V: 0, Z: 48, nm: 60, gamma: 1, dt: 0.004, spf: 4, trace: 'ec', overlay: 'phi' },
    strong: { lB: 5, zp: 3, zm: 1, salt: 150, phi: 0.12, E0: 0, V: 0, Z: 48, nm: 60, gamma: 1, dt: 0.002, spf: 6, trace: 'pair', overlay: 'pairs' },
    crystal: { lB: 40, zp: 1, zm: 1, salt: 200, phi: 0.6, E0: 0, V: 0, Z: 48, nm: 60, gamma: 1, dt: 0.0015, spf: 8, trace: 'psi4', overlay: 'phi' },
    electrode: { lB: 1.43, zp: 1, zm: 1, salt: 200, phi: 0.15, E0: 0, V: 20, Z: 48, nm: 60, gamma: 0.5, dt: 0.004, spf: 12, trace: 'J', overlay: 'phi' },
    polymer: { lB: 5, zp: 1, zm: 1, salt: 0, phi: 0.05, E0: 0, V: 0, Z: 48, nm: 60, gamma: 1, dt: 0.002, spf: 18, trace: 'rg', overlay: 'phi' },
    colloid: { lB: 3, zp: 3, zm: 1, salt: 40, phi: 0.08, E0: 0, V: 0, Z: 48, nm: 60, gamma: 1, dt: 0.003, spf: 12, trace: 'zeff', overlay: 'phi' }
  };

  function scenarioById(id) {
    for (var i = 0; i < SCENARIOS.length; i++) if (SCENARIOS[i].id === id) return SCENARIOS[i];
    return SCENARIOS[0];
  }

  // -------------------------------------------------------------------------
  // UI
  // -------------------------------------------------------------------------

  function logMap(lo, hi) {
    var r = Math.log(hi / lo);
    return {
      toVal: function (p) { return lo * Math.exp(r * p / 1000); },
      toPos: function (v) { return 1000 * Math.log(clamp(v, lo, hi) / lo) / r; }
    };
  }

  function epsFor(lB) { return EPS_WATER * (LB_WATER_NM / SIGMA_NM) / lB; }

  var SLIDERS = [
    { key: 'lB', label: 'Bjerrum length ℓB / σ', log: logMap(0.1, 60), rebuild: false,
      fmt: function (v) { return v.toFixed(v < 10 ? 2 : 1) + '  (εr ≈ ' + Math.round(epsFor(v)) + ')'; } },
    { key: 'zp', label: 'Cation valence', min: 1, max: 3, step: 1, rebuild: true, hide: ['crystal'],
      fmt: function (v) { return '+' + v; } },
    { key: 'zm', label: 'Anion valence', min: 1, max: 3, step: 1, rebuild: true, show: ['salt', 'strong', 'electrode'],
      fmt: function (v) { return '−' + v; } },
    { key: 'salt', label: 'Salt (formula units)', min: 0, max: 500, step: 5, rebuild: true,
      fmt: function (v) { return String(v); } },
    { key: 'nm', label: 'Chain length (monomers)', min: 12, max: 120, step: 6, rebuild: true, show: ['polymer'],
      fmt: function (v) { return String(v); } },
    { key: 'Z', label: 'Colloid charge |Z|', min: 12, max: 120, step: 6, rebuild: true, show: ['colloid'],
      fmt: function (v) { return '−' + v; } },
    { key: 'phi', label: 'Area fraction', min: 0.02, max: 0.45, step: 0.01, rescale: true, hide: ['crystal'],
      fmt: function (v) { return v.toFixed(2); } },
    { key: 'E0', label: 'Applied field E (kT / eσ)', min: -3, max: 3, step: 0.05, hide: ['electrode', 'crystal'],
      fmt: function (v) { return v.toFixed(2); } },
    { key: 'V', label: 'Voltage V (kT/e)', min: -120, max: 120, step: 1, show: ['electrode'],
      fmt: function (v) { return v.toFixed(0) + '  (' + (v * KT_MV / 1000).toFixed(2) + ' V)'; } },
    { key: 'Q', label: 'Wand charge |Q|', min: 1, max: 80, step: 1,
      fmt: function (v) { return String(v); } },
    { key: 'gamma', label: 'Friction γ (D = kT / mγ)', log: logMap(0.05, 20),
      fmt: function (v) { return v.toFixed(2); } },
    { key: 'dt', label: 'Time step', min: 0.0005, max: 0.01, step: 0.0005,
      fmt: function (v) { return v.toFixed(4); } },
    { key: 'spf', label: 'Steps / frame', min: 1, max: 24, step: 1,
      fmt: function (v) { return String(v); } }
  ];

  var OVERLAYS = [
    ['phi', 'Potential φ'],
    ['pairs', 'Potential + ion pairs'],
    ['field', 'Field arrows'],
    ['plain', 'Ions only']
  ];

  var TRACES = {
    ec: { label: 'Coulomb energy per ion (kT)', get: function (o) { return o.eC; } },
    T: { label: 'Kinetic temperature (kT)', get: function (o) { return o.T; }, ref: function () { return 1; } },
    pair: { label: 'Bjerrum-paired fraction', get: function (o) { return o.pairFrac; } },
    J: { label: 'Current density J (e/στ) · dashed: free ions', get: function (o) { return o.Jshort; }, ref: function (o) { return o.geom === 'electrode' ? NaN : o.Jideal; } },
    psi4: { label: 'Square bond order ψ4', get: function (o) { return o.psi4; }, only: ['crystal'] },
    Q: { label: 'Stored charge Q / Ly (e / σ)', get: function (o) { return o.stored; }, only: ['electrode'] },
    rg: { label: 'Radius of gyration Rg (σ)', get: function (o) { return o.rg; }, only: ['polymer'] },
    theta: { label: 'Neutralised fraction θ(2σ)  (dashed: Manning)', get: function (o) { return o.theta2; }, ref: function (o) { return o.manning; }, only: ['polymer'] },
    zeff: { label: 'Z_eff(R + 2.5σ) / |Z|', get: function (o) { return o.zeff2; }, ref: function () { return 0; }, only: ['colloid'] }
  };

  var TEMPLATE = [
    '<div class="cf-stage">',
    '  <div class="cf-panel">',
    '    <canvas class="cf-canvas" role="img" aria-label="Ions in a periodic plane"></canvas>',
    '    <div class="cf-hud">',
    '      <span>N <b data-hud="n">0</b></span>',
    '      <span>t <b data-hud="t">0</b></span>',
    '      <span>T <b data-hud="T">0</b></span>',
    '      <span><b data-hud="fps">0</b> fps</span>',
    '    </div>',
    '    <span class="cf-probe" data-probe></span>',
    '    <span class="cf-caption" data-cap="tool"></span>',
    '  </div>',
    '  <div class="cf-scope">',
    '    <div class="cf-plot"><canvas data-plot="0"></canvas><span class="cf-ptitle" data-ptitle="0"></span><span class="cf-legend" data-legend="0"></span></div>',
    '    <div class="cf-plot"><canvas data-plot="1"></canvas><span class="cf-ptitle" data-ptitle="1"></span><span class="cf-legend" data-legend="1"></span></div>',
    '    <div class="cf-plot"><canvas data-plot="2"></canvas><span class="cf-ptitle" data-ptitle="2"></span><span class="cf-legend" data-legend="2"></span></div>',
    '  </div>',
    '</div>',
    '<div class="cf-bar">',
    '  <button type="button" data-act="run">Pause</button>',
    '  <button type="button" data-act="reset">Reset</button>',
    '  <span class="cf-seg" role="group" aria-label="Pointer tool">',
    '    <button type="button" data-tool="wand+">Wand +</button>',
    '    <button type="button" data-tool="wand-">Wand −</button>',
    '    <button type="button" data-tool="add+">Drop +</button>',
    '    <button type="button" data-tool="add-">Drop −</button>',
    '    <button type="button" data-tool="erase">Erase</button>',
    '  </span>',
    '  <span class="cf-spacer"></span>',
    '  <button type="button" data-act="link">Copy link</button>',
    '  <span class="cf-status" role="status"></span>',
    '</div>',
    '<div class="cf-body">',
    '  <div class="cf-controls">',
    '    <div class="cf-field cf-field--wide"><label><span>Scenario</span></label><select data-sel="sc"></select></div>',
    '    <div class="cf-field"><label><span>Overlay</span></label><select data-sel="overlay"></select></div>',
    '    <div class="cf-field"><label><span>Trace</span></label><select data-sel="trace"></select></div>',
    '  </div>',
    '  <div class="cf-calc">',
    '    <h4>Calculations <span>live</span></h4>',
    '    <div class="cf-grid" data-calc></div>',
    '    <p class="cf-units">Units: σ = ion diameter, kT = 1, e = 1, τ = σ√(m/kT). Physical readings take σ = 0.5 nm and T = 298 K.</p>',
    '  </div>',
    '</div>'
  ].join('\n');

  var TOOL_HINT = {
    'wand+': 'hold to wield a +Q charge',
    'wand-': 'hold to wield a −Q charge',
    'add+': 'click to drop a cation',
    'add-': 'click to drop an anion',
    'erase': 'click an ion to remove it'
  };

  function fmtNum(v, d) {
    if (v === undefined || v === null || !isFinite(v)) return '—';
    return (v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
  }

  function mount(root) {
    if (root.dataset.cfMounted) return;
    root.dataset.cfMounted = '1';
    root.classList.add('cf');
    root.innerHTML = TEMPLATE;

    var canvas = root.querySelector('.cf-canvas');
    var ctx = canvas.getContext('2d', { alpha: false });
    if (!ctx) return;
    var plots = [0, 1, 2].map(function (k) {
      var c = root.querySelector('[data-plot="' + k + '"]');
      return {
        canvas: c, ctx: c.getContext('2d'),
        title: root.querySelector('[data-ptitle="' + k + '"]'),
        legend: root.querySelector('[data-legend="' + k + '"]')
      };
    });

    var conf = { sc: 'salt', Q: 30, overlay: 'phi', trace: 'ec', seed: 20261110 };
    var p0 = PRESETS.salt;
    for (var key in p0) if (Object.prototype.hasOwnProperty.call(p0, key)) conf[key] = p0[key];

    var statusEl = root.querySelector('.cf-status');
    var runBtn = root.querySelector('[data-act="run"]');
    var capEl = root.querySelector('[data-cap="tool"]');
    var probeEl = root.querySelector('[data-probe]');
    var calcEl = root.querySelector('[data-calc]');
    var hud = {
      n: root.querySelector('[data-hud="n"]'),
      t: root.querySelector('[data-hud="t"]'),
      T: root.querySelector('[data-hud="T"]'),
      fps: root.querySelector('[data-hud="fps"]')
    };
    var tool = 'wand+';

    // ---- selects ----------------------------------------------------------

    function fillSelect(sel, items) {
      sel.innerHTML = '';
      items.forEach(function (it) {
        var o = document.createElement('option');
        o.value = it[0];
        o.textContent = it[1];
        sel.appendChild(o);
      });
    }

    var selSc = root.querySelector('[data-sel="sc"]');
    var selOv = root.querySelector('[data-sel="overlay"]');
    var selTr = root.querySelector('[data-sel="trace"]');
    fillSelect(selSc, SCENARIOS.map(function (s) { return [s.id, s.label]; }));
    fillSelect(selOv, OVERLAYS);

    function availableTraces() {
      var out = [];
      for (var k in TRACES) {
        if (!Object.prototype.hasOwnProperty.call(TRACES, k)) continue;
        var tr = TRACES[k];
        if (tr.only && tr.only.indexOf(sys.geom) < 0) continue;
        out.push([k, tr.label]);
      }
      return out;
    }

    selSc.addEventListener('change', function () { loadScenario(selSc.value, true); });
    selOv.addEventListener('change', function () { conf.overlay = selOv.value; });
    selTr.addEventListener('change', function () { conf.trace = selTr.value; trace.length = 0; });

    // ---- sliders ----------------------------------------------------------

    var controlsHost = root.querySelector('.cf-controls');
    var readouts = {};
    var rebuildTimer = 0;
    SLIDERS.forEach(function (ctl) {
      var wrap = document.createElement('div');
      wrap.className = 'cf-field';
      var label = document.createElement('label');
      var name = document.createElement('span');
      name.textContent = ctl.label;
      var value = document.createElement('span');
      label.appendChild(name);
      label.appendChild(value);
      var input = document.createElement('input');
      input.type = 'range';
      if (ctl.log) {
        input.min = '0'; input.max = '1000'; input.step = '1';
      } else {
        input.min = String(ctl.min); input.max = String(ctl.max); input.step = String(ctl.step);
      }
      wrap.appendChild(label);
      wrap.appendChild(input);
      controlsHost.appendChild(wrap);
      readouts[ctl.key] = { input: input, value: value, ctl: ctl, wrap: wrap, name: name };

      input.addEventListener('input', function () {
        var raw = parseFloat(input.value);
        if (!isFinite(raw)) return;
        var v = ctl.log ? ctl.log.toVal(raw) : raw;
        conf[ctl.key] = v;
        value.textContent = ctl.fmt(v);
        if (ctl.rebuild) {
          clearTimeout(rebuildTimer);
          rebuildTimer = setTimeout(function () { rebuild(); }, 180);
        } else if (ctl.rescale) {
          sys.rescale(targetL());
          applyParams();
          sys.resetStats();
        } else {
          applyParams();
          if (ctl.key === 'lB' || ctl.key === 'E0' || ctl.key === 'V' || ctl.key === 'gamma') sys.resetStats();
        }
      });
    });

    function sliderVisible(ctl) {
      if (ctl.show) return ctl.show.indexOf(conf.sc) >= 0;
      if (ctl.hide) return ctl.hide.indexOf(conf.sc) < 0;
      return true;
    }

    function syncControls() {
      SLIDERS.forEach(function (ctl) {
        var r = readouts[ctl.key];
        r.input.value = String(ctl.log ? ctl.log.toPos(conf[ctl.key]) : conf[ctl.key]);
        r.value.textContent = ctl.fmt(conf[ctl.key]);
        r.wrap.hidden = !sliderVisible(ctl);
      });
      var poly = sys.geom === 'polymer' || sys.geom === 'colloid';
      readouts.zp.name.textContent = poly ? 'Counterion valence' : 'Cation valence';
      readouts.salt.name.textContent = poly ? 'Added salt (pairs)' : (conf.sc === 'crystal' ? 'Ion pairs (lattice)' : 'Salt (formula units)');
      readouts.salt.input.min = poly || conf.sc === 'crystal' ? '0' : '10';
      selSc.value = conf.sc;
      selOv.value = conf.overlay;
      fillSelect(selTr, availableTraces());
      var ok = false;
      for (var i = 0; i < selTr.options.length; i++) if (selTr.options[i].value === conf.trace) ok = true;
      if (!ok) conf.trace = selTr.options.length ? selTr.options[0].value : 'ec';
      selTr.value = conf.trace;
    }

    // ---- engine -----------------------------------------------------------

    var sys = new System(conf.seed);
    var running = true;
    var trace = [];
    var frames = 0, lastStamp = 0;
    var potScale = 1;
    var mapBack = document.createElement('canvas');
    var mapImage = null;
    var pointer = { inside: false, x: 0, y: 0, down: false, dragMacro: false };

    function targetL() {
      var area = 0;
      for (var i = 0; i < sys.n; i++) area += sys.kind[i] === 2 ? Math.PI * Math.pow(sys.macroR + 0.5, 2) : Math.PI / 4;
      var L = Math.sqrt(area / clamp(conf.phi, 0.02, 0.45));
      if (sys.geom === 'polymer') L = Math.max(L, 24);
      if (sys.geom === 'colloid') L = Math.max(L, 2 * sys.macroR + 18);
      return Math.max(L, 16);
    }

    function applyParams() {
      sys.lB = conf.lB;
      sys.gamma = conf.gamma;
      sys.dt = conf.dt;
      sys.wand.Q = (tool === 'wand-' ? -1 : 1) * conf.Q;
      sys.E0 = sys.geom === 'electrode' ? conf.V / sys.L : (sys.geom === 'crystal' ? 0 : conf.E0);
    }

    function rebuild() {
      sys.seed = conf.seed >>> 0;
      applyParams();
      sys.build(scenarioById(conf.sc).geom, conf);
      applyParams();
      sys.collect = true;                 // measure the t = 0 state too, so a
      sys.forces();                       // paused Reset shows e.g. the exact
      sys.collect = false;                // Madelung energy of the lattice
      sys.measure();
      trace.length = 0;
      syncControls();
      drawCalc();
    }

    function loadScenario(id, usePreset) {
      conf.sc = scenarioById(id).id;
      if (usePreset) {
        var p = PRESETS[conf.sc];
        for (var k in p) if (Object.prototype.hasOwnProperty.call(p, k)) conf[k] = p[k];
      }
      sys.geom = scenarioById(conf.sc).geom;
      rebuild();
    }

    // ---- canvas helpers ---------------------------------------------------

    function fitCanvas(c) {
      var dpr = Math.min(2, window.devicePixelRatio || 1);
      var rect = c.getBoundingClientRect();
      var w = Math.max(1, Math.round(rect.width * dpr));
      var hgt = Math.max(1, Math.round(rect.height * dpr));
      if (w !== c.width || hgt !== c.height) { c.width = w; c.height = hgt; }
      return dpr;
    }

    function ionColor(qv, k) {
      if (k === 1) return MONO_COL;
      var z = Math.min(3, Math.round(Math.abs(qv))) - 1;
      if (z < 0) return '#888';
      return qv > 0 ? CAT_COL[z] : ANI_COL[z];
    }

    // ---- main view --------------------------------------------------------

    function drawMain() {
      var dpr = fitCanvas(canvas);
      var W = canvas.width, H = canvas.height, L = sys.L, s = W / L;
      var n = sys.n, x = sys.x, y = sys.y, q = sys.q, kind = sys.kind, i;
      var pme = sys.pme, Mx = pme.Mx, My = pme.My;

      ctx.fillStyle = INK;
      ctx.fillRect(0, 0, W, H);

      if (conf.overlay === 'phi' || conf.overlay === 'pairs') {
        var phi = pme.potential();
        if (!mapImage || mapBack.width !== Mx || mapBack.height !== My) {
          mapBack.width = Mx; mapBack.height = My;
          mapImage = mapBack.getContext('2d').createImageData(Mx, My);
        }
        var px = mapImage.data, E0 = sys.geom === 'electrode' ? sys.E0 : 0;
        var rms = 0, mean = 0, S = Mx * My, id;
        for (id = 0; id < S; id++) mean += phi[id] - E0 * (id % Mx) * pme.hx;
        mean /= S;
        for (id = 0; id < S; id++) {
          var v0 = phi[id] - E0 * (id % Mx) * pme.hx - mean;
          rms += v0 * v0;
        }
        rms = Math.sqrt(rms / S);
        potScale += (Math.max(0.25, 2.2 * rms) - potScale) * 0.05;
        for (id = 0; id < S; id++) {
          var v = (phi[id] - E0 * (id % Mx) * pme.hx - mean) / potScale;
          var tt = Math.abs(v) > 6 ? (v > 0 ? 1 : -1) : Math.tanh(v);
          var col = tt > 0 ? POS_RGB : NEG_RGB, a = 0.5 * Math.abs(tt);
          px[4 * id] = INK_RGB[0] + (col[0] - INK_RGB[0]) * a;
          px[4 * id + 1] = INK_RGB[1] + (col[1] - INK_RGB[1]) * a;
          px[4 * id + 2] = INK_RGB[2] + (col[2] - INK_RGB[2]) * a;
          px[4 * id + 3] = 255;
        }
        mapBack.getContext('2d').putImageData(mapImage, 0, 0);
        ctx.imageSmoothingEnabled = true;
        ctx.drawImage(mapBack, 0, 0, W, H);
      }

      if (conf.overlay === 'field') {
        var ex = pme.ex, ey = pme.ey, step = Math.max(2.5, L / 22);
        var gn = Math.floor(L / step), gsp = L / gn;
        ctx.strokeStyle = 'rgba(255,255,255,0.42)';
        ctx.lineWidth = Math.max(1, 1.1 * dpr);
        ctx.beginPath();
        for (var gy = 0; gy < gn; gy++) {
          for (var gx = 0; gx < gn; gx++) {
            var cx = (gx + 0.5) * gsp, cy = (gy + 0.5) * gsp;
            var Ex = pme.sample(ex, cx, cy) + sys.E0, Ey = pme.sample(ey, cx, cy);
            var mag = Math.sqrt(Ex * Ex + Ey * Ey);
            if (mag < 1e-6) continue;
            var len = 0.45 * gsp * Math.min(1, Math.log(1 + mag) / 1.6);
            var ux = Ex / mag, uy = Ey / mag;
            var x0 = (cx - ux * len) * s, y0 = (cy - uy * len) * s;
            var x1 = (cx + ux * len) * s, y1 = (cy + uy * len) * s;
            ctx.moveTo(x0, y0); ctx.lineTo(x1, y1);
            var hl = Math.min(6 * dpr, 0.4 * len * s);
            ctx.moveTo(x1, y1); ctx.lineTo(x1 - hl * (ux - 0.5 * uy), y1 - hl * (uy + 0.5 * ux));
            ctx.moveTo(x1, y1); ctx.lineTo(x1 - hl * (ux + 0.5 * uy), y1 - hl * (uy - 0.5 * ux));
          }
        }
        ctx.stroke();
      }

      if (sys.wall) {
        var ew = Math.max(3, 0.35 * s);
        var Vs = conf.V;
        ctx.fillStyle = Vs >= 0 ? 'rgba(255,93,115,0.85)' : 'rgba(111,157,255,0.85)';
        ctx.fillRect(0, 0, ew, H);
        ctx.fillStyle = Vs >= 0 ? 'rgba(111,157,255,0.85)' : 'rgba(255,93,115,0.85)';
        ctx.fillRect(W - ew, 0, ew, H);
      }

      // chain backbone
      if (sys.nMono > 1) {
        ctx.strokeStyle = 'rgba(181,123,255,0.55)';
        ctx.lineWidth = Math.max(1, 0.3 * s);
        ctx.lineCap = 'round';
        ctx.beginPath();
        for (i = 1; i < sys.nMono; i++) {
          var bx = x[i] - x[i - 1], by = y[i] - y[i - 1];
          if (Math.abs(bx) > L / 2 || Math.abs(by) > L / 2) continue;
          ctx.moveTo(x[i - 1] * s, y[i - 1] * s);
          ctx.lineTo(x[i] * s, y[i] * s);
        }
        ctx.stroke();
      }

      // colloid
      if (sys.macro >= 0) {
        var M = sys.macro, mR = sys.macroR * s, mx = x[M] * s, my = y[M] * s;
        var grd = ctx.createRadialGradient(mx - 0.3 * mR, my - 0.3 * mR, 0.1 * mR, mx, my, mR);
        grd.addColorStop(0, '#3a4f9a');
        grd.addColorStop(1, '#16204a');
        ctx.fillStyle = grd;
        ctx.beginPath();
        ctx.arc(mx, my, mR, 0, TWO_PI);
        ctx.fill();
        ctx.strokeStyle = 'rgba(111,157,255,0.8)';
        ctx.lineWidth = Math.max(1, 1.5 * dpr);
        ctx.stroke();
        ctx.fillStyle = 'rgba(255,255,255,0.8)';
        ctx.font = Math.round(Math.max(10, 0.9 * s) ) + 'px ui-monospace, Menlo, Consolas, monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('−' + Math.abs(q[M]), mx, my);
      }

      // ions, batched by colour
      var groups = {};
      for (i = 0; i < n; i++) {
        if (kind[i] === 2) continue;
        var key = ionColor(q[i], kind[i]);
        (groups[key] || (groups[key] = [])).push(i);
      }
      ctx.strokeStyle = 'rgba(11,13,18,0.9)';
      ctx.lineWidth = Math.max(1, 0.06 * s);
      for (var gk in groups) {
        if (!Object.prototype.hasOwnProperty.call(groups, gk)) continue;
        var list = groups[gk];
        ctx.fillStyle = gk;
        ctx.beginPath();
        for (var li = 0; li < list.length; li++) {
          var j = list[li], rr = 0.5 * s * (1 + 0.1 * (Math.min(3, Math.abs(q[j])) - 1));
          ctx.moveTo(x[j] * s + rr, y[j] * s);
          ctx.arc(x[j] * s, y[j] * s, rr, 0, TWO_PI);
        }
        ctx.fill();
        if (s > 5) ctx.stroke();
      }

      // Bjerrum pair bonds, on top of the discs so contact pairs show
      if (conf.overlay === 'pairs' && sys.nPairs) {
        ctx.strokeStyle = 'rgba(255,255,255,0.7)';
        ctx.lineWidth = Math.max(1, 0.16 * s);
        ctx.lineCap = 'round';
        ctx.beginPath();
        for (var pp = 0; pp < sys.nPairs; pp++) {
          var a = sys.pairs[2 * pp], b = sys.pairs[2 * pp + 1];
          var ddx = x[b] - x[a], ddy = y[b] - y[a];
          if (Math.abs(ddx) > L / 2 || Math.abs(ddy) > L / 2) continue;
          ctx.moveTo(x[a] * s, y[a] * s);
          ctx.lineTo(x[b] * s, y[b] * s);
        }
        ctx.stroke();
      }

      // charge glyphs
      if (s >= 7) {
        ctx.strokeStyle = 'rgba(11,13,18,0.85)';
        ctx.lineWidth = Math.max(1, 0.1 * s);
        ctx.lineCap = 'butt';
        ctx.beginPath();
        var g = 0.22 * s;
        for (i = 0; i < n; i++) {
          if (kind[i] === 2 || q[i] === 0) continue;
          var gx2 = x[i] * s, gy2 = y[i] * s;
          ctx.moveTo(gx2 - g, gy2); ctx.lineTo(gx2 + g, gy2);
          if (q[i] > 0) { ctx.moveTo(gx2, gy2 - g); ctx.lineTo(gx2, gy2 + g); }
        }
        ctx.stroke();
      }

      // wand
      if (sys.wand.on) {
        var wx = sys.wand.x * s, wy = sys.wand.y * s, wr = Math.max(8 * dpr, 1.2 * s);
        var wpos = sys.wand.Q > 0;
        var halo = ctx.createRadialGradient(wx, wy, 0, wx, wy, 3 * wr);
        halo.addColorStop(0, wpos ? 'rgba(255,93,115,0.55)' : 'rgba(111,157,255,0.55)');
        halo.addColorStop(1, 'rgba(0,0,0,0)');
        ctx.fillStyle = halo;
        ctx.beginPath(); ctx.arc(wx, wy, 3 * wr, 0, TWO_PI); ctx.fill();
        ctx.strokeStyle = wpos ? '#ff5d73' : '#6f9dff';
        ctx.lineWidth = Math.max(1.5, 2 * dpr);
        ctx.beginPath(); ctx.arc(wx, wy, wr, 0, TWO_PI); ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(wx - 0.5 * wr, wy); ctx.lineTo(wx + 0.5 * wr, wy);
        if (wpos) { ctx.moveTo(wx, wy - 0.5 * wr); ctx.lineTo(wx, wy + 0.5 * wr); }
        ctx.stroke();
      }
    }

    // ---- plots ------------------------------------------------------------

    function drawPlot(P, spec) {
      var dpr = fitCanvas(P.canvas);
      var c = P.ctx, W = P.canvas.width, H = P.canvas.height;
      c.fillStyle = INK;
      c.fillRect(0, 0, W, H);
      P.title.textContent = spec.title || '';
      P.legend.innerHTML = spec.legend || '';
      var ml = 38 * dpr, mr = 10 * dpr, mt = 24 * dpr, mb = 20 * dpr;
      var pw = W - ml - mr, ph = H - mt - mb;
      if (pw < 10 || ph < 10) return;
      var x0 = spec.xmin, x1 = spec.xmax, y0 = spec.ymin, y1 = spec.ymax;
      if (!(y1 > y0)) { y1 = y0 + 1; }
      var X = function (v) { return ml + (v - x0) / (x1 - x0) * pw; };
      var Y = function (v) { return mt + (1 - (v - y0) / (y1 - y0)) * ph; };

      c.font = Math.round(10 * dpr) + 'px ui-monospace, Menlo, Consolas, monospace';
      c.fillStyle = 'rgba(255,255,255,0.4)';
      c.strokeStyle = 'rgba(255,255,255,0.12)';
      c.lineWidth = 1;
      c.textBaseline = 'middle';
      c.textAlign = 'right';
      var yt = niceTicks(y0, y1, 4);
      c.beginPath();
      yt.forEach(function (v) {
        var yy = Math.round(Y(v)) + 0.5;
        c.moveTo(ml, yy); c.lineTo(ml + pw, yy);
        c.fillText(tickLabel(v, yt.step), ml - 4 * dpr, yy);
      });
      c.stroke();
      c.textAlign = 'center';
      c.textBaseline = 'top';
      var xt = niceTicks(x0, x1, 5);
      xt.forEach(function (v) {
        c.fillText(tickLabel(v, xt.step), X(v), mt + ph + 4 * dpr);
      });
      if (spec.xlabel) {
        c.textAlign = 'right';
        c.textBaseline = 'bottom';
        c.fillText(spec.xlabel, ml + pw, mt + ph - 2 * dpr);
      }

      c.save();
      c.beginPath();
      c.rect(ml, mt, pw, ph);
      c.clip();
      (spec.hlines || []).forEach(function (hl) {
        if (!isFinite(hl.y)) return;
        c.strokeStyle = hl.color || 'rgba(255,255,255,0.5)';
        c.setLineDash([4 * dpr, 4 * dpr]);
        c.lineWidth = 1.2 * dpr;
        c.beginPath();
        c.moveTo(ml, Y(hl.y)); c.lineTo(ml + pw, Y(hl.y));
        c.stroke();
      });
      c.setLineDash([]);
      (spec.series || []).forEach(function (sr) {
        var xs = sr.x, ys = sr.y, k;
        c.strokeStyle = sr.color;
        c.fillStyle = sr.color;
        c.lineWidth = (sr.width || 1.5) * dpr;
        c.setLineDash(sr.dash ? [5 * dpr, 4 * dpr] : []);
        if (sr.dots) {
          for (k = 0; k < xs.length; k++) {
            if (!isFinite(ys[k])) continue;
            c.beginPath();
            if (ys[k] > y1) {                       // off scale: caret at the top
              var cx0 = X(xs[k]), cy0 = mt + 2 * dpr;
              c.moveTo(cx0, cy0); c.lineTo(cx0 - 4 * dpr, cy0 + 7 * dpr); c.lineTo(cx0 + 4 * dpr, cy0 + 7 * dpr);
              c.closePath();
            } else {
              c.arc(X(xs[k]), Y(ys[k]), 1.8 * dpr, 0, TWO_PI);
            }
            c.fill();
          }
        } else {
          c.beginPath();
          var started = false;
          for (k = 0; k < xs.length; k++) {
            if (!isFinite(ys[k])) { started = false; continue; }
            if (!started) { c.moveTo(X(xs[k]), Y(ys[k])); started = true; }
            else c.lineTo(X(xs[k]), Y(ys[k]));
          }
          c.stroke();
        }
      });
      c.restore();
    }

    function niceTicks(lo, hi, count) {
      var span = hi - lo;
      if (!(span > 0)) return [lo];
      var raw = span / count, mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
      var f = raw / mag, stp = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
      var out = [];
      for (var v = Math.ceil(lo / stp) * stp; v <= hi + 1e-9 * span; v += stp) out.push(Math.abs(v) < 1e-12 * stp ? 0 : v);
      out.step = stp;
      return out;
    }

    /* Decimals follow the tick step, so neighbouring labels never collide. */
    function tickLabel(v, step) {
      if (v === 0) return '0';
      var d = step ? clamp(-Math.floor(Math.log(step) / Math.LN10 + 1e-9), 0, 5) : 2;
      return v.toFixed(d);
    }

    function lg(color, text, dashed) {
      return '<i style="' + (dashed ? 'border-top:2px dashed ' + color + ';background:none' : 'background:' + color) + '"></i>' + text;
    }

    function rangeOf(arrs, padFrac, floor0) {
      var lo = Infinity, hi = -Infinity;
      arrs.forEach(function (a) {
        for (var i = 0; i < a.length; i++) if (isFinite(a[i])) { if (a[i] < lo) lo = a[i]; if (a[i] > hi) hi = a[i]; }
      });
      if (!isFinite(lo)) { lo = 0; hi = 1; }
      if (floor0) lo = Math.min(lo, 0);
      var pad = (hi - lo) * (padFrac || 0.08) || 0.5;
      return [lo - (floor0 && lo === 0 ? 0 : pad), hi + pad];
    }

    function drawPlots() {
      var o = sys.obs, i;
      // --- plot 0: g(r)
      var rs = [], gpm = [], gpp = [], gmm = [], dr = RDF_MAX / RDF_NB, ga = sys.rdfAvg;
      for (i = 0; i < RDF_NB; i++) {
        rs.push((i + 0.5) * dr);
        gpm.push(ga[i]); gpp.push(ga[RDF_NB + i]); gmm.push(ga[2 * RDF_NB + i]);
      }
      var gmax = 0;
      for (i = 0; i < RDF_NB; i++) gmax = Math.max(gmax, gpm[i], gpp[i], gmm[i]);
      drawPlot(plots[0], {
        title: 'g(r) — pair distribution',
        legend: lg('#b57bff', '+− ') + lg(CAT_COL[0], '++ ') + lg(ANI_COL[0], '−−'),
        xmin: 0, xmax: RDF_MAX, ymin: 0, ymax: Math.max(2, Math.min(gmax * 1.1, 60)),
        xlabel: 'r / σ',
        hlines: [{ y: 1, color: 'rgba(255,255,255,0.3)' }],
        series: [
          { x: rs, y: gpp, color: CAT_COL[0] },
          { x: rs, y: gmm, color: ANI_COL[0] },
          { x: rs, y: gpm, color: '#b57bff', width: 1.8 }
        ]
      });

      // --- plot 1: trace
      var tr = TRACES[conf.trace] || TRACES.ec;
      var tx = [], ty = [], tref = [];
      for (i = 0; i < trace.length; i++) { tx.push(trace[i][0]); ty.push(trace[i][1]); tref.push(trace[i][2]); }
      var tr0 = tx.length ? tx[0] : 0, tr1 = tx.length ? tx[tx.length - 1] : 1;
      if (tr1 - tr0 < 1) tr1 = tr0 + 1;
      var yr = rangeOf([ty, tref], 0.15, false);
      var last = ty.length ? ty[ty.length - 1] : NaN;
      drawPlot(plots[1], {
        title: tr.label,
        legend: '<b>' + fmtNum(last, Math.abs(last) < 10 ? 3 : 1) + '</b>',
        xmin: tr0, xmax: tr1, ymin: yr[0], ymax: yr[1], xlabel: 't / τ',
        series: [
          { x: tx, y: tref, color: 'rgba(255,255,255,0.55)', dash: true, width: 1.2 },
          { x: tx, y: ty, color: '#f6e7b2', width: 1.6 }
        ]
      });

      // --- plot 2: scenario-specific
      var P = plots[2], geom = sys.geom;
      if (geom === 'bulk' || geom === 'crystal') {
        var k1 = TWO_PI / sys.L, ka = [], sa = [], kd = [], sd = [], kth = [], sth = [];
        for (i = 1; i <= sys.nkMax; i++) {
          ka.push(i * k1); sa.push(sys.skAxis[i]);
          if (i * k1 * Math.SQRT2 <= SK_KMAX) { kd.push(i * k1 * Math.SQRT2); sd.push(sys.skDiag[i]); }
        }
        for (i = 0; i <= 120; i++) { var kk = i / 120 * SK_KMAX; kth.push(kk); sth.push(kk / (kk + o.kappa)); }
        var smax = 1.4, peak = 0, kpk = 0, cap = 5;
        for (i = 0; i < sa.length; i++) { smax = Math.max(smax, sa[i] * 1.1); if (sa[i] > peak) { peak = sa[i]; kpk = ka[i]; } }
        for (i = 0; i < sd.length; i++) { smax = Math.max(smax, sd[i] * 1.1); if (sd[i] > peak) { peak = sd[i]; kpk = kd[i]; } }
        var clipped = smax > cap;
        drawPlot(P, {
          title: 'Charge structure factor S_ZZ(k) / ⟨q²⟩' + (clipped ? '  · peak ' + peak.toFixed(0) + ' at kσ = ' + kpk.toFixed(2) + ' (▲)' : ''),
          legend: lg('#f6e7b2', 'axes ') + lg('#3fc8d8', 'diagonals ') + lg('rgba(255,255,255,0.6)', 'Debye–Hückel', true),
          xmin: 0, xmax: SK_KMAX, ymin: 0, ymax: Math.min(smax, cap), xlabel: 'k σ',
          hlines: [{ y: 1, color: 'rgba(255,255,255,0.25)' }],
          series: [
            { x: kth, y: sth, color: 'rgba(255,255,255,0.6)', dash: true, width: 1.2 },
            { x: kd, y: sd, color: '#3fc8d8', dots: true },
            { x: ka, y: sa, color: '#f6e7b2', dots: true }
          ]
        });
      } else if (geom === 'electrode') {
        var xb = [], rp = [], rm = [], dxb = sys.L / PROF_NB;
        for (i = 0; i < PROF_NB; i++) { xb.push((i + 0.5) * dxb); rp.push(sys.profP[i]); rm.push(sys.profM[i]); }
        var xp = [], php = [], Mx = sys.pme.Mx, pmax = 0;
        var mid = sys.profPhi[Mx >> 1] || 0;
        for (i = 0; i < Mx; i++) { xp.push(i * sys.pme.hx); php.push(sys.profPhi[i] - mid); pmax = Math.max(pmax, Math.abs(php[i])); }
        var dmax = 0;
        for (i = 0; i < PROF_NB; i++) dmax = Math.max(dmax, rp[i], rm[i]);
        // potential rescaled onto the density axis, labelled in the legend
        var phs = pmax > 1e-9 ? dmax / pmax : 1, phScaled = php.map(function (v) { return v * phs * 0.9; });
        drawPlot(P, {
          title: 'Profiles across the cell',
          legend: lg(CAT_COL[0], 'ρ₊ ') + lg(ANI_COL[0], 'ρ₋ ') + lg('#f6e7b2', 'φ (±' + pmax.toFixed(1) + ' kT/e)'),
          xmin: 0, xmax: sys.L, ymin: -dmax * 1.0 - 1e-6, ymax: dmax * 1.05 + 1e-6, xlabel: 'x / σ',
          hlines: [{ y: 0, color: 'rgba(255,255,255,0.25)' }],
          series: [
            { x: xb, y: rp, color: CAT_COL[0] },
            { x: xb, y: rm, color: ANI_COL[0] },
            { x: xp, y: phScaled, color: '#f6e7b2', width: 1.4 }
          ]
        });
      } else if (geom === 'polymer') {
        var dd = [], th = [], db = THETA_MAX / THETA_NB;
        for (i = 0; i < THETA_NB; i++) { dd.push((i + 1) * db); th.push(sys.theta[i]); }
        var tr2 = rangeOf([th, [o.manning || 0, 1]], 0.05, true);
        drawPlot(P, {
          title: 'Chain charge neutralised within d',
          legend: lg('#f6e7b2', 'θ(d) ') + lg('rgba(255,255,255,0.6)', 'Manning 1−1/zξ', true),
          xmin: 0, xmax: THETA_MAX, ymin: tr2[0], ymax: tr2[1], xlabel: 'd / σ',
          hlines: [{ y: o.manning, color: 'rgba(255,255,255,0.6)' }, { y: 1, color: 'rgba(255,255,255,0.2)' }],
          series: [{ x: dd, y: th, color: '#f6e7b2', width: 1.8 }]
        });
      } else if (geom === 'colloid') {
        var zr = [], zz = [], zb = ZEFF_SPAN / ZEFF_NB;
        for (i = 0; i < ZEFF_NB; i++) { zr.push(sys.macroR + (i + 1) * zb); zz.push(sys.zeff[i]); }
        var zrg = rangeOf([zz, [0]], 0.08, false);
        drawPlot(P, {
          title: 'Integrated charge Z_eff(r) / |Z|',
          legend: lg('#f6e7b2', 'Z_eff ') + '<span class="cf-note">&gt; 0 = charge inversion</span>',
          xmin: sys.macroR, xmax: sys.macroR + ZEFF_SPAN, ymin: Math.max(-1.05, zrg[0]), ymax: Math.max(0.15, zrg[1]),
          xlabel: 'r / σ',
          hlines: [{ y: 0, color: 'rgba(255,255,255,0.55)' }],
          series: [{ x: zr, y: zz, color: '#f6e7b2', width: 1.8 }]
        });
      }
    }

    // ---- calculations table ----------------------------------------------

    function row(label, value, note) {
      // a div grid, not a <table>: the theme styles every table in posts
      return '<span class="cf-k">' + label + '</span><span class="cf-v">' + value +
        '</span><span class="cf-n">' + (note || '') + '</span>';
    }

    function drawCalc() {
      var o = sys.obs, h = '', geom = sys.geom, lB = sys.lB;
      var comp = o.nPos + ' cations, ' + o.nNeg + ' anions' + (sys.nMono ? ', ' + sys.nMono + ' monomers' : '');
      h += row('Particles', String(sys.n), comp);
      h += row('Box', fmtNum(sys.L, 1) + ' σ', 'PME mesh ' + sys.pme.Mx + '², α = ' + sys.alpha.toFixed(2) + ' σ⁻¹, rᴄ = ' + sys.rc.toFixed(1) + ' σ');
      h += row('Bjerrum length ℓB', fmtNum(lB, 2) + ' σ', '= ' + (lB * SIGMA_NM).toFixed(2) + ' nm; εr ≈ ' + Math.round(epsFor(lB)) + ' at 298 K');
      h += row('Coupling Γ', fmtNum(o.Gamma, 2), 'ℓB⟨q²⟩ / a, a = ' + fmtNum(o.aws, 2) + ' σ (Wigner–Seitz)');
      h += row('Screening length 1/κ', o.kappa > 0 ? (1 / o.kappa < 0.1 ? (1 / o.kappa).toPrecision(2) : fmtNum(1 / o.kappa, 2)) + ' σ' : '—', 'κ = 2πℓB Σn q² (2-D Debye–Hückel)');
      h += row('Temperature', fmtNum(o.T, 3) + ' kT', Math.abs(sys.E0) > 1e-4 ?
        'kinetic, drift removed; a strong field heats the ions through collisions' : 'kinetic; thermostat at 1');
      h += row('Coulomb energy / ion', fmtNum(o.eC, 3) + ' kT', 'real ' + fmtNum(sys.eReal, 0) + ' + mesh ' + fmtNum(sys.eLong, 0) + ' + self ' + fmtNum(sys.eSelf, 0));
      var rB = lB * Math.max(1, conf.zp) * (geom === 'bulk' || geom === 'electrode' ? Math.max(1, conf.zm) : 1);
      h += row('Bjerrum-paired ions', (100 * o.pairFrac).toFixed(1) + ' %', rB > 1 ?
        '+/− closer than ' + (rB < 2 ? 'ℓB|q₊q₋| = ' + rB.toFixed(2) + ' σ' : '2 σ (ℓB|q₊q₋| = ' + rB.toFixed(1) + ' σ, capped)') :
        'none possible: ℓB|q₊q₋| < σ');
      if (Math.abs(o.netQ) > 1e-9) h += row('Net charge', fmtNum(o.netQ, 0) + ' e', 'neutralised by a uniform background');

      if (geom !== 'crystal' && Math.abs(sys.E0) > 1e-4) {
        h += row('Current density J', fmtNum(o.Jlong, 4), 'e / (σ τ), time-averaged');
        if (geom !== 'electrode') {
          h += row('Conductivity Λ / Λₙₑ', fmtNum(o.lamRatio, 3), 'measured / Nernst–Einstein Σq²D/kTA');
        }
      }
      if (geom === 'crystal') {
        h += row('Eᴄ per ion ÷ (ℓB/a)', fmtNum(o.madelung, 4), 'perfect lattice: −M/2 = ' + (-MADELUNG_2D / 2).toFixed(4));
        h += row('Bond order ψ₄', fmtNum(o.psi4, 3), '1 = square crystal, ≈ 0 = liquid');
      }
      if (geom === 'electrode') {
        h += row('Applied voltage', fmtNum(o.V, 1) + ' kT/e', '= ' + (o.V * KT_MV / 1000).toFixed(3) + ' V across one period');
        h += row('Bulk field / applied', fmtNum(o.screen, 3), 'middle third of the cell; 0 = fully screened');
        h += row('Stored charge Q / Ly', fmtNum(o.stored, 3) + ' e/σ', 'net ion charge in the left half');
        h += row('Capacitance per length', fmtNum(o.cap, 4), '(e/σ) / (kT/e)');
      }
      if (geom === 'polymer') {
        h += row('Radius of gyration', fmtNum(o.rg, 2) + ' σ', 'end-to-end ' + fmtNum(o.re, 2) + ' σ');
        h += row('Manning parameter ξ', fmtNum(o.xi, 3), 'ℓB / b, b = ' + fmtNum(o.bond, 3) + ' σ');
        h += row('Neutralised within 2σ', fmtNum(o.theta2, 3), 'Manning: 1 − 1/(zξ) = ' + fmtNum(o.manning, 3) + ' (z = ' + o.zc + ')');
      }
      if (geom === 'colloid') {
        h += row('Colloid', 'Z = ' + fmtNum(o.Z, 0), 'radius ' + sys.macroR + ' σ, pinned (drag it)');
        h += row('Z_eff(R + 1.5σ) / |Z|', fmtNum(o.zeff1, 3), 'integrated charge, colloid included');
        h += row('Max Z_eff / |Z|', fmtNum(o.zeffMax, 3), o.zeffMax > 0.02 ? 'above zero: overcharged (charge inversion)' : 'below zero: not overcharged');
      }
      calcEl.innerHTML = h;
    }

    // ---- HUD + probe ------------------------------------------------------

    function updateHUD(stamp) {
      var o = sys.obs;
      hud.n.textContent = String(sys.n);
      hud.t.textContent = sys.t.toFixed(1);
      hud.T.textContent = fmtNum(o.T, 2);
      frames++;
      if (stamp - lastStamp > 500) {
        hud.fps.textContent = String(Math.round(frames * 1000 / (stamp - lastStamp)));
        frames = 0;
        lastStamp = stamp;
      }
      if (pointer.inside) {
        var pme = sys.pme, phi = pme.potential();
        var pv = pme.sample(phi, pointer.x, pointer.y) - (sys.geom === 'electrode' ? sys.E0 * pointer.x : 0);
        var Ex = pme.sample(pme.ex, pointer.x, pointer.y) + sys.E0, Ey = pme.sample(pme.ey, pointer.x, pointer.y);
        probeEl.textContent = 'φ ' + fmtNum(pv, 2) + ' kT/e (' + fmtNum(pv * KT_MV, 0) + ' mV)  |E| ' + Math.sqrt(Ex * Ex + Ey * Ey).toFixed(2) + ' kT/eσ';
      } else {
        probeEl.textContent = '';
      }
    }

    // ---- loop -------------------------------------------------------------

    var statusTimer = 0;
    function say(msg) {
      statusEl.textContent = msg;
      clearTimeout(statusTimer);
      statusTimer = setTimeout(function () { statusEl.textContent = ''; }, 3000);
    }

    var calcTick = 0;
    function frame(stamp) {
      if (running) {
        applyParams();
        var ok = true, spf = Math.round(conf.spf);
        for (var s = 0; s < spf && ok; s++) {
          sys.collect = s === spf - 1;
          ok = sys.step();
        }
        sys.collect = false;
        if (!ok) {
          rebuild();
          say('Numerical blow-up — rebuilt (lower dt or ℓB)');
        } else {
          sys.measure();
          if (sys.obs.capHits > 0) say('Speed cap active — overlaps or dt too large');
          var tr = TRACES[conf.trace] || TRACES.ec;
          trace.push([sys.t, tr.get(sys.obs), tr.ref ? tr.ref(sys.obs) : NaN]);
          if (trace.length > TRACE_LEN) trace.shift();
        }
      }
      drawMain();
      drawPlots();
      updateHUD(stamp || 0);
      if (calcTick++ % 12 === 0 || !running) drawCalc();
      requestAnimationFrame(frame);
    }

    // ---- interaction ------------------------------------------------------

    function toSim(ev) {
      var rect = canvas.getBoundingClientRect();
      return [
        clamp((ev.clientX - rect.left) / rect.width, 0, 0.99999) * sys.L,
        clamp((ev.clientY - rect.top) / rect.height, 0, 0.99999) * sys.L
      ];
    }

    function setTool(t) {
      tool = t;
      var btns = root.querySelectorAll('[data-tool]');
      for (var i = 0; i < btns.length; i++) btns[i].setAttribute('aria-pressed', btns[i].dataset.tool === t ? 'true' : 'false');
      capEl.textContent = TOOL_HINT[t] + (sys.macro >= 0 ? ' · drag the colloid' : '');
      applyParams();
    }

    canvas.addEventListener('pointerdown', function (ev) {
      ev.preventDefault();
      var p = toSim(ev);
      pointer.down = true;
      pointer.x = p[0]; pointer.y = p[1];
      if (canvas.setPointerCapture) { try { canvas.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ } }
      if (sys.macro >= 0) {
        var dx = p[0] - sys.x[sys.macro], dy = p[1] - sys.y[sys.macro];
        if (dx * dx + dy * dy < Math.pow(sys.macroR + 0.5, 2)) { pointer.dragMacro = true; return; }
      }
      if (tool === 'wand+' || tool === 'wand-') {
        sys.wand.on = true;
        sys.wand.x = p[0]; sys.wand.y = p[1];
        applyParams();
      } else if (tool === 'add+' || tool === 'add-') {
        var poly = sys.geom === 'polymer' || sys.geom === 'colloid';
        var zq = tool === 'add+' ? Math.round(conf.zp) : -(poly || sys.geom === 'crystal' ? 1 : Math.round(conf.zm));
        sys.addIon(p[0], p[1], zq);
      } else if (tool === 'erase') {
        if (!sys.removeNear(p[0], p[1], 1.5)) say('No ion under the pointer');
      }
    });

    canvas.addEventListener('pointermove', function (ev) {
      var p = toSim(ev);
      pointer.inside = true;
      pointer.x = p[0]; pointer.y = p[1];
      if (pointer.dragMacro && sys.macro >= 0) {
        sys.x[sys.macro] = sys.wrap(p[0]);
        sys.y[sys.macro] = sys.wrap(p[1]);
      } else if (sys.wand.on) {
        sys.wand.x = p[0]; sys.wand.y = p[1];
      }
    });

    function release() {
      pointer.down = false;
      pointer.dragMacro = false;
      sys.wand.on = false;
    }
    canvas.addEventListener('pointerup', release);
    canvas.addEventListener('pointercancel', release);
    canvas.addEventListener('pointerleave', function () { pointer.inside = false; if (!pointer.down) release(); });

    root.addEventListener('click', function (ev) {
      var btn = ev.target.closest ? ev.target.closest('[data-act],[data-tool]') : null;
      if (!btn || !root.contains(btn)) return;
      if (btn.dataset.tool) { setTool(btn.dataset.tool); return; }
      var act = btn.dataset.act;
      if (act === 'run') {
        running = !running;
        runBtn.textContent = running ? 'Pause' : 'Run';
      } else if (act === 'reset') {
        conf.seed = (conf.seed + 1) >>> 0;
        rebuild();
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

    var URLKEYS = [
      ['lb', 'lB', 0.1, 60], ['zp', 'zp', 1, 3], ['zm', 'zm', 1, 3], ['s', 'salt', 0, 500],
      ['nm', 'nm', 12, 120], ['z', 'Z', 12, 120], ['ph', 'phi', 0.02, 0.45], ['e', 'E0', -3, 3],
      ['v', 'V', -120, 120], ['q', 'Q', 1, 80], ['gm', 'gamma', 0.05, 20], ['dt', 'dt', 0.0005, 0.01],
      ['sf', 'spf', 1, 24]
    ];

    function readHash() {
      var raw = (location.hash || '').replace(/^#/, '');
      if (raw.indexOf('cf=1') === -1) return false;
      var p;
      try { p = new URLSearchParams(raw); } catch (e) { return false; }
      if (p.has('sc')) {
        var sc = scenarioById(p.get('sc'));
        conf.sc = sc.id;
        var pr = PRESETS[sc.id];
        for (var k in pr) if (Object.prototype.hasOwnProperty.call(pr, k)) conf[k] = pr[k];
      }
      URLKEYS.forEach(function (u) {
        if (!p.has(u[0])) return;
        var v = parseFloat(p.get(u[0]));
        if (isFinite(v)) conf[u[1]] = clamp(v, u[2], u[3]);
      });
      ['zp', 'zm', 'salt', 'nm', 'Z', 'spf'].forEach(function (k2) { conf[k2] = Math.round(conf[k2]); });
      if (p.has('ov')) OVERLAYS.forEach(function (ov) { if (ov[0] === p.get('ov')) conf.overlay = ov[0]; });
      if (p.has('tr') && TRACES[p.get('tr')]) conf.trace = p.get('tr');
      if (p.has('sd')) {
        var seed = parseInt(p.get('sd'), 10);
        if (isFinite(seed)) conf.seed = seed >>> 0;
      }
      return true;
    }

    function buildURL() {
      var p = new URLSearchParams();
      p.set('cf', '1');
      p.set('sc', conf.sc);
      URLKEYS.forEach(function (u) {
        var v = conf[u[1]];
        p.set(u[0], u[1] === 'dt' ? v.toFixed(4) : (Math.round(v * 1000) / 1000).toString());
      });
      p.set('ov', conf.overlay);
      p.set('tr', conf.trace);
      p.set('sd', String(conf.seed));
      return location.origin + location.pathname + location.search + '#' + p.toString();
    }

    // ---- boot -------------------------------------------------------------

    readHash();
    sys.seed = conf.seed;
    sys.geom = scenarioById(conf.sc).geom;
    rebuild();
    setTool('wand+');

    requestAnimationFrame(frame);
    root._coulombFluids = { getSystem: function () { return sys; }, conf: conf };
  }

  // -------------------------------------------------------------------------

  var api = {
    System: System, PME: PME, FFT: FFT, mulberry32: mulberry32, erfc: erfc,
    SCENARIOS: SCENARIOS, PRESETS: PRESETS, MADELUNG_2D: MADELUNG_2D
  };
  if (typeof window !== 'undefined') window.CoulombFluids = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  if (typeof document !== 'undefined') {
    var boot = function () {
      var nodes = document.querySelectorAll('[data-sim="coulomb-fluids"]');
      for (var i = 0; i < nodes.length; i++) mount(nodes[i]);
    };
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
    else boot();
  }
})();
