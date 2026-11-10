// Correctness tests: FFT, PME vs an independent Ewald sum, Madelung constant,
// NVE energy conservation, thermostat + stability for every scenario.
//   node test_core.js
// Headless correctness tests for coulomb-fluids.js
const CF = require(require('path').join(__dirname, '../../static/sims/coulomb-fluids/coulomb-fluids.js'));
const { System, PME, FFT, mulberry32, erfc } = CF;
let fails = 0;
function check(name, ok, info) { console.log((ok ? 'PASS ' : 'FAIL ') + name + (info ? '  ' + info : '')); if (!ok) fails++; }

// ---- 1. FFT vs naive DFT
{
  const n = 16, f = new FFT(n), rnd = mulberry32(7);
  const re = new Float64Array(n), im = new Float64Array(n);
  for (let i = 0; i < n; i++) { re[i] = rnd() - 0.5; im[i] = rnd() - 0.5; }
  const r0 = Float64Array.from(re), i0 = Float64Array.from(im);
  f.transform(re, im, 0, 1, -1);
  let err = 0;
  for (let k = 0; k < n; k++) {
    let sr = 0, si = 0;
    for (let m = 0; m < n; m++) { const a = -2 * Math.PI * k * m / n; sr += r0[m] * Math.cos(a) - i0[m] * Math.sin(a); si += r0[m] * Math.sin(a) + i0[m] * Math.cos(a); }
    err = Math.max(err, Math.abs(sr - re[k]), Math.abs(si - im[k]));
  }
  f.transform(re, im, 0, 1, 1);
  let rt = 0; for (let i = 0; i < n; i++) rt = Math.max(rt, Math.abs(re[i] / n - r0[i]), Math.abs(im[i] / n - i0[i]));
  check('FFT matches naive DFT', err < 1e-12 && rt < 1e-12, 'err=' + err.toExponential(2) + ' roundtrip=' + rt.toExponential(2));
}

const refEwald = require('./ref_ewald.js');
// ---- 2. PME vs reference on a random neutral electrolyte (charges only, no WCA)
for (const L of [16, 24, 40]) {
  const sys = new System(3);
  sys.lB = 1.0;
  sys.ensure(200);
  sys.setBox(L);
  const rnd = mulberry32(11 + L);
  const N = Math.round(0.15 * L * L / 0.785 / 4) * 2;
  for (let i = 0; i < N; i++) {
    let px, py, ok = false;
    while (!ok) { px = rnd() * L; py = rnd() * L; ok = true;
      for (let j = 0; j < sys.n; j++) { let dx = px - sys.x[j], dy = py - sys.y[j]; dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); if (dx * dx + dy * dy < 0.81) { ok = false; break; } } }
    sys.add(px, py, (i % 2 ? -1 : 1) * (i % 3 === 0 ? 2 : 1), 0);
  }
  // make neutral
  let Q = 0; for (let i = 0; i < sys.n; i++) Q += sys.q[i];
  sys.q[sys.n - 1] -= Q;
  // remove WCA influence: spread out overlapping pairs is irrelevant -- compare Coulomb only
  const WCA = sys.forces.bind(sys);
  // compute coulomb-only forces: temporarily push WCA cutoff out of play by scaling positions? simpler: subtract WCA part manually
  sys.forces();
  const fx = Float64Array.from(sys.fx.subarray(0, sys.n)), fy = Float64Array.from(sys.fy.subarray(0, sys.n));
  // WCA part
  const n = sys.n, h = L / 2;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    let dx = sys.x[j] - sys.x[i], dy = sys.y[j] - sys.y[i];
    if (dx > h) dx -= L; else if (dx < -h) dx += L; if (dy > h) dy -= L; else if (dy < -h) dy += L;
    const r2 = dx * dx + dy * dy; if (r2 >= Math.pow(2, 1 / 3)) continue;
    const ir2 = 1 / r2, s6 = ir2 ** 3, fr = 24 * s6 * (2 * s6 - 1) * ir2;
    fx[i] += fr * dx; fy[i] += fr * dy; fx[j] -= fr * dx; fy[j] -= fr * dy;
  }
  const ref = refEwald(n, sys.x, sys.y, sys.q, L, 1.0);
  let num = 0, den = 0, mx = 0;
  for (let i = 0; i < n; i++) {
    const ex = fx[i] - ref.fx[i], ey = fy[i] - ref.fy[i];
    num += ex * ex + ey * ey; den += ref.fx[i] ** 2 + ref.fy[i] ** 2; mx = Math.max(mx, Math.hypot(ex, ey));
  }
  const Epme = sys.eReal + sys.eLong + sys.eSelf;
  let sq2 = 0; for (let i = 0; i < n; i++) sq2 += sys.q[i] * sys.q[i];
  const relF = Math.sqrt(num / den), relE = Math.abs(Epme - ref.E) / Math.abs(ref.E);
  check('PME vs Ewald L=' + L + ' N=' + n + ' M=' + sys.pme.Mx + ' alpha=' + sys.alpha.toFixed(3) + ' rc=' + sys.rc.toFixed(2),
    relF < 1e-3 && Math.abs(Epme - ref.E) / sq2 < 1e-4, 'dE/sumq2=' + (Math.abs(Epme - ref.E) / sq2).toExponential(1) + ' rmsF/rmsF=' + relF.toExponential(2) + ' maxErr=' + mx.toExponential(2) + ' E=' + Epme.toFixed(5) + ' ref=' + ref.E.toFixed(5) + ' relE=' + relE.toExponential(2));
}

// ---- 3. Madelung: perfect square NaCl net
{
  const sys = new System(1);
  sys.lB = 1.0;
  sys.build('crystal', { salt: 200 });
  sys.lB = 1.0; sys.forces();
  const N = sys.n, a = sys.latticeA;
  const Ec = (sys.eReal + sys.eLong + sys.eSelf) / N;
  const want = -CF.MADELUNG_2D / 2 / a;
  check('Madelung constant (2-D square, 1/r)', Math.abs(Ec / want - 1) < 1e-3, 'E/N=' + Ec.toFixed(6) + ' want=' + want.toFixed(6) + ' M_measured=' + (-2 * Ec * a).toFixed(6));
  // forces on perfect lattice should vanish by symmetry
  let fmax = 0; for (let i = 0; i < N; i++) fmax = Math.max(fmax, Math.hypot(sys.fx[i], sys.fy[i]));
  check('Perfect lattice is force-free', fmax < 1e-3, 'max|F|=' + fmax.toExponential(2));
}

// ---- 4. Energy conservation at gamma = 0
{
  const sys = new System(5);
  sys.lB = 1.43; sys.gamma = 0; sys.dt = 0.002;
  sys.build('bulk', { zp: 1, zm: 1, salt: 100, phi: 0.15 });
  sys.gamma = 1;
  for (let s = 0; s < 500; s++) sys.step();     // equilibrate first
  sys.gamma = 0;
  const etot = () => {
    let ke = 0; for (let i = 0; i < sys.n; i++) ke += 0.5 * (sys.vx[i] ** 2 + sys.vy[i] ** 2);
    return ke + sys.eReal + sys.eLong + sys.eSelf + sys.eBg + sys.eWca + sys.eBond + sys.eWall;
  };
  const E0 = etot(); let emin = E0, emax = E0;
  const t0 = Date.now();
  for (let s = 0; s < 4000; s++) { sys.step(); const e = etot(); emin = Math.min(emin, e); emax = Math.max(emax, e); }
  const ms = (Date.now() - t0) / 4000;
  const E1 = etot();
  let ke = 0; for (let i = 0; i < sys.n; i++) ke += 0.5 * (sys.vx[i] ** 2 + sys.vy[i] ** 2);
  check('NVE energy conservation (N=' + sys.n + ', dt=0.002, 4000 steps)', Math.abs(E1 - E0) / ke < 0.02,
    'E0=' + E0.toFixed(3) + ' E1=' + E1.toFixed(3) + ' spread=' + (emax - emin).toFixed(4) + ' KE=' + ke.toFixed(1) + '  ' + ms.toFixed(3) + ' ms/step');
}

// ---- 5. Thermostat + timing per scenario + no NaN
for (const sc of CF.SCENARIOS) {
  const p = CF.PRESETS[sc.id];
  const sys = new System(9);
  sys.lB = p.lB; sys.gamma = p.gamma; sys.dt = p.dt;
  sys.build(sc.geom, p);
  sys.E0 = sc.geom === 'electrode' ? p.V / sys.L : p.E0;
  let ok = true, Tsum = 0, cnt = 0;
  const t0 = Date.now();
  const steps = sc.geom === 'polymer' ? 25000 : 5000;
  for (let s = 0; s < steps && ok; s++) {
    sys.collect = s % 6 === 5;
    ok = sys.step();
    if (sys.collect) { sys.measure(); if (s > steps / 2) { Tsum += sys.obs.T; cnt++; } }
  }
  const ms = (Date.now() - t0) / steps;
  const o = sys.obs;
  check('scenario ' + sc.id + ' N=' + sys.n + ' L=' + sys.L.toFixed(1) + ' M=' + sys.pme.Mx, ok && Math.abs(Tsum / cnt - 1) < 0.1,
    'T=' + (Tsum / cnt).toFixed(3) + ' eC=' + (o.eC || 0).toFixed(3) + ' pair=' + (o.pairFrac || 0).toFixed(3) + ' ' + ms.toFixed(3) + ' ms/step capHits=' + o.capHits +
    (sc.geom === 'crystal' ? ' psi4=' + o.psi4.toFixed(3) + ' mad=' + o.madelung.toFixed(4) : '') +
    (sc.geom === 'polymer' ? ' Rg=' + o.rg.toFixed(2) + ' theta2=' + o.theta2.toFixed(3) + ' manning=' + o.manning.toFixed(3) : '') +
    (sc.geom === 'colloid' ? ' zeffMax=' + o.zeffMax.toFixed(3) + ' zeff2=' + o.zeff2.toFixed(3) : '') +
    (sc.geom === 'electrode' ? ' screen=' + o.screen.toFixed(3) + ' stored=' + o.stored.toFixed(3) + ' J=' + o.Jshort.toFixed(4) : ''));
}

console.log(fails ? fails + ' FAILURES' : 'ALL PASS');
