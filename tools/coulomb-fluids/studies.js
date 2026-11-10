// Reproduces the numbers quoted under "Things worth trying" on
// content/simulations/coulomb-fluids/index.md. Headless, seeded, so reruns
// give the same values on the same machine/JS engine.
//
//   node studies.js [study]      study = conductivity | crystal | electrode |
//                                        polymer | colloid | structure | all
//
// Long-time averages: the whole set takes ~15-20 minutes.
const path = require('path');
const SIM = path.join(__dirname, '../../static/sims/coulomb-fluids/coulomb-fluids.js');
const CF = require(SIM);
const which = process.argv[2] || 'all';
const want = s => which === 'all' || which === s;

function make(geom, p, lB, dt, gamma) {
  const sys = new CF.System(p.seed || 9);
  sys.lB = lB; sys.dt = dt; sys.gamma = gamma || 1;
  sys.build(geom, p);
  return sys;
}
function run(sys, tau, every, onSample) {
  const steps = Math.round(tau / sys.dt);
  for (let s = 0; s < steps; s++) {
    sys.collect = s % every === every - 1;
    sys.step();
    if (sys.collect) { sys.measure(); if (onSample) onSample(s, steps); }
  }
  sys.collect = false;
}
function secondHalfMean(sys, tau, every, keys) {
  const acc = {}; let c = 0;
  run(sys, tau, every, (s, steps) => {
    if (s < steps / 2) return;
    keys.forEach(k => { acc[k] = (acc[k] || 0) + sys.obs[k]; }); c++;
  });
  const out = {}; keys.forEach(k => { out[k] = +(acc[k] / c).toFixed(3); });
  return out;
}

if (want('conductivity')) {
  // Lambda / Lambda_NE: crowding (lB -> 0), single-species control, Wien effect.
  for (const [phi, lB, E] of [[0.12, 0.1, 0.25], [0.03, 3, 0.25], [0.03, 3, 1.0]]) {
    const sys = make('bulk', { zp: 1, zm: 1, salt: phi > 0.1 ? 300 : 150, phi, seed: 5 }, lB, 0.004);
    sys.E0 = E; run(sys, 10, 1000); sys.resetStats();
    let pair = 0, c = 0; run(sys, 40, 10, () => { pair += sys.obs.pairFrac; c++; });
    console.log('conductivity phi=' + phi + ' lB=' + lB + ' E=' + E + '  Lambda/Lambda_NE=' + sys.obs.lamRatio.toFixed(3) + '  paired=' + (pair / c).toFixed(3));
  }
  const one = make('bulk', { zp: 1, zm: 1, salt: 300, phi: 0.12, seed: 5 }, 0, 0.004);
  for (let i = 0; i < one.n; i++) one.q[i] = 1;          // lB = 0: charge couples only to E
  one.E0 = 0.5; run(one, 10, 1000); one.resetStats(); run(one, 40, 10);
  console.log('conductivity single species, lB=0, E=0.5  Lambda/Lambda_NE=' + one.obs.lamRatio.toFixed(3));
}

if (want('crystal')) {
  for (const lB of [40, 25, 15, 5]) {
    const sys = make('crystal', { salt: 200, seed: 4 }, lB, 0.0015);
    const m = secondHalfMean(sys, 40, 10, ['psi4', 'T']);
    console.log('crystal lB=' + lB + '  psi4=' + m.psi4 + '  T=' + m.T);
  }
  const c = make('crystal', { salt: 200, seed: 4 }, 1, 0.0015);
  c.forces();
  console.log('crystal perfect lattice: Ec per ion / (lB/a) = ' + ((c.eReal + c.eLong + c.eSelf) / c.n * c.latticeA).toFixed(4) + '  (exact -M/2 = ' + (-CF.MADELUNG_2D / 2).toFixed(4) + ')');
}

if (want('electrode')) {
  const sys = make('electrode', { zp: 1, zm: 1, salt: 200, phi: 0.15 }, 1.43, 0.004, 0.5);
  sys.E0 = 20 / sys.L;
  let line = 'electrode gamma=0.5 V=20: ';
  run(sys, 80, 6, s => { if ((s + 1) % Math.round(10 / sys.dt) < 6 && s > 6) line += 't' + sys.t.toFixed(0) + ' bulk/applied=' + sys.obs.screen.toFixed(2) + ' | '; });
  console.log(line);
}

if (want('polymer')) {
  for (const zp of [1, 3]) {
    const sys = make('polymer', { zp, salt: 0, phi: 0.05, nm: 60 }, 5, 0.002);
    const m = secondHalfMean(sys, 150, 10, ['rg', 'theta2', 'manning', 'T']);
    console.log('polymer lB=5 z=' + zp + '  Rg=' + m.rg + '  theta(2sigma)=' + m.theta2 + '  Manning=' + m.manning + '  T=' + m.T);
  }
}

if (want('colloid')) {
  for (const [zp, salt, lB] of [[1, 40, 3], [3, 40, 3], [3, 40, 8], [3, 60, 12]]) {
    const sys = make('colloid', { zp, salt, phi: 0.08, Z: 48 }, lB, lB > 5 ? 0.0015 : 0.003);
    const m = secondHalfMean(sys, 60, 10, ['zeff1', 'zeffMax']);
    console.log('colloid Z=-48 z=' + zp + ' salt=' + salt + ' lB=' + lB + '  Zeff(R+1.5)/|Z|=' + m.zeff1 + '  max Zeff/|Z|=' + m.zeffMax);
  }
}

if (want('structure')) {
  const sys = make('bulk', { zp: 1, zm: 1, salt: 300, phi: 0.12, seed: 5 }, 5, 0.003);
  run(sys, 18, 1000); sys.resetStats(); run(sys, 60, 20);
  const k1 = 2 * Math.PI / sys.L, kap = sys.obs.kappa;
  let line = 'S_ZZ(k)/<q^2> measured vs Debye-Huckel k/(k+kappa), lB=5:';
  for (const m of [1, 2, 3, 5, 8, 12, 20, 30]) { const k = m * k1; line += '  k=' + k.toFixed(2) + ' ' + sys.skAxis[m].toFixed(3) + '/' + (k / (k + kap)).toFixed(3); }
  console.log(line);
}
