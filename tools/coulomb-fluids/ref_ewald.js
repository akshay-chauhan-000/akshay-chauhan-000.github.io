const { erfc } = require(require('path').join(__dirname, '../../static/sims/coulomb-fluids/coulomb-fluids.js'));
// ---- reference Ewald (quasi-2D, 1/r), high accuracy, independent alpha
function refEwald(N, x, y, q, L, lB) {
  const a = 6.0 / L * 1.0;           // real-space sum within minimum image + 1 shell
  const A = L * L;
  let E = 0; const fx = new Float64Array(N), fy = new Float64Array(N);
  // real space over images |n| <= 1
  const rcut = 1.5 * L;
  for (let i = 0; i < N; i++) for (let j = i; j < N; j++) {
    for (let nx = -2; nx <= 2; nx++) for (let ny = -2; ny <= 2; ny++) {
      if (i === j && nx === 0 && ny === 0) continue;
      const dx = x[j] - x[i] + nx * L, dy = y[j] - y[i] + ny * L, r2 = dx * dx + dy * dy, r = Math.sqrt(r2);
      if (r > rcut) continue;
      const fac = (i === j) ? 0.5 : 1;
      const ec = erfc(a * r);
      E += fac * lB * q[i] * q[j] * ec / r;
      if (i !== j) {
        const fr = lB * q[i] * q[j] * (ec / r + 2 * a / Math.sqrt(Math.PI) * Math.exp(-a * a * r2)) / r2;
        fx[i] -= fr * dx; fy[i] -= fr * dy; fx[j] += fr * dx; fy[j] += fr * dy;
      }
    }
  }
  // k space, direct sum
  const kmax = 2 * a * 5.0, nmax = Math.ceil(kmax * L / (2 * Math.PI));
  for (let mx = -nmax; mx <= nmax; mx++) for (let my = -nmax; my <= nmax; my++) {
    if (mx === 0 && my === 0) continue;
    const kx = 2 * Math.PI * mx / L, ky = 2 * Math.PI * my / L, k = Math.hypot(kx, ky);
    if (k > kmax) continue;
    const G = 2 * Math.PI * erfc(k / (2 * a)) / k;
    let Sr = 0, Si = 0;
    for (let i = 0; i < N; i++) { const ph = kx * x[i] + ky * y[i]; Sr += q[i] * Math.cos(ph); Si += q[i] * Math.sin(ph); }
    E += lB / (2 * A) * G * (Sr * Sr + Si * Si);
    for (let i = 0; i < N; i++) {
      const ph = kx * x[i] + ky * y[i];
      // dE/dx_i = (lB/2A) G * 2 q_i kx (-sin ph Sr + cos ph Si)
      const g = lB / A * G * q[i] * (-Math.sin(ph) * Sr + Math.cos(ph) * Si);
      fx[i] -= g * kx; fy[i] -= g * ky;
    }
  }
  let s2 = 0; for (let i = 0; i < N; i++) s2 += q[i] * q[i];
  E -= lB * a / Math.sqrt(Math.PI) * s2;
  return { E, fx, fy };
}


module.exports = refEwald;
