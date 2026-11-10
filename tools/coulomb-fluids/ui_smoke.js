//   node ui_smoke.js
process.on("uncaughtException", e => { console.log("UNCAUGHT", e && e.stack); process.exit(2); });
// Drives the real mount() UI under a minimal DOM stub: every scenario, many
// frames, pointer + tool + slider + select interactions. Flags exceptions and
// any frame slower than 2 s (which would indicate a hang/infinite loop).
const path = require('path').join(__dirname, '../../static/sims/coulomb-fluids/coulomb-fluids.js');
let rafQueue = [];
function noop() {}
function makeCtx() {
  const ctx = new Proxy({}, {
    get(t, k) {
      if (k in t) return t[k];
      if (k === 'createImageData') return (w, h) => ({ data: new Uint8ClampedArray(4 * w * h), width: w, height: h });
      if (k === 'createRadialGradient' || k === 'createLinearGradient') return () => ({ addColorStop: noop });
      if (k === 'measureText') return () => ({ width: 10 });
      return noop;
    },
    set(t, k, v) { t[k] = v; return true; }
  });
  return ctx;
}
function El(tag) {
  const el = {
    tagName: (tag || 'div').toUpperCase(), children: [], dataset: {}, style: {}, hidden: false,
    textContent: '', value: '', width: 300, height: 150, _listeners: {}, _cache: {}, attrs: {},
    classList: { add: noop, remove: noop, contains: () => false },
    set innerHTML(v) { this._html = v; if (v === '') this.children = []; },
    get innerHTML() { return this._html || ''; },
    get options() { return this.children; },
    appendChild(c) { this.children.push(c); if (this.tagName === 'SELECT' && this.children.length === 1) this.value = c.value; return c; },
    addEventListener(type, fn) { (this._listeners[type] = this._listeners[type] || []).push(fn); },
    dispatch(type, ev) { (this._listeners[type] || []).forEach(fn => fn(ev)); },
    setAttribute(k, v) { this.attrs[k] = v; }, getAttribute(k) { return this.attrs[k]; },
    getContext() { return this._ctx || (this._ctx = makeCtx()); },
    getBoundingClientRect() { return { left: 0, top: 0, width: 600, height: 600 }; },
    setPointerCapture: noop, contains: () => true,
    closest(sel) { return this._closest ? this._closest(sel) : null; },
    querySelector(sel) {
      if (!this._cache[sel]) {
        let tag = 'div';
        if (/canvas|data-plot/.test(sel)) tag = 'canvas';
        if (/data-sel/.test(sel)) tag = 'select';
        const e = El(tag);
        const m = /data-(\w+)="([^"]+)"/.exec(sel);
        if (m) e.dataset[m[1]] = m[2];
        this._cache[sel] = e;
      }
      return this._cache[sel];
    },
    querySelectorAll(sel) {
      if (sel === '[data-tool]') return ['wand+', 'wand-', 'add+', 'add-', 'erase'].map(t => this.querySelector('[data-tool="' + t + '"]'));
      return [];
    }
  };
  return el;
}
global.window = { devicePixelRatio: 1 };
global.navigator = {};
global.location = { hash: '', origin: 'http://x', pathname: '/', search: '' };
global.history = { replaceState: noop };
global.requestAnimationFrame = fn => rafQueue.push(fn);
const root = El('div');
global.document = {
  readyState: 'complete',
  createElement: El,
  querySelectorAll: sel => (sel === '[data-sim="coulomb-fluids"]' ? [root] : []),
  addEventListener: noop
};
require(path);
const api = root._coulombFluids;
if (!api) { console.log('FAIL mount did not run'); process.exit(1); }
let t = 0, worst = 0, fails = 0;
function frames(k) {
  for (let i = 0; i < k; i++) {
    const q = rafQueue; rafQueue = [];
    const t0 = Date.now();
    q.forEach(fn => fn(t += 16));
    const dt = Date.now() - t0;
    worst = Math.max(worst, dt);
    if (dt > 2000) { console.log('SLOW FRAME', dt, 'ms'); fails++; }
  }
}
const canvas = root.querySelector('.cf-canvas');
const sel = root.querySelector('[data-sel="sc"]');
const ov = root.querySelector('[data-sel="overlay"]');
const trSel = root.querySelector('[data-sel="trace"]');
function click(attr, val) {
  const target = root.querySelector('[data-' + attr + '="' + val + '"]');
  target.dataset[attr] = val;
  root.dispatch('click', { target: { closest: () => target } });
}
const scenarios = ['salt', 'strong', 'crystal', 'electrode', 'polymer', 'colloid'];
try {
  for (const sc of scenarios) {
    sel.value = sc; sel.dispatch('change', {});
    frames(60);
    for (const o of ['phi', 'pairs', 'field', 'plain']) { ov.value = o; ov.dispatch('change', {}); frames(5); }
    for (const tr of trSel.options) { trSel.value = tr.value; trSel.dispatch('change', {}); frames(4); }
    // wand drag
    canvas.dispatch('pointerdown', { preventDefault: noop, clientX: 300, clientY: 300, pointerId: 1 });
    for (let k = 0; k < 10; k++) { canvas.dispatch('pointermove', { clientX: 300 + 10 * k, clientY: 300 }); frames(2); }
    canvas.dispatch('pointerup', {});
    // tools
    for (const tool of ['wand-', 'add+', 'add-', 'erase']) {
      click('tool', tool);
      canvas.dispatch('pointerdown', { preventDefault: noop, clientX: 150, clientY: 420, pointerId: 1 });
      canvas.dispatch('pointerup', {});
      frames(3);
    }
    click('tool', 'wand+');
    // pause / run / reset / link
    click('act', 'run'); frames(3); click('act', 'reset'); frames(2); click('act', 'run');
    click('act', 'link'); frames(2);
    const sys = api.getSystem(), o = sys.obs;
    console.log('OK', sc.padEnd(9), 'N=' + sys.n, 't=' + sys.t.toFixed(2), 'T=' + (o.T || 0).toFixed(2),
      'calc chars=' + root.querySelector('[data-calc]').innerHTML.length, 'worst frame ' + worst + ' ms');
    worst = 0;
  }
  // sliders: fire input on every slider at a few positions (incl. rebuild + rescale paths)
  const inputs = [];
  const host = root.querySelector('.cf-controls');
  (function walk(e) { (e.children || []).forEach(c => { if (c.tagName === 'INPUT') inputs.push(c); walk(c); }); })(host);
  for (const sc of ['salt', 'polymer', 'colloid', 'electrode']) {
    sel.value = sc; sel.dispatch('change', {});
    for (const inp of inputs) {
      for (const frac of [0, 1, 0.5]) {
        const mn = parseFloat(inp.min), mx = parseFloat(inp.max);
        inp.value = String(mn + frac * (mx - mn));
        inp.dispatch('input', {});
        frames(2);
      }
    }
  }
  // let debounced rebuilds fire
  setTimeout(() => {
    frames(30);
    const sys = api.getSystem();
    console.log('sliders exercised:', inputs.length, '| after slider sweep N=' + sys.n + ' finite=' + isFinite(sys.x[0]));
    console.log(fails ? fails + ' FAILURES' : 'UI SMOKE PASS');
  }, 400);
} catch (e) {
  console.log('EXCEPTION', e && e.stack);
  process.exit(1);
}
