# Coulomb Fluids — verification scripts

Headless checks for `static/sims/coulomb-fluids/coulomb-fluids.js`. This folder
is outside `content/` and `static/`, so Hugo does not publish it.

| Script | What it checks | Runtime |
| --- | --- | --- |
| `test_core.js` | FFT against a naive DFT. PME forces and energy against an independent brute-force Ewald sum (`ref_ewald.js`). The 2-D Madelung constant. A perfect lattice being force-free. NVE energy conservation. Thermostat and stability for every scenario preset. | ~1 min |
| `ui_smoke.js` | Drives the real `mount()` UI under a DOM stub. Covers every scenario, overlay, trace, tool and slider extreme, and flags exceptions and hung frames. | ~30 s |
| `studies.js [name]` | Reproduces the numbers quoted under "Things worth trying" on the page: `conductivity`, `crystal`, `electrode`, `polymer`, `colloid`, `structure`, or `all`. | 1–20 min |
| `harness.html` | Opens the sim outside Hugo. `#cf=1&sc=electrode&warm=4000` pre-runs 4000 steps and records JS errors on the page. | — |

## Running

These are plain Node scripts:

```sh
node tools/coulomb-fluids/test_core.js
```

Without a system Node, VS Code's bundled Electron can run them:

```sh
ELECTRON_RUN_AS_NODE=1 "$LOCALAPPDATA/Programs/Microsoft VS Code/Code.exe" tools/coulomb-fluids/test_core.js
```

Run `test_core.js` and `ui_smoke.js` after editing the engine. If you change
physics or presets, rerun the relevant study and update the numbers on the page.
