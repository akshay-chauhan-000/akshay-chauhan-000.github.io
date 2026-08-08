---
title: "Double-Pendulum Lattice"
date: 2026-08-09
summary: "An n×n grid of chaotic double pendulums, spring-coupled to their neighbours on a torus. One slider takes the system from independent chaos to global synchronization — watched live on a phase heatmap."
math: true
tags: ["simulation", "chaos", "synchronization", "coupled oscillators"]
---

A single double pendulum is the textbook chaotic system: two rods, two
angles, and sensitive dependence that destroys predictability within a few
swings. This simulation asks what happens when you build a whole material
out of them — an $n \times n$ lattice, one double pendulum per site, each
coupled to its four nearest neighbours by torsional springs, with periodic
boundaries so the lattice is a torus with no edges.

Individually, every site wants to be chaotic. Collectively, they have to
negotiate. The result is a competition between chaos and coupling that
produces synchronized domains, travelling fronts, and phase textures that
none of the individual pendulums knows anything about.

{{< sim "pendulum-lattice" >}}

The left panel shows the pendulums themselves, coloured by the phase of the
upper angle. The right panel is a live heatmap of the same lattice —
choose the field: either phase angle, the site energy, or how well each
site agrees with its neighbours. **Click any pendulum to kick it** (shift
kicks the other way), and drag the *Lattice size* slider to resize the
lattice mid-flight — the current pattern is resampled onto the new grid
rather than discarded.

## The model

Each site carries one double pendulum with unit masses and rod lengths,
described by two angles from the vertical, $\theta_1$ (upper) and
$\theta_2$ (lower). Writing $\Delta = \theta_1 - \theta_2$, the standard
Lagrangian equations of motion take the mass-matrix form
$M(q)\,\ddot q = F(q, \dot q)$ with

<div>
$$
M = \begin{pmatrix} 2 & \cos\Delta \\ \cos\Delta & 1 \end{pmatrix},
\qquad
F = \begin{pmatrix}
-\dot\theta_2^2 \sin\Delta - 2g\sin\theta_1 + \tau_1 \\[2pt]
\;\;\,\dot\theta_1^2 \sin\Delta - g\sin\theta_2 + \tau_2
\end{pmatrix}
$$
</div>

Neighbour interactions and damping enter as the generalized torques

<div>
$$
\tau_1 = K_1 \sum_{j \in \mathcal{N}(i)} \sin\!\big(\theta_1^{(j)} - \theta_1^{(i)}\big) - \gamma\,\dot\theta_1,
\qquad
\tau_2 = K_2 \sum_{j \in \mathcal{N}(i)} \sin\!\big(\theta_2^{(j)} - \theta_2^{(i)}\big) - \gamma\,\dot\theta_2
$$
</div>

— torsional springs that pull corresponding joints of adjacent sites into
alignment, exactly the coupling of the XY and Kuramoto models. Because the
torques are applied through the inertia matrix rather than added to the
accelerations by hand, the undamped lattice is genuinely Hamiltonian: the
sum of site energies plus the bond potential
$K\sum (1 - \cos\delta\theta)$ is a conserved quantity, which is also how
the implementation is tested.

The HUD tracks the Kuramoto order parameter for each joint,

<div>
$$
R_k = \frac{1}{n^2} \left| \sum_i e^{\,i\,\theta_k^{(i)}} \right|,
$$
</div>

which reads 1 when every pendulum swings in unison and near 0 when the
phases are incoherent.

## What to look for

- **The synchronization transition.** Set gravity to 0 (free rotors — this
  matters: with gravity on, damping alone drags every site to the hanging
  state and they end up aligned without ever talking to each other). Press
  *Randomize* with damping around 0.2 and sweep the coupling. At
  $K \approx 0$ the phases freeze wherever damping caught them and $R_1$
  stays near zero; with coupling on, domains of agreement nucleate,
  coarsen, and swallow the lattice. The *neighbour sync* heatmap shows the
  domain walls as dark seams while they last.
- **Chaos spreading.** *Reset* to the quiet lattice, then kick the centre
  and watch the *site energy* heatmap: the disturbance spreads outward as a
  front, and — because the boundaries are periodic — the fronts leaving
  opposite edges re-enter from the other side and collide.
- **Two-tier coupling.** The upper and lower joints have separate coupling
  constants. Set $K_1$ high and $K_2 = 0$: upper angles lock while the
  lower rods stay chaotic — an ordered lattice wearing a disordered fringe,
  visible by flipping the heatmap between $\theta_1$ and $\theta_2$ phase.
- **Gravity off.** With $g = 0$ the pendulums become free rotors and the
  model is essentially a two-layer XY dynamics; phase textures wind and
  anneal without a preferred direction to hang toward.
- **Size matters.** Slide $n$ up mid-run. The pattern survives the
  resampling, but the balance shifts: the same kick that dominated an
  $8 \times 8$ lattice is a local ripple in a $32 \times 32$ one.

**Copy link** captures every slider, the heatmap field, and the seed in a
shareable URL.

## Implementation notes

Plain JavaScript, no dependencies —
[`static/sims/pendulum-lattice/pendulum-lattice.js`](https://github.com/akshay-chauhan-000/akshay-chauhan-000.github.io/blob/main/static/sims/pendulum-lattice/pendulum-lattice.js)
is the entire thing, with the DOM-free engine exported on
`window.PendulumLattice`.

The full lattice — a $4n^2$-dimensional ODE system, up to 4096 variables at
$n = 32$ — is integrated with classical RK4, evaluating the complete
coupled derivative at all four stages. The trick that keeps this cheap:
each stage computes $\sin$ and $\cos$ of both angles once per site, and
every other trigonometric quantity — $\sin\Delta$, $\cos\Delta$, and all
eight coupling terms — is reconstructed from that cache via
angle-difference identities. The derivative costs exactly four trig calls
per site regardless of coupling.

Rendering avoids per-site canvas state changes: pendulums are binned into
24 hue buckets by phase and drawn as one stroked path per bucket, and the
heatmap is an $n \times n$ `ImageData` scaled up without smoothing, so each
lattice site is one crisp pixel block. Changing $n$ mid-run does
nearest-neighbour resampling of the state onto the new grid instead of
resetting it.

RK4 is not symplectic, and a chaotic Hamiltonian lattice will eventually
drift; at the default time step the drift is far below what the eye can
see, and a guard resets the state cleanly if a large time step is pushed
into a blow-up.

## References

The double pendulum as a laboratory for chaos: T. Shinbrot, C. Grebogi,
J. Wisdom, J. A. Yorke, *Chaos in a double pendulum*, Am. J. Phys. **60**,
491 (1992). Synchronization of coupled oscillators: Y. Kuramoto, *Chemical
Oscillations, Waves, and Turbulence* (Springer, 1984); S. H. Strogatz,
*From Kuramoto to Crawford*, Physica D **143**, 1 (2000); J. A. Acebrón
et al., *The Kuramoto model: a simple paradigm for synchronization
phenomena*, Rev. Mod. Phys. **77**, 137 (2005).
