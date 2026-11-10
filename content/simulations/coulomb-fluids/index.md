---
title: "Coulomb Fluids"
date: 2026-11-10
summary: "Hundreds of ions with the full long-range Coulomb law, solved by particle-mesh Ewald on every step: salt, an ionic crystal that melts, electrodes that charge, a polyelectrolyte that collapses. The theory is checked live as it runs."
math: true
tags: ["simulation", "electrostatics", "electrolytes", "polyelectrolytes"]
---

Coulomb's law is the one interaction in soft matter you are not allowed to
cut off. A Lennard-Jones tail is negligible beyond a few diameters; $1/r$
never is. Every ion feels every other ion, and in a periodic box it feels
all of their images too, out to infinity. Truncate the sum and you get
artefacts that look like physics. Do it properly and you get the real
thing: screening, ion pairing, crystals, double layers, counterion
condensation, the whole repertoire that makes charged matter strange.

This page does it properly, in your browser, every time step. The ions are
discs confined to a plane, and the long-range sum is a smooth particle-mesh
Ewald solver, from the same family as LAMMPS's `kspace_style pppm` and
GROMACS's PME, written from scratch in plain JavaScript. Six scenarios share the one engine, and the panel under the
viewport computes what theory says each one should do. You can watch
theory and simulation agree, and watch them fall apart.

{{< sim "coulomb-fluids" >}}

**Hold the mouse on the viewport** to wield a charged wand (*Wand +* /
*Wand −*, strength set by the *Wand charge* slider) and drag counterion
clouds around. *Drop* adds an ion of the current valence where you click,
*Erase* removes one. In the colloid scenario you can grab the colloid
itself and pull it through its own ion cloud. Hovering shows the potential
and field under the cursor (the smooth Ewald part, averaged over about
$1/\alpha \approx 2\sigma$, which is also what the potential map draws). The *Scenario* menu loads a preset,
*Overlay* switches between the potential map, Bjerrum pair bonds and field
arrows, and *Trace* picks what the strip chart follows. **Copy link**
encodes every parameter and the random seed in the URL.

## The model

Each particle is a disc of diameter $\sigma$ carrying a valence $q_i$. Pairs
interact through the bare Coulomb law plus a purely repulsive WCA core
(Lennard-Jones cut at its minimum, so ions touch but never overlap):

<div>
$$
\frac{U_{ij}(r)}{k_BT} = \ell_B\,\frac{q_i q_j}{r}
+ \Big[4\big(r^{-12} - r^{-6}\big) + 1\Big]_{\,r < 2^{1/6}\sigma},
\qquad
\ell_B = \frac{e^2}{4\pi\varepsilon_0\varepsilon_r k_BT}.
$$
</div>

The Bjerrum length $\ell_B$ is the distance at which two unit charges
interact with energy $k_BT$. It is the one coupling constant that matters,
and the main slider. In water at room temperature $\ell_B = 0.71$ nm;
taking $\sigma = 0.5$ nm puts water at $\ell_B \approx 1.43\\,\sigma$, and the
readout converts any $\ell_B$ back to the equivalent dielectric constant
$\varepsilon_r$. Raising $\ell_B$ is the same as lowering $\varepsilon_r$ or
the temperature.

The ions live in a plane, but the field lines do not: the interaction is
the three-dimensional $1/r$, as for ions bound to a membrane, an interface
or a thin film. That choice has consequences you can see, as described below.
The solvent is implicit, a viscous continuum entering through Langevin
dynamics with friction $\gamma$, so every ion has free diffusion constant
$D = k_BT/m\gamma$. The integrator is BAOAB (Leimkuhler and Matthews), which
samples configurations accurately at large steps for one force evaluation
per step.

The scenarios add structure on top of the same engine:

- **Polyelectrolyte**: a Kremer–Grest bead-spring chain (FENE bonds plus
  WCA) with charge $-1$ on every monomer, and its counterions.
- **Colloid**: a pinned disc of radius $4\sigma$ and charge $-Z$; ions see
  a WCA wall shifted to its surface.
- **Electrode**: an impenetrable wall at $x = 0 \equiv L$ and a uniform
  applied field $E_0$, so the potential drop across one period is
  $V = E_0 L$. With the periodic (tin-foil) Ewald sum, this is a cell
  held at fixed voltage.

### Ewald in a plane

The trick is to split $1/r$ into a short-ranged piece and a smooth piece,

<div>
$$
\frac{1}{r} = \frac{\operatorname{erfc}(\alpha r)}{r} + \frac{\operatorname{erf}(\alpha r)}{r}.
$$
</div>

The first piece dies within a few $1/\alpha$ and is summed directly over
neighbours with a cell list. The second is smooth everywhere, so it is
summed in Fourier space. For charges confined to a plane that is periodic
in $x$ and $y$ but not in $z$, the in-plane transform at $z=0$ is
$2\pi\\,\operatorname{erfc}(k/2\alpha)/k$, and the total energy is

<div>
$$
\frac{E}{k_BT} = \ell_B\!\!\sum_{i \lt j,\ r \lt r_c}\!\! q_iq_j\frac{\operatorname{erfc}(\alpha r_{ij})}{r_{ij}}
+ \frac{\ell_B}{2A}\sum_{\mathbf k\neq 0}\frac{2\pi\operatorname{erfc}(k/2\alpha)}{k}\,\big|S(\mathbf k)\big|^2
- \frac{\ell_B\alpha}{\sqrt\pi}\sum_i q_i^2 ,
$$
</div>

where $S(\mathbf k) = \sum_j q_j e^{i\mathbf k\cdot\mathbf r_j}$ is the
Fourier-transformed charge density.

This is the quasi-2-D Ewald sum. Particle-mesh Ewald makes the $k$-space
part cheap: spread the charges onto a mesh, take one FFT, multiply by the
influence function, and transform back. The cost drops from $O(N^2)$ to
$O(N\log N)$.

### What the calculations panel computes

Every number in the panel is measured from the running simulation, and
where a theory exists its prediction sits next to the measurement. The
theory has to be the theory for *this* geometry, and some textbook
formulas change.

**Screening is not exponential here.** Linearise Poisson–Boltzmann for a
sheet of mobile charge with $1/r$ interactions and the screened potential
of a test charge is, in Fourier space,

<div>
$$
\hat\phi(k) = \frac{2\pi\ell_B q}{k + \kappa},
\qquad
\kappa = 2\pi\ell_B\sum_i n_i q_i^2 ,
$$
</div>

with $n_i$ the areal densities. The non-analytic $k$ in the denominator
makes the potential decay as $1/r^3$ at long range, not as
$e^{-\kappa r}$: screening in a plane is algebraic. The panel reports
$1/\kappa$ as the crossover length and does not claim a Debye length. The same
theory predicts the charge–charge structure factor
$S\_{ZZ}(k)/\langle q^2\rangle = k/(k+\kappa)$, the dashed curve in the
right-hand plot, measured by direct Fourier sums along the axes and
diagonals. As $k \to 0$ it must vanish: long-wavelength charge fluctuations
are suppressed (perfect screening, the content of the Stillinger–Lovett
condition). Debye–Hückel holds at small $k$ even at strong coupling, and
fails at large $k$, where a peak grows as ions order locally.

**Bjerrum pairs.** Bjerrum's criterion counts a $+/-$ pair as bound when it
sits inside the minimum of the radial probability $P(r)$. In three dimensions
$P(r) \propto r^2 e^{\ell_B|q_+q_-|/r}$, minimal at $\ell_B|q_+q_-|/2$.
In a plane the measure is $r\\,dr$, and the minimum moves to
$r = \ell_B|q_+q_-|$. That is the cutoff used here, capped at $2\sigma$ so the
count stays meaningful for multivalent ions. Pairing is impossible when the
cutoff falls below contact, which is why the readout says so for weakly
coupled monovalent salt.

**Conductivity.** Under a field the panel measures the current density
$J = \sum q_i \dot x_i / A$ and compares it with the Nernst–Einstein value
for independent ions, $\Lambda\_{\rm NE} = \sum_i q_i^2 D / k_BT A$.

**Lattice energy.** The ionic crystal starts as a perfect square
(rock-salt) net, whose Coulomb energy per ion is $-M\ell_B/2a$ with the
two-dimensional Madelung constant $M = 1.6155426\ldots$

**Counterion condensation.** For a line charge with spacing $b$ between
unit charges, Manning's parameter is $\xi = \ell_B/b$. When $z\xi > 1$,
counterions of valence $z$ condense until the fraction $1 - 1/(z\xi)$ of
the chain's charge is neutralised. The panel measures the net small-ion
charge within $2\sigma$ of the chain.

## Things worth trying

The numbers below come from headless runs of exactly this code. They are
what you should see, give or take thermal noise. Some are long-time
averages: the polyelectrolyte needs a minute or two of wall-clock time
(about 100–150 $\tau$) to settle, so let the trace flatten before you read
the panel.

- **Turn up the coupling.** In *Salt solution*, drag $\ell_B$ from 1.4
  (water) towards 5 or 10. The potential map breaks into finer patches,
  the $+-$ peak of $g(r)$ grows, and the paired fraction climbs. Switch the
  overlay to *Potential + ion pairs* to see the bonds. *Strong coupling*
  loads a 3:1 salt at $\ell_B = 5$. More than 90% of the ions end up bound
  in pairs, and most trivalent cations collect anions into small,
  nearly neutral clusters.
- **Ohm's law, and where it breaks.** Give the salt a field. At area fraction
  0.12, even with Coulomb switched off ($\ell_B = 0.1$), $\Lambda/\Lambda\_{\rm NE}$
  sits near 0.55. Cations and anions stream through each other, and their
  collisions pass momentum to the implicit solvent. A single species
  drifting alone gives 0.99, which is how I know it is collisions and not a
  bug. Dilute to 0.03 and turn the coupling up to $\ell_B = 3$: the ratio is
  0.65 at $E = 0.25$ but 0.75 at $E = 1$, while the paired fraction drops
  from 42% to 33%. A strong field tears ion pairs apart and the
  conductivity rises: Onsager's second Wien effect.
- **Melt a crystal.** *Ionic crystal* starts as a perfect square net.
  Pause, hit *Reset*, and read the Madelung row: the PME energy of the
  untouched lattice is $-0.8075\\,\ell_B/a$ per ion against the exact
  $-M/2 = -0.8078$. Then run it and lower $\ell_B$.
  The square bond order $\psi_4$ is about 0.87 at $\ell_B = 40$, 0.63 at 25,
  0.32 at 15, and has fallen to noise by $\ell_B \approx 5$. The $S\_{ZZ}$
  Bragg peak at $(\pi/a, \pi/a)$ goes off-scale and is marked ▲.
- **Charge a capacitor.** *Electrode double layers* applies 20 $k_BT/e$
  (half a volt). Anions pile up against the left electrode and cations
  against the right, the current decays, and the bulk field falls to about a
  tenth of the applied value within ~70 $\tau$. The profile plot shows the
  potential dropping at the electrodes and going flat in the bulk. In this
  planar geometry the screening charge does not stay in thin layers: a line
  of charge in a plane has a logarithmic potential, so the induced charge
  spreads across the cell with long tails. That is the profile of a
  conducting strip, not of Gouy–Chapman.
- **Collapse a polyelectrolyte.** *Polyelectrolyte* starts at
  $\ell_B = 5$ with monovalent counterions. The chain stays extended
  ($R_g \approx 10\sigma$), and the neutralised fraction within $2\sigma$
  settles at 0.82 against Manning's 0.80. Now switch the counterion valence
  to +3 (this rebuilds the system). The condensed trivalent ions bridge
  distant parts of the chain, and over the next minute or so it folds up
  to $R_g \approx 5\sigma$. Like-charge attraction mediated by
  multivalent ions is how spermidine and cobalt hexammine condense DNA.
- **Dress a colloid.** *Charged colloid* surrounds a $-48$ disc with
  trivalent counterions. Compare $Z\_{\rm eff}$ at contact with monovalent
  counterions: roughly three quarters of the colloid's charge is
  neutralised within $1.5\sigma$ of its surface, against about 40%. In three dimensions, multivalent
  counterions can *overcharge* a colloid (charge inversion). In the runs I
  tried in this planar geometry (trivalent ions, $\ell_B$ up to 12, added
  salt), $Z\_{\rm eff}$ approaches zero from below but never crosses it:
  the salt ions prefer to cluster among themselves. If you find a setting
  where it crosses, the readout will say so.
- **Play with the wand.** In any scenario, a +30 wand gathers an anion
  cloud and pushes cations away. Watch the potential map: the cloud
  neutralises the wand within a couple of $\sigma$, which is screening seen
  directly.

## Implementation notes

Plain JavaScript, no dependencies, no build step:
[`static/sims/coulomb-fluids/coulomb-fluids.js`](https://github.com/akshay-chauhan-000/akshay-chauhan-000.github.io/blob/main/static/sims/coulomb-fluids/coulomb-fluids.js)
is the whole thing, with the DOM-free engine exported on
`window.CoulombFluids`.

**The mesh.** Charges are spread with cubic cardinal B-splines, and the mesh
goes through a radix-2 FFT. The field comes from $ik$-differentiation.
Because $E_x$ and $E_y$ are both real, they are packed into one complex
inverse transform as $E_x + iE_y$, so each force evaluation costs exactly
two FFTs. The smooth-PME Euler-spline factor
$|b(k)|^2 = 36/(4 + 2\cos 2\pi m/M)^2$ per axis corrects the interpolation.
The splitting adapts to the mesh spacing $h$: $\alpha = \pi/7h$ puts the mesh
Nyquist frequency at $\operatorname{erfc}(3.5)$, and $\alpha r_c = 3.5$ does
the same in real space. A lone charge on a mesh also interacts spuriously
with its own smeared image. Its average over where the particle sits in a
cell is computed analytically at set-up, as an alias sum over the
B-spline's Fourier images, and subtracted (cf. Ballenegger, Cerdà and Holm).

**Verified, not assumed.** The test harness compares the solver against an
independent, brute-force Ewald sum with a different $\alpha$ and direct
$k$-sums. The rms force error is $2$–$6\times10^{-4}$ of the rms force, and
the energy error is $\lesssim 2.5\times10^{-5}\\,k_BT$ per unit $q^2$. It
recovers the Madelung constant to $4\times10^{-4}$, and keeps a perfect
lattice force-free to $10^{-5}$. With the thermostat off, total energy over
4000 steps drifts by $3\times10^{-5}$ of the kinetic energy and fluctuates
by $10^{-4}$. A full step for 600 ions takes about a millisecond.

**Real space** uses a cell list built by counting sort, with a half-shell
stencil, so each pair is visited once. The cutoff is capped at $L/3$, with
the mesh refined to compensate, so small boxes stay exact. The kinetic
temperature is measured about each species' mean velocity, so drift in a
field is not mistaken for heat. Genuine heating from collisions under a
strong field does show up, and is real. A speed cap catches overlaps
created by dropping ions or compressing the box; the status line says when
it engages.

## References

Ewald summation: P. P. Ewald, *Die Berechnung optischer und
elektrostatischer Gitterpotentiale*, Ann. Phys. **369**, 253 (1921).
Particle-mesh Ewald: T. Darden, D. York, L. Pedersen, J. Chem. Phys.
**98**, 10089 (1993); U. Essmann et al., *A smooth particle mesh Ewald
method*, J. Chem. Phys. **103**, 8577 (1995); V. Ballenegger, J. J. Cerdà,
C. Holm, *Removal of spurious self-interactions in particle–mesh methods*,
Comput. Phys. Commun. **182**, 1919 (2011). The integrator: B. Leimkuhler,
C. Matthews, Appl. Math. Res. Express **2013**, 34 (2013). The models:
J. D. Weeks, D. Chandler, H. C. Andersen, J. Chem. Phys. **54**, 5237
(1971); K. Kremer, G. S. Grest, J. Chem. Phys. **92**, 5057 (1990).
The physics: F. H. Stillinger, R. Lovett, J. Chem. Phys. **49**, 1991
(1968); L. Onsager, *Deviations from Ohm's law in weak electrolytes*,
J. Chem. Phys. **2**, 599 (1934); G. S. Manning, J. Chem. Phys. **51**, 924
(1969); V. A. Bloomfield, *DNA condensation by multivalent cations*,
Biopolymers **44**, 269 (1997); W. M. Gelbart, R. F. Bruinsma,
P. A. Pincus, V. A. Parsegian, *DNA-inspired electrostatics*, Phys. Today
**53**(9), 38 (2000); A. Yu. Grosberg, T. T. Nguyen, B. I. Shklovskii,
*The physics of charge inversion in chemical and biological systems*, Rev.
Mod. Phys. **74**, 329 (2002). For the big picture, Y. Levin,
*Electrostatic correlations: from plasma to biology*, Rep. Prog. Phys.
**65**, 1577 (2002), is the review to read.
