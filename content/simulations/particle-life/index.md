---
title: "Particle Life"
date: 2026-08-07
summary: "A few hundred lines of physics: N particles, K species, one asymmetric interaction matrix. Break Newton's third law and the system starts building things."
math: true
tags: ["simulation", "active matter", "self-assembly"]
---

Give a few thousand point particles a short-ranged force law, sort them into a
handful of species, and let the strength of the interaction depend on *which*
species is looking at *which*. Then do one slightly illegal thing: let the
matrix of interaction strengths be **asymmetric**, so that species 1 can chase
species 2 while species 2 runs from species 1.

That single change is the whole model. There is no chemistry in it, no bonding
rule, no template, nothing that encodes "make a cell membrane." What comes out
anyway are cells with skins, crawling worms, rotating mills, and clusters that
divide when they get too large.

{{< sim "particle-life" >}}

Drag on the viewport to stir; hold shift to pull instead of push. The matrix on
the right is live — drag any cell up or down and watch the structures reorganise
within a second or two.

## The model

Each particle carries a position, a velocity, and a species label
$s_i \in \\{1,\dots,K\\}$. Every pair closer than a cutoff $r_\max$
interacts through a radial force that depends only on the reduced separation
$r = |\mathbf{r}\_{ij}| / r_\max$ and on the matrix entry
$A\_{s_i s_j}$:

<div>
$$
f(r) = \begin{cases}
\dfrac{r}{\beta} - 1, & 0 \le r < \beta \\[6pt]
A_{s_i s_j}\left(1 - \dfrac{\left|2r - 1 - \beta\right|}{1 - \beta}\right), & \beta \le r < 1 \\[6pt]
0, & r \ge 1
\end{cases}
$$
</div>

The first branch is a species-independent core repulsion: it rises linearly
from $-1$ at contact to $0$ at $r = \beta$, and it is the only thing
stopping the system from collapsing to a point. The second is a symmetric tent
peaking midway between $\beta$ and $1$, scaled by the matrix entry —
positive entries attract, negative entries repel. Beyond the cutoff, nothing.

Dynamics are overdamped, which is the right limit for anything moving through a
viscous medium at this scale. Momentum is not conserved and is not meant to be:

<div>
$$
\mathbf{v}_i \leftarrow \mathbf{v}_i\, 2^{-\Delta t / \tau} + \mathbf{F}_i\, \Delta t,
\qquad
\mathbf{x}_i \leftarrow \mathbf{x}_i + \mathbf{v}_i\, \Delta t \pmod{\mathbf{L}}
$$
</div>

with $\tau$ the velocity half-life (the *Velocity half-life* slider) and
periodic boundaries in both directions, so the domain is a torus with no walls
and no corners for structures to hide in.

### Why the asymmetry matters

If $A$ is symmetric, the force law derives from a pair potential. The system
has a well-defined energy, it relaxes downhill into that energy's minima, and
what you get is condensed-matter-flavoured: droplets, phase separation,
occasionally a lattice. Everything stops moving eventually.

If $A$ is asymmetric, $f_{ij} \neq -f_{ji}$, no potential exists, and
there is no energy to minimise. Each pair is a tiny motor injecting momentum
into the system, and structures persist because they are continuously driven,
not because they are stable. This is the defining feature of **active matter** —
the same reason a school of fish and a flock of birds have dynamics no
equilibrium statistical mechanics will reproduce.

The *Matrix symmetry* slider interpolates between the two. Set it to 1 and watch
the system go quiet.

## Things worth trying

- **Chase loops.** Set a 3-species cycle by hand: 1 attracted to 2, 2 attracted
  to 3, 3 attracted to 1, with the reverse entries negative. You get travelling
  packets that never settle.
- **Self-enclosing cells.** Make one species self-attracting and a second
  species attracted to the first but self-repelling. The second forms a shell
  around the first — a membrane, from two numbers.
- **Cutoff sweep.** Fix a matrix you like, then walk $r_\max$ from 0.02
  upward. Structures have a preferred size set by $r_\max$, and there is a
  fairly sharp point where local clusters percolate into a single network.
- **Core radius.** Small $\beta$ gives dense, liquid-like blobs; large
  $\beta$ gives open, chain-like structures, because the attractive tent
  gets squeezed into a narrow shell.
- **Turn the force off.** Set the gain to 0 and the system is a collisionless
  gas coasting to a halt under drag — a useful sanity check that nothing but the
  matrix is doing the work.

Every control is captured in the URL, matrix included. **Copy link** gives you a
permalink that reproduces the exact configuration, quantised matrix and random
seed and all.

## Implementation notes

The simulation is plain JavaScript with no dependencies and no build step —
[`static/sims/particle-life/particle-life.js`](https://github.com/akshay-chauhan-000/akshay-chauhan-000.github.io/blob/main/static/sims/particle-life/particle-life.js)
is the whole thing. Three details make it fast enough to run at a few thousand
particles inside a blog post:

**Uniform hash grid.** Naive pair enumeration is $O(N^2)$; at $N = 4000$
that is 16 million distance checks per step, far too slow for 60 fps. Instead
particles are binned into a grid whose cell size is at least $r_\max$, so
every interacting pair lies within the $3\times3$ block around a particle's
own cell. Binning uses a counting sort into preallocated typed arrays, which
costs one pass to count, one prefix sum, and one pass to place — no allocation
and no per-frame garbage. Cost per step drops to $O(N)$ at fixed density.

**Grid-order traversal.** The force loop iterates particles in the sorted grid
order rather than by index. Neighbours in space end up adjacent in memory, which
turns most of the inner loop's reads into cache hits.

**Batched drawing.** Species are stored in contiguous blocks, so rendering emits
one canvas path per species — $K$ fill calls per frame instead of $N$.

Everything lives in `Float32Array`s, the integrator is semi-implicit Euler
(positions updated only after every velocity is known), and the random number
generator is a seeded `mulberry32`, so a given seed reproduces a run exactly.
The engine object is DOM-free and exported on `window.ParticleLife`, which makes
it straightforward to drive headlessly if you want to measure something rather
than watch it.

One deliberate compromise: there is a speed cap. Push the force gain and time
step high enough and semi-implicit Euler will go unstable, so velocities are
clamped and non-finite values reset. It is a guard rail, not physics — if you
are near it, shorten the half-life or lower the gain rather than trusting what
you see.

## References

The model was popularised by Jeffrey Ventrella as *Clusters* and by Tom Mohr's
*Particle Life* implementation; the force law above follows the formulation in
Tom Mohr's write-up. The physics it gestures at is the literature on active
matter and non-reciprocal interactions — a good entry point is Fruchart,
Hanai, Littlewood and Vitelli, *Non-reciprocal phase transitions*, Nature **592**,
363 (2021), which works out what non-reciprocity does to phase behaviour far
more carefully than a canvas demo can.
