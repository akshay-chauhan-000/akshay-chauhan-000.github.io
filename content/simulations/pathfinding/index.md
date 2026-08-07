---
title: "Pathfinding on Energy Landscapes"
date: 2026-08-08
summary: "From a maze to the Müller–Brown surface with one cost function: Dijkstra and A* hunting the lowest route between two minima, with temperature deciding what 'lowest' means."
math: true
tags: ["simulation", "algorithms", "energy landscapes"]
---

How do you get from one minimum of a landscape to another? In chemistry this
is the reaction-path problem — a system sitting in one basin of a potential
energy surface has to find a saddle point to reach the next basin, and the
route it takes sets the rate. Stripped of the physics, it is graph search:
discretise the landscape, assign each step a cost, and ask for the cheapest
walk from A to B.

This demo runs that search live, on landscapes from a plain maze up to the
Müller–Brown surface — the standard test case of the reaction-path
literature. The unifying idea: **a maze is just an energy landscape whose
barriers are infinite.** Same problem, same algorithms, different terrain.

{{< sim "pathfinding" >}}

Press **Run** to watch the frontier flood outward — that animation *is* the
algorithm, one expanded cell at a time. Drag the green and red markers to
move the endpoints (the path re-solves instantly as you drag). On the maze
you can paint walls with the mouse (shift erases); on the smooth landscapes
the same gesture raises hills (shift carves valleys).

## One cost function, two limits

The landscape is a grid of cells with energies $E$ normalised to $[0, 1]$.
Moving into a cell costs its geometric step length times an Arrhenius-style
penalty:

<div>
$$
c(a \to b) \;=\; \ell_{ab}\, e^{E_b / T},
\qquad \ell_{ab} \in \{1, \sqrt{2}\}
$$
</div>

where $T$ is the temperature slider. The two limits are the whole point:

- **$T \to \infty$:** the exponential flattens to 1, cost reduces to
  distance, and the optimal route is the geometrically shortest path — it
  marches straight over any barrier in the way.
- **$T \to 0$:** the exponential blows up on high ground, the maximum energy
  along the route dominates its total cost, and the optimal path becomes the
  one that crosses the landscape at its **lowest saddle**, hugging valleys
  the rest of the way. That is a discrete cousin of the minimum energy path
  that transition-state theory cares about.

The **barrier** readout in the corner shows the highest energy the found
path visits. Watching it drop as you lower $T$ — while the path visibly
lengthens — is the entire tradeoff of activated dynamics in one number.

A maze is the $E \to \infty$ caricature: walls are simply impassable, every
open cell costs the same, and the temperature slider goes quiet.

## The algorithms

All four searches grow a frontier outward from the start; they differ only
in which frontier cell they expand next.

**Dijkstra** always expands the cell with the smallest accumulated cost
$g$. It is provably optimal and completely agnostic about where the goal
is — the flood fills isotropically until it happens to touch the target.

**A\*** expands the cell minimising $f = g + w\,h$, where $h$ is the octile
distance to the goal (the exact cost of the best conceivable route on an
empty grid, so it never overestimates) and $w$ is the heuristic-weight
slider. At $w = 0$ you recover Dijkstra; at $w = 1$ the heuristic is
admissible and the result is still provably optimal, found with far fewer
expansions — you can see the flood lean toward the goal. Past $w = 1$ the
search gets greedier and faster but the optimality guarantee is gone; the
status line says so.

**Greedy best-first** expands whatever looks closest to the goal, ignoring
accumulated cost entirely. Fast, single-minded, and cheerfully wrong on
landscapes that punish it.

**BFS** expands in order of hop count, ignoring weights altogether. On the
maze it is exact (all steps cost the same); on a smooth landscape it
blithely walks over the ridge — a useful reminder that the fewest-steps
path and the cheapest path are different questions.

## Things worth trying

- **The temperature switch.** Load *Double well*. The two basins are
  separated by a ridge pierced by two saddles of different heights. At high
  $T$ the path shoots nearly straight across; drop $T$ below ~0.15 and it
  detours through the lower saddle. One slider, two mechanisms — watch the
  barrier readout snap down.
- **Müller–Brown at low temperature.** The canonical result: the good route
  from the upper-left minimum to the lower-right one is not remotely
  straight — it descends through the intermediate basin along the curved
  valley floor. Compare against BFS, which cuts straight across the high
  ground.
- **Dijkstra vs A\* on the maze.** Same optimal path, wildly different
  expansion counts. Then push $w$ past 1 and find a maze where weighted A*
  returns a genuinely longer path.
- **Seal the goal.** Paint a wall all the way around the red marker. The
  frontier exhausts and the status reports no path — exhaustion is a proof,
  not a timeout.
- **Sculpt your own.** Start from *Blank canvas*, paint a mountain range
  with a gap in it, and watch the path thread the gap at low $T$ and stop
  caring at high $T$.

**Copy link** captures landscape, seed, endpoints, and all parameters in
the URL, so a configuration can be shared and reproduced exactly
(hand-painted edits are the one thing that stays local).

## Implementation notes

Plain JavaScript, no dependencies, no build step —
[`static/sims/pathfinding/pathfinding.js`](https://github.com/akshay-chauhan-000/akshay-chauhan-000.github.io/blob/main/static/sims/pathfinding/pathfinding.js)
is the whole thing, and the engine (field, terrains, search) is DOM-free
and exported on `window.PathfindingSim` for headless use.

The open set is a **binary min-heap with lazy deletion**: decreasing a
cell's cost just pushes a duplicate entry, and stale entries are discarded
when popped. With every per-node quantity (cost, parent, state, Arrhenius
multiplier) in preallocated typed arrays, a full solve on a $201 \times
101$ grid is sub-millisecond — cheap enough to re-run the entire search on
every mouse move while you drag a marker, which is why the path feels like
it is being dragged rather than recomputed.

The terrain is baked once into an offscreen canvas at grid resolution and
scaled up (smoothly for landscapes, crisply for the maze); the animated
frontier is a second RGBA buffer composited on top. Mazes come from a
recursive backtracker on the odd sublattice, which guarantees every
corridor cell is reachable — so "no path" can only ever mean you painted
the blockage yourself. Diagonal moves are allowed but corner-cutting is
not: a diagonal step requires both adjacent orthogonal cells open.

## Where this goes

This is deliberately the low-dimensional warm-up. The serious version of
"get from one minimum to another" on a molecular potential energy surface
is the domain of the nudged elastic band and string methods, which relax an
entire chain of states toward the minimum energy path in thousands of
dimensions — grid search stops scaling long before that. The plan is to
walk there in stages: this grid world, then chain-of-states methods on
these same 2-D surfaces (where you can watch the band converge), then real
molecular landscapes. A different follow-up lives closer to the project's
original phrasing: with *many* minima on one landscape, asking for the
cheapest tour that visits all of them is a genuine travelling-salesman
problem, with pairwise costs given by the searches on this page.

## References

Dijkstra's algorithm: E. W. Dijkstra, *A note on two problems in connexion
with graphs*, Numer. Math. **1**, 269 (1959). A\* and admissibility:
P. E. Hart, N. J. Nilsson, B. Raphael, IEEE Trans. SSC **4**, 100 (1968).
The test surface: K. Müller, L. D. Brown, *Location of saddle points and
minimum energy paths…*, Theor. Chim. Acta **53**, 75 (1979). For the
methods this page is warming up to: H. Jónsson, G. Mills, K. W. Jacobsen,
*Nudged elastic band method* (in *Classical and Quantum Dynamics in
Condensed Phase Simulations*, 1998); W. E, W. Ren, E. Vanden-Eijnden,
*String method for the study of rare events*, Phys. Rev. B **66**, 052301
(2002).
