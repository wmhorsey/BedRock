# Aptik Simulation Specification

This document is reserved for simulation-specific assumptions, closures, numerical choices, diagnostics, and validation criteria.

The ontology belongs in [Axioms.md](Axioms.md). This file should only contain implementation-facing items needed to build and falsify simulations of that ontology.

## Status

This is a draft simulation-facing specification. Any rule stated here that is not already an ontology axiom should be treated as a closure assumption rather than as settled physics.

When possible, closure assumptions should remain explicit, configurable, and easy to disable for falsification runs.

## Hierarchy Handling Constraint

The simulator should not explicitly track shell generation numbers, hierarchical depth, or parent-child relationships.

It should only compute local field quantities such as:

- $A$
- $\nabla A$
- $\nabla^2 A$
- local curvature
- local compression
- local ambient level

Nested shells should arise from these quantities alone. Hierarchy is therefore an emergent geometric property, not an explicit data-model feature.

## Quantization Handling Constraint

The simulator should not introduce primitive quantized packets as fundamental entities unless a run is explicitly testing such an approximation.

Localized quanta should be treated as emergent, persistent geometric states of a continuous substrate field. Implementation choices should therefore prefer continuous field energy density, continuous transport, and geometry-derived confinement over hard-coded particle identities whenever the model permits it.

## Resistance Handling Constraint

The simulator should not introduce resistance as an independent primitive force law unless such a term is being added explicitly as numerical scaffolding or as a tested alternative hypothesis.

What appears as resistance in this ontology should arise from emergent tension: competing attraction across shared substrate, geometric incompatibility of local pulls, finite propagation, and limited replenishment from the surrounding field.

## 1. State Variables

### Physical State

- $A(\mathbf{x}, t)$: substrate scalar field.
- $A_{\text{amb}}(\mathbf{x}, t)$: local ambient level, defined as the locally flat plateau surrounding a taper.
- $T(\mathbf{x}, t) = \frac{1}{\alpha} \nabla^2 A(\mathbf{x}, t)$: local attraction density.
- $\mathbf{G}(\mathbf{x}, t) = \nabla A(\mathbf{x}, t)$: attraction field.
- $C(\mathbf{x}, t)$: local compression factor.
- $v_{\text{eff}}(\mathbf{x}, t) = \frac{v_0}{C(\mathbf{x}, t)}$: locally reduced effective propagation speed.
- $V(\mathbf{x})$: true-void geometry mask or equivalent void boundary representation.

### Representation State

The ontology does not force a single discretization. A simulation may use either:

- an Eulerian field representation,
- a Lagrangian particle or sample representation,
- or a hybrid representation.

Regardless of representation, the simulator should only store quantities needed to approximate the physical state and causal neighborhood. It should not store explicit hierarchy metadata.

## 2. Derived Quantities

The following quantities should be derived from the local field and geometry rather than introduced as independent labels.

- **Core contrast**

	$$
	\Delta A_{\text{core}} = A_{\text{amb}} - A_{\text{core}}
	$$

	used to distinguish void, depression, and spike behavior relative to the local ambient scale.

- **Shell surface**

	The shell is the locus of local maxima of $|\nabla A|$ relative to the outward normal from a core or composite interior.

- **Shell bias**

	$$
	\Delta A_{\text{shell}} = A_{\text{out}} - A_{\text{in}}
	$$

	evaluated across the shell surface.

- **Taper region**

	The taper is the monotonic falloff from the shell toward the local ambient.

- **Local curvature**

	$\kappa$ is derived from the shell isosurface or, in a discrete approximation, from local geometric reconstruction of the shell neighborhood.

- **Wake contrast**

	$$
	\Delta A_{\text{wake}} = A_{\text{shared}} - A_{\text{ambient}}
	$$

	in overlapping taper regions.

- **Radiance / release front**

	Radiance is the outward propagation of slowdown caused by local supply deficit. Operationally, it appears when inward attraction into a spike-like or collapsing region exceeds the rate at which surrounding substrate can replenish that flow.

	In that regime, a depression forms outside the spike or collapse surface, and the associated slowdown propagates outward along the local direction of flow. Simulations should therefore treat radiance as a propagated field effect or release front, not as a primitive emitted particle species unless such a particle approximation is explicitly being tested.

- **Void missing-interaction fraction**

	$m_{\text{void}}(\mathbf{x})$ is the fraction of the local causal neighborhood or solid angle occluded by a true void. This is the key geometric quantity that captures the loss of support near a non-interacting boundary.

	For simple circular or spherical void tests, the simulator may use geometry-specific approximations such as:

	$$
	m_{\text{void}}^{2D} = \frac{\arcsin(r / d)}{\pi}
	$$

	and

	$$
	m_{\text{void}}^{3D} = \frac{1 - \cos\theta}{2}, \qquad \sin\theta = \frac{r}{d},
	$$

	where $r$ is void radius and $d$ is distance from the void center.

- **Shell overlap diagnostic**

	$O_{\text{shell}}(\mathbf{x})$ measures incompatible shell occupancy in the same support region. It should be treated as a diagnostic of invalid or reforming structure, not as an explicit hierarchy indicator.

## 3. Evolution Laws

These laws describe the working simulation closure, not the ontology itself.

- **Directional rule**

	Flows and test particles move toward increasing $A$, in the direction of $\mathbf{G} = \nabla A$.

- **Causal propagation rule**

	No local update may depend on sources farther away than the causal horizon reachable within one step:

	$$
	\ell_{\text{causal}} \le v_{\text{eff}}(\mathbf{x}, t)\, \Delta t.
	$$

	In discrete implementations, this means interaction neighborhoods must be truncated by the local causal horizon.

- **Continuity-style transport**

	Away from true void interiors, the substrate should obey a continuity-like transport law. A minimal working form is

	$$
	\partial_t A + \nabla \cdot (A\, \mathbf{u}) = S_{\text{release}} - S_{\text{lock}},
	$$

	where $\mathbf{u}$ is the local flow field and the source terms represent internal redistribution rather than creation from nothing.

- **Geodesic closure**

	The local approximation to the Aptik Geodesic Principle should bias flow toward directions that best preserve dense, coherent attraction density within the local causal neighborhood.

	In minimal closures, this may be approximated using directions aligned with $\nabla A$, corrected by shell bias, void occlusion geometry, and local curvature.

- **Compression law**

	$C(\mathbf{x}, t)$ must increase monotonically with local compression indicators such as $|\nabla A|$, $T$, or a validated combination of both. The exact closure remains a simulation choice and should be benchmarked separately.

- **Conservation rule**

	Total substrate should change only through explicit boundary fluxes or declared release terms. Numerical stabilization must not silently create or destroy substrate without being labeled as a non-physical safeguard.

- **Non-physical numerical terms**

	Any drag, damping, smoothing, viscosity, relaxation, overlap clamp, or similar stabilizer must be documented as numerical scaffolding and remain independently toggleable.

## 4. Void and Boundary Handling

- **True void interior**

	A true void is a region with $A = 0$ and no substrate interaction. It is not a low-field region that still participates weakly; it is a non-interacting boundary.

- **Void boundary effect**

	A void is not a direct force source. Its effect arises from missing interaction channels on the void-facing side of the local neighborhood. The substrate near the void therefore experiences a support asymmetry.

- **Boundary-derived shell formation**

	Shell formation around a true void should arise from the geometry of missing support, not from an explicitly painted radial ring profile. A valid closure may use a shell-localized kernel multiplied by the missing-interaction fraction $m_{\text{void}}$.

- **Boundary response direction**

	The net geometric response near a true void should point away from the void, because the remaining substrate support comes from the non-void side. A sign convention that points the response inward is incorrect for this ontology.

- **Far-field boundary condition for isolated tests**

	For isolated finite-domain experiments, the preferred default is an ambient-fixed far boundary so the remote substrate continues to provide the missing support that the void-adjacent region pulls against.

- **Shell non-overlap rule**

	Stable shells should not occupy the same local support region. If shell overlap grows large, the simulator should interpret that as breach, reform, spike transition, or tension release, not as a stable stacked-shell configuration.

- **Spike transition rule**

	If a structure's core ceases to behave like a depression relative to its local exterior and instead becomes equal to or greater than the surrounding field level, it should be treated as entering spike-like behavior. In that regime, the structure should pull substrate inward and shed new depressions or low-tension releases from its surface back into the field.

## 9. Diagnostics and Falsification Tests

Falsification tests must produce quantitative, inspectable numbers, not just
plausible-looking motion. Each test states a hypothesis derived from the
ontology, a controlled setup, the observables, and explicit pass/fail criteria.

### 9.1 Cone Test (causal front)

**Hypothesis.** Axiom 10 and the causal propagation rule in section 3 assert a
finite propagation speed enforced by a per-step causal horizon
$\ell_{\text{causal}} \le v_{\text{eff}}\,\Delta t$. A localized disturbance must
therefore stay inside a causal cone whose radius grows at most linearly in time.

**Setup.** Fill the domain with a quiescent uniform substrate (no voids, ambient
tension, zero velocity). Inject one compactly-supported tension pulse at the
center. Evolve with the unmodified update loop. Run a **control** copy with
identical seed and layout but no pulse, and measure every quantity as the
perturbed-minus-control difference so static lattice noise cancels and only the
causal signal remains.

**Observables (per step).**

- **Leading edge** $r_{\text{edge}}$: outermost radius with a force difference
  above a tiny floor. This is the causal front.
- **Signal front** $r_{\text{sig}}$: outermost radius above a strong-signal
  threshold. This is the detectable amplitude front and is expected to be
  diffusion-limited (slower than the causal edge).
- **Causal bound** $r_{\text{causal}} = r_{\text{pulse}} + n\,\ell_{\text{causal}}$.
- **Isotropy** via the coefficient of variation of the per-angle leading edge.

**Predictions and pass/fail.**

- *Causality (sharp falsification).* $r_{\text{edge}}(n) \le r_{\text{causal}}(n)$
  for all $n$, within one lattice spacing. Any response beyond the cone falsifies
  the causal-horizon closure.
- *Isotropy.* Leading-edge coefficient of variation should be small (target
  $< 0.15$); a square or diamond front indicates a grid/lattice artifact.
- *Emergent speed.* The fitted cone-edge speed is an emergent invariant. The
  ratio $v_{\text{edge}}/\ell_{\text{causal}}$ near $1$ means the substrate
  transmits signals near its causal limit (wave-like); a ratio near $0$ means the
  closure is diffusion-dominated and does not carry signals at the declared
  speed. The current closure is strongly diffusive (ratio $\ll 1$); raising it
  toward $O(1)$ is an open target that requires an inertial or restoring term
  rather than a purely relaxational tension update.

**Knob-robustness protocol (this is what makes it physics).** Sweep the cosmetic
gains (`attraction_strength`, `shell_overlap_gain`, `spike_threshold`,
`merge_influence_gain`, ...) and confirm the causality verdict and the cone-edge
speed are invariant. A causal speed set by the physical `signal_speed` and
insensitive to the cosmetic knobs is a law; one that drifts with every gain is
just a tuning artifact.

**Numerical caveat.** The horizon is `min(signal_speed*dt, interaction_radius)`
floored at the substrate spacing. If `signal_speed*dt` exceeds
`interaction_radius`, the causal speed is silently set by the neighbor cutoff,
not by physics; the harness reports which term binds. Keep
`spacing < signal_speed*dt < interaction_radius` for a physically meaningful run.

**Runner.** `cargo cone` (or `cargo run --release -- --cone-test`). Per-step data
is written to CSV for external analysis.

### 9.2 Constitutive law: conservative wave vs relaxational fluid

The 9.1 cone test showed the default relaxational closure does not propagate a
signal: the field front decays in place (edge/causal ~ 0.02), so `signal_speed`
acts only as the neighbor-list reach. A real signal needs a conservative
restoring response, not first-order relaxation plus velocity damping.

The `--cone-model wave` closure supplies one, using quantities the ontology
already names. With $A$ the field and $\mathbf{G} = \nabla A$ (Axiom 2):

$$
\partial_t A = -c\,\nabla\!\cdot\mathbf{G}, \qquad \partial_t \mathbf{G} = -c\,\nabla A
\quad\Rightarrow\quad \partial_t^2 A = c^2 \nabla^2 A.
$$

This is a wave at speed $c$ that conserves $E = \tfrac12\sum(A^2 + |\mathbf{G}|^2)$
(Axiom 0: no drag, no resistance) and needs no particle inertia. The mesh-free
discretization uses symmetric pair weights so the discrete divergence is the
exact negative transpose of the gradient, making the semi-discrete system
Hamiltonian; a symplectic leapfrog then has no systematic energy drift.

**Acceptance criteria (wave).**

- The $|A|$ amplitude front radiates at $\sim c$ (so `signal_speed` is a real
  speed, not a neighbor cutoff).
- Total energy stays bounded (no damping) — the direct refutation of the
  relaxational decay.
- The wavefront is isotropic.
- Finite propagation holds: nothing outruns the operator stencil.

**Caveats and open items.**

- The discrete scheme carries a faint dispersive precursor ahead of the physical
  front — a lattice artifact bounded by the stencil reach, not true acausality.
  Tightening it needs a smaller stencil or a dispersion-controlled operator.
- The wave closure currently governs the field only. Coupling it back to the
  particle-motion engine and to true voids is future work.
- With a real wave in hand, a local compression factor $C(\mathbf{x})$ (Axiom 10)
  becomes testable: it should slow the wave near high-tension cores, an emergent
  time-dilation / lensing analog.

## Remaining Sections

5. Propagation and compression rules
6. Shell formation, breach, and reform criteria
7. Composite structure rules
8. Numerical safeguards and toggles
10. Initial conditions and parameter presets