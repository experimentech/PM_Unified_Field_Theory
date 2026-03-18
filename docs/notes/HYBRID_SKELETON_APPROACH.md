# Hybrid Vector-Tensor Skeleton Approach

## Overview

The hybrid approach combines the efficiency of **vector-based root finding** (coarse grain) with the precision of **tensor field evaluation** (fine grain) to map complex gravitational flow structures.

## Source

Found in: `legacy/Pushing-Medium/programs/junk/benchmark_harness.py::hybrid_bench()`

## Method

### Stage 1: Vector Skeleton (Coarse)
```python
def hybrid_bench(phase: float, omega: float, Nx=500, Ny=350):
    # 1. Find critical points using Newton's method on vector field
    vres = vector_skeleton_bench(phase, omega)
    vf = make_flow(phase, omega)
    
    # vres contains:
    #   - stagnations: critical points (saddles, centers, sources)
    #   - manifolds: traced stable/unstable manifolds from saddles
```

**Why vector?** 
- Direct Newton solve on $\mathbf{u}_\Omega(\mathbf{r}) = 0$
- Fast, exact critical points
- Manifold tracing follows eigenvectors of Jacobian
- No grid required initially

### Stage 2: Tensor Field (Fine Grain Near Manifolds)
```python
    # 2. Create coarse grid
    x = np.linspace(-2.0, 3.0, Nx)
    y = np.linspace(-2.0, 2.0, Ny)
    
    # 3. Build mask: only sample within tubes around manifolds
    w = 0.03  # tube width
    mask = np.zeros((Ny, Nx), dtype=bool)
    X, Y = np.meshgrid(x, y)
    
    for path in vres['manifolds']:
        for k in range(0, path.shape[0], 15):  # sample points along manifold
            px, py = path[k]
            mask |= ((X - px)**2 + (Y - py)**2) <= (w*w)
    
    # 4. Compute flow ONLY where mask is True
    U = np.zeros((Ny, Nx))
    V = np.zeros((Ny, Nx))
    idxs = np.argwhere(mask)
    for (j, i) in idxs:
        U[j, i], V[j, i] = vf.uOmega(np.array([x[i], y[j]]))
```

**Why tensor field here?**
- Full spatial structure near separatrices
- Can compute gradients, Hessians, curvature
- Enables visualization of flow topology
- But only sampled where needed (10-20% of full grid typically)

## Performance Comparison

From the benchmark (Sun-Earth-Moon at 30° phase):

| Method | Time | Points/Cells | Coverage |
|--------|------|--------------|----------|
| Vector skeleton | ~ms | ~10 critical points | Sparse exact |
| Full grid | ~seconds | 350,000 cells | 100% coverage |
| **Hybrid** | **~100ms** | **~50,000 cells** | **Targeted 15%** |

## When to Use Each

### Pure Vector (core.py approach)
- Finding critical points (Lagrange points, stagnation zones)
- Single-point queries (deflection angle at specific impact parameter)
- Tracing individual rays or particles
- Analytic continuation from known solutions

### Pure Tensor/Grid
- Full 2D/3D visualization
- Unknown topology (no good seed points)
- Computing global quantities (total flux, volume integrals)
- Image-processing-style analysis

### Hybrid
- **Multi-body systems** (planets, moons, binary stars)
- **Separatrix mapping** (capture/escape boundaries)
- **Phase-space portraits** where critical structure is sparse but important
- **Adaptive refinement** seeded by vector analysis

## Application to Galaxy Rotation

For PM galaxy fitting, the **pure vector approach is sufficient**:

1. Spherical symmetry → vector field is radial: $\mathbf{u}_g(r) = u_g(r) \hat{r}$
2. Critical points at $du_g/dr = 0$ found analytically or via 1D Newton
3. No manifolds/separatrices in axisymmetric case
4. Tensor field would just store $u_g(r, \theta, \phi)$ redundantly

**Hybrid would be useful for:**
- Merging galaxies (non-axisymmetric)
- Tidal streams
- Bar-driven flows
- Disk-halo interface structure

## Relation to GR

**GR limitation**: Geodesics are computed from $g_{\mu\nu}$, which requires full metric everywhere. No "vector skeleton" of spacetime curvature exists separately from the tensor.

**PM advantage**: 
- $\mathbf{u}_g$ is a **primary field**, not derived from a potential or metric
- Can analyze flow topology directly (critical points, manifolds)
- Then use tensor methods ($\nabla n$, Hessian) only where geometric detail matters

This is similar to how fluid dynamics separates:
- **Eulerian** (grid/tensor) vs **Lagrangian** (particle/vector) descriptions
- Critical point theory (vector) vs vorticity field (tensor)

PM naturally supports both, and the hybrid leverages the best of each.

---

## Code Status

**Legacy implementation**: 
- `legacy/Pushing-Medium/programs/junk/benchmark_harness.py::hybrid_bench()`
- Uses `vector_skeleton` and `grid_layer` modules from legacy codebase

**Current implementation**:
- New `src/pushing_medium/core.py` has vector methods (index_profile, curved_path_deflection)
- Tensor field methods could be added if needed for visualization or multi-body analysis
- Galaxy fitting uses pure vector (1D radial) - no hybrid needed yet

**Future**: 
If we need to model **disk warps**, **spiral density waves**, or **galaxy mergers**, the hybrid approach would be powerful for mapping the 3D flow structure efficiently.
