# Analytical Shapes Design and Usage

## Purpose
This document summarizes how `AnalyticalShape3D`-based shapes are defined, transformed, combined, and intensity-scaled in this repository. It also documents memory-management patterns for time-varying simulations.

## Class Hierarchy
- `AnalyticalShape3D` (abstract base)
- Basic concrete shapes:
- `AnalyticalBox3D`
- `AnalyticalCylinder3D`
- `AnalyticalEllipsoid3D`
- `AnalyticalEllipticalCylinder3D`
- `AnalyticalSphere3D` (wrapper around ellipsoid with equal axes)
- Composite/container shapes:
- `CompositeAnalyticalShape3D` (add/subtract components)
- `MultipleMaterialPhantom` (sum multiple materials)
- Domain phantoms built on top of these:
- `BreastPhantom`, `EnhancingVessel`, `BeatingHeart`, `BreathingLung`

## Core Shape Model
Each analytical shape has two core inputs:
- `shapeParameters`: geometry plus optional pose
- `shapeIntensity`: intensity scaling

Base constructor signature:
- `obj = SomeShape(intensity, shapeParameters)`

### Required pose schema
All shapes use a common pose struct (defaults to zero transform if omitted):

```matlab
params.pose.center.x_mm
params.pose.center.y_mm
params.pose.center.z_mm
params.pose.roll_deg
params.pose.pitch_deg
params.pose.yaw_deg
```

Units:
- Position in `mm`
- Rotation in `deg`
- k-space coordinates expected in `cycles/mm`

## Coordinate Frames and Transform Order
The code distinguishes:
- WORLD frame: simulation coordinates (`x,y,z` or `kx,ky,kz` passed by caller)
- BODY frame: shape-local coordinates where formulas are defined

### Image-domain path (`estimateImage`)
`AnalyticalShape3D` does:
1. subtract pose center (`WORLD -> centered`)
2. rotate world-to-body using roll/pitch/yaw
3. evaluate `percentInsideBody(...)`
4. multiply by evaluated intensity

### k-space path (`kspace`)
`AnalyticalShape3D` does:
1. rotate WORLD k-vectors into BODY frame
2. evaluate BODY analytic FT via `kspaceBaseShape(...)`
3. apply translation phase `exp(-i 2 pi (k dot center))`
4. multiply by evaluated intensity

So intensity is always a post-geometry multiplicative scale in the base class.

## Shape Parameterization (All Analytical Shapes)

### `AnalyticalBox3D`
- Required fields: `Lx_mm`, `Ly_mm`, `Lz_mm`
- Occupancy: axis-aligned bounds in BODY frame
- k-space: separable sinc product with volume scale

### `AnalyticalCylinder3D`
- Required fields: `radius_mm`, `length_mm`
- Occupancy: radial disk (`x^2+y^2 <= r^2`) and axial bound (`|z| <= L/2`)
- k-space: Bessel radial term (`J1`) times axial sinc

### `AnalyticalEllipsoid3D`
- Required fields: `a_mm`, `b_mm`, `c_mm` (semi-axes)
- Occupancy: `(x/a)^2 + (y/b)^2 + (z/c)^2 <= 1`
- k-space: closed-form ellipsoid FT using spherical-Bessel-style expression

### `AnalyticalEllipticalCylinder3D`
- Required fields: `a_mm`, `b_mm`, `length_mm`
- Occupancy: elliptical cross-section plus axial bound
- k-space: elliptical Bessel radial term times axial sinc

### `AnalyticalSphere3D`
- Accepts either:
- `radius_mm`
- or `a_mm,b_mm,c_mm` (must all match)
- Internally normalizes to ellipsoid axes `a=b=c=radius`

## Intensity Scaling Behavior

### Base behavior (`AnalyticalShape3D`)
Intensity may be:
- scalar numeric
- numeric array matching evaluation array size
- function handle returning numeric scalar/array

Applied as:
- image domain: `estimateImage = percentInsideShape .* intensity`
- k-space: `kspace = kspaceWorldPlacedShape .* intensity`

Validation requires numeric, real, finite, non-empty, and either scalar or size-matched.

### In composites and phantoms
`CompositeAnalyticalShape3D`:
- `kspace` path uses composite `shapeIntensity` (through base-class `kspace`)
- component k-space contributions are geometry-only (`kspaceWorldPlacedShape`) in `kspaceBaseShape`
- component own intensities are intentionally ignored in k-space composition

`MultipleMaterialPhantom`:
- sums each contained shape via `shape.kspace(...)` (so each child keeps its own intensity)
- parent phantom intensity is effectively fixed/unused in overridden `kspace`

## Nesting: Pose and Center Composition in `MultipleMaterialPhantom`
You can nest shapes with both child and parent poses.

Conceptually:
1. Parent phantom maps WORLD k-space to parent BODY via parent rotation
2. Each child shape then applies its own pose/intensity in that parent BODY coordinate system
3. Parent translation phase is applied after summing children

This enables hierarchical placement:
- local shape center/orientation inside a material group
- group-level center/orientation for the whole cluster

This is how complex phantoms (for example `BreastPhantom`) place organs inside a thorax and then place thorax-level structures.

## Parameter Source Patterns and Memory Management
The framework supports three storage patterns for parameter/intensity sources.

### 1) Scalars (lowest memory)
Use scalar geometry/intensity when static in time and space.

Example:
```matlab
params = struct('a_mm', 10, 'b_mm', 8, 'c_mm', 6, 'pose', pose);
shape = AnalyticalEllipsoid3D(1.0, params);
```

### 2) Vectors/arrays (precomputed trajectories)
Use vectors/arrays when values vary over time or samples and you can afford storage.

Example:
```matlab
params = struct('radius_mm', radiusVec_mm, 'length_mm', lenVec_mm, 'pose', poseVec);
shape = AnalyticalCylinder3D(1.0, params);
```

All non-scalar numeric parameter fields must share the same size.

### 3) Function handles (on-demand computation)
Use function handles to avoid storing large arrays and/or to compute lazily.

Supported places:
- full `shapeParameters` can be a function handle returning a struct
- individual parameter fields can also be function handles
- intensity can be a function handle

Pattern:
```matlab
paramsFcn = @() struct('a_mm', aNow, 'b_mm', bNow, 'c_mm', cNow, 'pose', poseNow);
shape = AnalyticalEllipsoid3D(@() intensityNow, paramsFcn);
```

In this repo, this pattern is used in time-varying classes (for example `EnhancingVessel`) to keep memory low and defer recomputation.

### Chunking for long waveforms
`cardiac_ellipsoid_waveform` explicitly processes in chunks (15,625 samples) to keep temporary memory around ~0.001 GB instead of allocating very large matrices at once.

## Usage Examples

### Basic shape with pose
```matlab
pose = struct('center', struct('x_mm', 0, 'y_mm', 0, 'z_mm', 0), ...
              'roll_deg', 0, 'pitch_deg', 10, 'yaw_deg', 20);
params = struct('Lx_mm', 60, 'Ly_mm', 40, 'Lz_mm', 20, 'pose', pose);
box = AnalyticalBox3D(1.0, params);
```

### Nested phantom
```matlab
child1 = AnalyticalEllipsoid3D(0.8, childParams1);
child2 = AnalyticalCylinder3D(1.2, childParams2);
parentPose = struct('pose', parentPoseStruct);
phantom = MultipleMaterialPhantom([child1, child2], parentPose);
```

## Important Implementation Details and Caveats
- Classes are MATLAB handle classes. Reusing the same component handle in multiple composites means later edits affect all owners.
- `shapeParameters` defaults are filled by `ensurePoseFields(...)`; missing pose fields become zero.
- Size checks are strict: non-scalar parameter arrays must be mutually size-compatible.
- `shapeIntensity` function handles are currently called with no arguments; if time is needed, pass it through a closure (for example `@() f(t_s_current)`).
- `CompositeAnalyticalShape3D.estimateImage(...)` currently sums/subtracts component `estimateImage(...)` directly and does not apply composite-level pose/intensity the same way as its k-space path.
- `MultipleMaterialPhantom.percentInsideShape(...)` unions child occupancies directly; this bypasses the parent pose transform path used by the base `AnalyticalShape3D.percentInsideShape(...)` implementation.

## Recommended Patterns
- Prefer scalar parameters/intensity for static phantoms.
- Prefer function handles for long time series or large sampled trajectories.
- If trajectories are needed repeatedly and fit memory comfortably, precompute vectors once and reuse.
- For nested phantoms, keep child poses local and use parent pose only for group transforms.
- For reproducibility/performance, use chunked updates (`kspaceAtTime` style) for long runs.
