
Draft AGENTS.md for MATLAB coding guidelines
Dec 2
·
MRIsimulator
·
dev/draft-agents.md-for-matlab-coding-guidelines
·
+45
-0

Archive

Share
View PR


Help me draft an AGENTS.md file that explains how to code well in MATLAB for this project
Copy
Summary

Added repository-wide AGENTS.md outlining MATLAB coding, documentation, validation, numerical, testing, dependency, performance, and version control guidelines for MRIsimulator contributors.

Testing

⚠️ Not run (not requested)


File (1)

AGENTS.md
New



Diff
Logs

AGENTS.md
New
+45
-0

# MATLAB Coding Guidelines for MRIsimulator

## Scope
These guidelines apply to all MATLAB (.m) files in the repository.

## Style and Structure
- Favor clear, descriptive function and variable names (e.g., `signalAmplitude`, `rfPulseDuration`).
- Keep functions short and single-purpose; break complex logic into helper functions.
- Place one top-level function per file, with any local helper functions at the end of the same file.
- Use consistent indentation of four spaces; avoid tabs.
- Prefer logical sectioning with `%%` to delineate major steps for scripts.
- Order file contents as: documentation block, input validation, core logic, and output assembly.

## Documentation and Comments
- Include a MATLAB help comment block at the top of each function describing purpose, inputs, outputs, and usage examples when helpful.
- Document units for all physical quantities and clearly note coordinate systems or reference frames.
- Use inline comments sparingly to clarify non-obvious decisions; avoid restating code.

## Error Handling and Validation
- Validate inputs with `arguments` blocks when using recent MATLAB versions; otherwise use `narginchk`/`validateattributes`.
- Provide actionable error messages that describe both the issue and the expected input format or range.
- Avoid suppressing warnings unless absolutely necessary; when suppressed, add a comment explaining why.

## Numerical Practices
- Prefer vectorized operations over explicit loops when it improves readability and performance.
- Use preallocation for arrays in performance-sensitive sections.
- Be explicit about units and scaling factors; avoid hidden unit conversions.
- When using randomness, expose seeds or random streams for reproducibility.

## Testing and Examples
- Add simple, runnable examples in documentation blocks when feasible.
- Where behavior is critical, create minimal test scripts or functions under `src/tests` (or a similar folder) and keep them deterministic.

## File and Dependency Management
- Keep MATLAB paths minimal; avoid modifying the global path from within library code.
- Use relative paths based on `mfilename('fullpath')` to locate resources when needed.
- Do not add external dependencies without discussion; favor built-in MATLAB functionality.

## Performance and Memory
- Profile code (`profile`, `timeit`) for hotspots before optimizing.
- Use `parfor` or vectorized operations when beneficial, but guard parallel code with checks that fall back gracefully.

## Version Control Practices
- Keep changes focused and well-scoped; avoid mixing formatting-only edits with logic changes.
- Update this AGENTS.md if repository-wide MATLAB standards evolve.
