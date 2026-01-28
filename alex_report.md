# Report on profiling of `Yfd3TET4JOINT`

## Overview of profiling

- Consumes x % of wall time
- Scales poorly (RESULTS)

## Suggestions

- Use OpenMP `reduction` rather than atomics or critical sections
- Watch out for false sharing
- Cache results of transcendental functions, if possible
- Data-oriented design: make your data structures vectorisable and cache-friendly

## Misc advice

- Too many variables being passed around. Define structs with meaningful names to collect variables
  together.
- Old style of C with variable declarations at start of function block. Don't do this.
