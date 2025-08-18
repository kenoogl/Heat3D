# Tolerance Comparison Report

**Generated:** 2025-08-18T12:12:25.084
**Problem:** 3D Heat Diffusion (Mode 3, Grid 240x240x31)
**Solver:** PBiCGSTAB with Gauss-Seidel smoother

## Experimental Setup

This study compares the numerical accuracy of different convergence tolerance values:
- **Tolerance values tested:** 1E-03, 1E-04, 1E-05, 1E-06, 1E-07, 1E-08
- **Accuracy metric:** Temperature differences along Z-direction lines
- **Reference solution:** Most strict tolerance (smallest ε)

## Center Line Results

### Temperature Statistics

| Tolerance | Min [K] | Max [K] | Mean [K] | L2 Norm | Data Points |
|-----------|---------|---------|----------|---------|-------------|
| 1E-08 | 300.000000 | 360.116950 | 351.657772 | 1.960E+03 | 31 |
| 1E-07 | 300.000000 | 360.116949 | 351.657771 | 1.960E+03 | 31 |
| 1E-06 | 300.000000 | 360.116955 | 351.657776 | 1.960E+03 | 31 |
| 1E-05 | 300.000000 | 360.116613 | 351.657535 | 1.960E+03 | 31 |
| 1E-04 | 300.000000 | 360.110641 | 351.653122 | 1.960E+03 | 31 |
| 1E-03 | 300.000000 | 359.934699 | 351.531376 | 1.959E+03 | 31 |

### Accuracy Comparison (vs Reference: ε=1E-08)

| Tolerance | Max Abs Diff [K] | Mean Abs Diff [K] | Max Rel Diff [%] | Mean Rel Diff [%] | RMS Diff [K] |
|-----------|------------------|-------------------|------------------|-------------------|--------------|
| 1E-08 | 0.000000 | 0.000000 | 0.0000 | 0.000000 | 0.000000 |
| 1E-07 | 0.000000 | 0.000000 | 0.0000 | 0.000000 | 0.000000 |
| 1E-06 | 0.000006 | 0.000004 | 0.0000 | 0.000001 | 0.000005 |
| 1E-05 | 0.000359 | 0.000237 | 0.0001 | 0.000066 | 0.000259 |
| 1E-04 | 0.006951 | 0.004650 | 0.0019 | 0.001303 | 0.005043 |
| 1E-03 | 0.194286 | 0.126396 | 0.0540 | 0.035383 | 0.139945 |

### Convergence Analysis

**Recommended tolerance for high precision:** 1E-05
- Achieves less than 0.001K maximum temperature difference

**Recommended tolerance for engineering analysis:** 1E-04
- Achieves less than 0.01K maximum temperature difference

## Convergence History Files

## Conclusions

### Key Findings

1. **Convergence Tolerance Impact**
   - Stricter tolerances provide more accurate solutions
   - Temperature differences decrease significantly with tighter tolerances
   - Computational cost increases with stricter tolerances

2. **Practical Recommendations**
   - For engineering analysis: Use tolerance 1.0E-4 to 1.0E-5
   - For high-precision research: Use tolerance 1.0E-6 to 1.0E-7
   - For production simulations: Balance accuracy vs computational cost

3. **Solution Quality Assessment**
   - Temperature field shows good convergence behavior
   - Relative differences decrease with stricter tolerances
   - Both center line and TSV line show consistent patterns

---
*Report generated automatically by Heat3D tolerance comparison analysis*
