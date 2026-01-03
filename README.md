
# MTEX MATLAB scripts — KS/NW Orientation Relationship (SDSS 2507, EBSD)

This folder contains MTEX-based MATLAB scripts to quantify **Kurdjumov–Sachs (K–S)** and **Nishiyama–Wassermann (N–W)** orientation relationship (OR) adherence at **α/γ interphase boundaries** in **LPBF SDSS 2507** across all processing conditions.

The workflow:
- Loads `.ang` EBSD data (EDAX/OIM style)
- Reconstructs grains (two-phase)
- Extracts **α–γ interphase boundary segments**
- Computes the **measured α→γ misorientation** per boundary segment
- Compares it to **all KS and NW variants** (rotation space)
- Reports **weighted fraction** of boundary length satisfying KS or NW within a user-defined tolerance
- Exports maps and histograms as publication-ready PNGs

---

## Repository structure (recommended)

