"""
readers/read_pyamares.py
Basis Set Converter — PyAMARES prior knowledge CSV reader

Reads a PyAMARES prior knowledge CSV file into a peak params dict.
This dict can then be passed to write_oxsa() to produce an OXSA
MATLAB struct, or used for inspection/display.

PyAMARES CSV format (from generate_pyamares_pk_from_basis.py):

    Index, Asc, Asp, Cho, CrCHtwo, ...       ← metabolite names
    Initial Values, , , , ...
    amplitude, 0.028, 0.040, 0.270, ...
    chemicalshift, 3.742, 2.734, 3.198, ...
    linewidth, 6.0, 6.0, 5.0, ...
    phase, 0.0, 0.0, 0.0, ...
    g, 0.0, 0.0, 0.0, ...
    Bounds, , , , ...
    amplitude, (0,, (0,, (0,, ...             ← open upper bound
    chemicalshift, (3.66,3.82), (2.65,2.81), ...
    linewidth, (1,80), (1,80), ...
    phase, (-180,180), (-180,180), ...
    g, (0,1), (0,1), ...

Peak params dict returned:
    {
        'names'         : list of str
        'amplitude'     : list of float  — initial values
        'chemicalshift' : list of float  — initial values in ppm
        'linewidth'     : list of float  — initial values in Hz
        'phase'         : list of float  — initial values in degrees
        'g'             : list of float  — 0=Lorentzian, 1=Gaussian
        'amp_bounds'    : list of tuple  — (lower, upper) or (lower, None)
        'cs_bounds'     : list of tuple  — (lower_ppm, upper_ppm)
        'lw_bounds'     : list of tuple  — (lower_hz, upper_hz)
        'phase_bounds'  : list of tuple  — (lower_deg, upper_deg)
        'g_bounds'      : list of tuple  — (0, 1) always
        'source'        : str, path to .csv file
    }

Entry point:
    read_pyamares(path) -> peak params dict
"""

import os
import csv
import numpy as np


def read_pyamares(path):
    """
    Read a PyAMARES prior knowledge CSV file.

    Parameters
    ----------
    path : str
        Full path to a PyAMARES .csv file.

    Returns
    -------
    peak params dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r', newline='') as f:
        rows = list(csv.reader(f))

    if not rows:
        raise RuntimeError(f"Empty CSV file: {os.path.basename(path)}")

    ############## Parse header row — metabolite names #############
    header = rows[0]
    if header[0].strip().lower() != 'index':
        raise RuntimeError(
            f"Expected 'Index' as first cell in {os.path.basename(path)}, "
            f"got '{header[0]}'"
        )
    names = [n.strip() for n in header[1:] if n.strip()]
    n_metabs = len(names)

    ############## Separate initial values and bounds sections #############
    init_rows   = {}
    bounds_rows = {}
    section     = None

    for row in rows[1:]:
        if not row:
            continue
        label = row[0].strip().lower()

        if label == 'initial values':
            section = 'init'
            continue
        if label == 'bounds':
            section = 'bounds'
            continue

        if section == 'init' and label in ('amplitude', 'chemicalshift',
                                            'linewidth', 'phase', 'g'):
            vals = row[1:n_metabs + 1]
            init_rows[label] = [_parse_float(v) for v in vals]

        if section == 'bounds' and label in ('amplitude', 'chemicalshift',
                                              'linewidth', 'phase', 'g'):
            vals = row[1:n_metabs + 1]
            bounds_rows[label] = [_parse_bound(v) for v in vals]

    ############## Validate required initial value rows are present #############
    required = ['amplitude', 'chemicalshift', 'linewidth', 'phase', 'g']
    missing  = [r for r in required if r not in init_rows]
    if missing:
        raise RuntimeError(
            f"Missing initial value rows in {os.path.basename(path)}: "
            f"{', '.join(missing)}"
        )

    ############## Build result dict #############
    result = {
        'names'         : names,
        'amplitude'     : init_rows['amplitude'],
        'chemicalshift' : init_rows['chemicalshift'],
        'linewidth'     : init_rows['linewidth'],
        'phase'         : init_rows['phase'],
        'g'             : init_rows['g'],
        'amp_bounds'    : bounds_rows.get('amplitude',
                            [(0.0, None)] * n_metabs),
        'cs_bounds'     : bounds_rows.get('chemicalshift',
                            [(cs - 0.08, cs + 0.08)
                             for cs in init_rows['chemicalshift']]),
        'lw_bounds'     : bounds_rows.get('linewidth',
                            [(1.0, 80.0)] * n_metabs),
        'phase_bounds'  : bounds_rows.get('phase',
                            [(-180.0, 180.0)] * n_metabs),
        'g_bounds'      : bounds_rows.get('g',
                            [(0.0, 1.0)] * n_metabs),
        'source'        : path,
    }

    # Suppress verbose output — GUI displays this instead
    return result


############## Internal helpers #############

def _parse_float(val):
    """
    Parse a float from a CSV cell.
    Returns 0.0 if empty or unparseable.
    """
    val = str(val).strip()
    if not val:
        return 0.0
    try:
        return float(val)
    except ValueError:
        return 0.0


def _parse_bound(val):
    """
    Parse a bound string into a (lower, upper) tuple.

    Handles:
        '(0,'          → (0.0, None)    open upper bound
        '(0.0, None)'  → (0.0, None)
        '(1,80)'       → (1.0, 80.0)
        '(-180,180)'   → (-180.0, 180.0)
        '(3.66,3.82)'  → (3.66, 3.82)
        ''             → (None, None)
    """
    val = str(val).strip()
    if not val:
        return (None, None)

    # Remove surrounding parentheses
    val = val.strip('()')

    parts = [p.strip() for p in val.split(',')]

    lower = _safe_float(parts[0]) if len(parts) > 0 else None
    upper = _safe_float(parts[1]) if len(parts) > 1 else None

    return (lower, upper)


def _safe_float(val):
    """Convert to float, returning None for empty or 'None' strings."""
    val = str(val).strip()
    if not val or val.lower() == 'none':
        return None
    try:
        return float(val)
    except ValueError:
        return None


############## Quick test #############

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python read_pyamares.py <path_to_csv>")
        sys.exit(1)

    result = read_pyamares(sys.argv[1])

    print(f"\nNames : {result['names']}")
    print(f"CS    : {result['chemicalshift'][:3]} ...")
    print(f"Amp   : {result['amplitude'][:3]} ...")
    print(f"LW    : {result['linewidth'][:3]} ...")
    print(f"\nCS bounds sample: {result['cs_bounds'][:3]}")
    print(f"LW bounds sample: {result['lw_bounds'][:3]}")
