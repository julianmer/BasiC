"""
readers/read_gava.py
Basis Set Converter — GAVA Text reader

Reads GAVA text prior files.

Format:
    ;comment lines (start with ';')
    abbr \\t dim1 \\t dim2 \\t dim3 \\t line_idx \\t ppm \\t area \\t phase

Returns synthesized FIDs from peak lists.

Entry point:
    read_gava(path) -> list of core struct dicts
"""

import os
import numpy as np

DEFAULT_SW = 4000.0
DEFAULT_SF = None
DEFAULT_N  = 2048
PPM_CALIB  = 4.65


def read_gava(path):
    """
    Read a GAVA text prior file.

    Returns list of core struct dicts, one per metabolite.
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    met_peaks = {}   # name → [(ppm, area, phase), ...]

    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue   # skip comments and blank lines

            parts = line.split('\t')
            if len(parts) < 8:
                continue

            abbr  = parts[0].strip()
            # parts[1,2,3] = dims, parts[4] = line_idx
            ppm   = float(parts[5])
            area  = float(parts[6])
            phase = float(parts[7])

            if abbr not in met_peaks:
                met_peaks[abbr] = []
            met_peaks[abbr].append((ppm, area, phase))

    if not met_peaks:
        raise RuntimeError(
            "No metabolite data found in GAVA file. "
            "Expected tab-separated lines: abbr\\tdim1\\tdim2\\tdim3\\tidx\\tppm\\tarea\\tphase"
        )

    # Build core structs
    results = []
    for name, peaks in met_peaks.items():
        fid = _synthesize_fid(peaks)
        results.append({
            'fid'   : fid,
            'sw'    : DEFAULT_SW,
            'sf'    : DEFAULT_SF,
            'n'     : DEFAULT_N,
            'name'  : name,
            'peaks' : peaks,
            'source': path,
        })
        print(f"  Loaded: {name:15s}  {len(peaks)} peaks")

    return results


def _synthesize_fid(peaks, sw=DEFAULT_SW, sf=123.26, n=DEFAULT_N,
                    lw_hz=2.0, ppm_calib=PPM_CALIB):
    """Synthesize a Lorentzian FID from peak list."""
    dt  = 1.0 / sw
    t   = np.arange(n) * dt
    T2  = 1.0 / (np.pi * lw_hz)
    fid = np.zeros(n, dtype=complex)
    for ppm, area, phase in peaks:
        freq_hz = (ppm - ppm_calib) * sf
        fid += area * np.exp(1j * (2*np.pi*freq_hz*t + np.radians(phase))) \
                    * np.exp(-t / T2)
    return fid


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python read_gava.py <gava_prior.txt>")
        sys.exit(1)
    results = read_gava(sys.argv[1])
    print(f"\nLoaded {len(results)} metabolites")
    for r in results:
        print(f"  {r['name']:15s}  peaks={len(r['peaks'])}  n={r['n']}")
