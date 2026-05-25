"""
readers/read_midas.py
Basis Set Converter — MIDAS Prior XML reader

Reads MIDAS FITT_Generic_XML prior files.

Format:
    <FITT_Generic_XML>
      <PRIOR_METABOLITE_INFORMATION>
        <param name="fitt_PriorLine00001"
               value="NAA++0++0++0++0++2.008++1.0++0.0"/>
        ...
      </PRIOR_METABOLITE_INFORMATION>
    </FITT_Generic_XML>

Each param value: abbr++dim1++dim2++dim3++line_idx++ppm++area++phase

Returns synthesized FIDs from peak lists (same as read_vespa.py).

Entry point:
    read_midas(path) -> list of core struct dicts
"""

import os
import numpy as np
import xml.etree.ElementTree as ET

DEFAULT_SW = 4000.0
DEFAULT_SF = None
DEFAULT_N  = 2048
PPM_CALIB  = 4.65


def read_midas(path):
    """
    Read a MIDAS FITT_Generic_XML prior file.

    Returns list of core struct dicts, one per metabolite.
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    try:
        tree = ET.parse(path)
        root = tree.getroot()
    except ET.ParseError as ex:
        raise RuntimeError(f"Could not parse XML: {ex}")

    if root.tag != 'FITT_Generic_XML':
        raise RuntimeError(
            f"Not a MIDAS prior XML file: root tag is '{root.tag}', "
            f"expected 'FITT_Generic_XML'"
        )

    prior_elem = root.find('PRIOR_METABOLITE_INFORMATION')
    if prior_elem is None:
        raise RuntimeError(
            "No PRIOR_METABOLITE_INFORMATION element found. "
            "This may not be a MIDAS prior file."
        )

    # Parse all param elements into per-metabolite peak lists
    met_peaks = {}   # name → [(ppm, area, phase), ...]

    for param in prior_elem.findall('param'):
        value = param.get('value', '')
        parts = value.split('++')
        if len(parts) < 8:
            continue
        abbr      = parts[0]
        # parts[1,2,3] = dims, parts[4] = line_idx
        ppm   = float(parts[5])
        area  = float(parts[6])
        phase = float(parts[7])

        if abbr not in met_peaks:
            met_peaks[abbr] = []
        met_peaks[abbr].append((ppm, area, phase))

    if not met_peaks:
        raise RuntimeError("No metabolite peaks found in MIDAS prior file.")

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
        print("Usage: python read_midas.py <midas_prior.xml>")
        sys.exit(1)
    results = read_midas(sys.argv[1])
    print(f"\nLoaded {len(results)} metabolites")
    for r in results:
        print(f"  {r['name']:15s}  peaks={len(r['peaks'])}  n={r['n']}")
