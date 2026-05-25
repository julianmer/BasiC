"""
readers/read_vespa.py
Basis Set Converter — VeSPA reader

Reads VeSPA VIFF XML files into a list of core struct dicts.

Handles two VeSPA XML file types:

1. VeSPA Analysis Prior XML (written by write_vespa.py or exported from VeSPA):
   <vespa_export>
     <experiment>
       <metabolite>
         <name>NAA</name>
         <peak_ppm_location>2.008</peak_ppm_location>
         <peak_area_scale>1.0</peak_area_scale>
         <peak_search_width>0.1</peak_search_width>
       </metabolite>
       ...
     </experiment>
   </vespa_export>

2. VeSPA Simulation Export (from Brian Soher style):
   <vespa_export>
     <experiment>
       <simulation>
         <metabolite_id>...</metabolite_id>
         <ppms>...</ppms>     ← base64 encoded float64
         <areas>...</areas>   ← base64 encoded float64
         <phases>...</phases> ← base64 encoded float64
       </simulation>
     </experiment>
   </vespa_export>

Note: VeSPA files contain peak lists, NOT FIDs.
The 'fid' field in the returned core struct is synthesized
from the peak list using Lorentzian line shapes.

Core struct returned per metabolite:
    {
        'fid'    : np.ndarray, complex, synthesized from peaks
        'sw'     : float (Hz) — default 4000 if not in file
        'sf'     : float (MHz) — None if not in file
        'n'      : int — default 2048
        'name'   : str
        'peaks'  : list of (ppm, area, phase) tuples
        'source' : str
    }

Entry point:
    read_vespa(path) → list of core struct dicts
    read_vespa_folder(folder) → list of core struct dicts
"""

import os
import numpy as np
import xml.etree.ElementTree as ET


def _log(*args): pass  # GUI handles display


# Default spectral parameters when not stored in the file
DEFAULT_SW = 4000.0
DEFAULT_SF = None
DEFAULT_N  = 2048
DEFAULT_PPM_CALIB = 4.65
PROTON_GAMMA = 42.577   # MHz/T


def read_vespa(path):
    """
    Read a VeSPA VIFF XML file.

    Parameters
    ----------
    path : str
        Full path to a VeSPA .xml file.

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    try:
        tree = ET.parse(path)
        root = tree.getroot()
    except ET.ParseError as ex:
        raise RuntimeError(f"Could not parse XML: {ex}")

    # Detect file type
    tag = root.tag.lower()
    if 'vespa' not in tag and root.tag != 'vespa_export':
        raise RuntimeError(
            f"Not a VeSPA VIFF XML file: root tag is '{root.tag}'"
        )

    # Check for raw_fidsum (data file — not a prior)
    if root.find('.//raw_fidsum') is not None:
        raise RuntimeError(
            f"{os.path.basename(path)} is a VeSPA Analysis data file, "
            f"not a prior information file. "
            f"Please export the prior from VeSPA Simulation."
        )

    # Check for pulse_sequence without prior (simulation recipe only)
    if root.find('.//pulse_sequence') is not None and        root.find('.//prior') is None and        root.find('.//prior_metabolite') is None:
        raise RuntimeError(
            f"{os.path.basename(path)} is a VeSPA Simulation pulse sequence "
            f"export, not a prior information file. "
            f"Run the simulation in VeSPA first, then use "
            f"Simulation → Third Party Export → Vespa-Analysis."
        )

    ############## Try reading as our written prior format #############
    results = _read_prior_format(root, path)
    if results:
        return results

    ############## Try reading VeSPA Simulation experiment results #############
    results = _read_simulation_format(root, path)
    if results:
        return results

    raise RuntimeError(
        f"Could not extract metabolite prior information from "
        f"{os.path.basename(path)}. "
        f"Expected either a VeSPA Analysis prior file or a "
        f"VeSPA Simulation experiment export."
    )


def _read_prior_format(root, path):
    """
    Read VeSPA Analysis prior format (confirmed from marss_to_vespa_basis.xml):
    <prior><prior_metabolite><name>...<line><ppm>...<area>...<phase>...
    """
    results = []

    for prior in root.findall('.//prior'):
        # Get TE if present
        te_elem = prior.find('seqte')
        te = float(te_elem.text.strip()) if te_elem is not None              and te_elem.text else None

        sw = DEFAULT_SW
        sf = DEFAULT_SF
        n  = DEFAULT_N

        for met_elem in prior.findall('prior_metabolite'):
            name_elem = met_elem.find('name')
            if name_elem is None:
                continue
            name = name_elem.text.strip()

            # Parse all <line> elements
            peaks = []
            for line in met_elem.findall('line'):
                ppm_e   = line.find('ppm')
                area_e  = line.find('area')
                phase_e = line.find('phase')
                if ppm_e is None:
                    continue
                ppm   = float(ppm_e.text.strip())
                area  = float(area_e.text.strip()) if area_e is not None else 1.0
                phase = float(phase_e.text.strip()) if phase_e is not None else 0.0
                peaks.append((ppm, area, phase))

            if not peaks:
                continue

            # Synthesize FID from peak list
            fid = _synthesize_fid(peaks, sw, sf or 123.26, n)

            results.append({
                'fid'   : fid,
                'sw'    : sw,
                'sf'    : sf,
                'n'     : n,
                'name'  : name,
                'peaks' : peaks,
                'te'    : te,
                'source': path,
            })
            _log(f"  Loaded: {name}  {len(peaks)} peaks")

    return results


def _read_simulation_format(root, path):
    """
    Read VeSPA Simulation experiment export with ppms/areas/phases arrays.
    These are base64+zlib encoded numpy arrays.
    """
    results = []

    # Build metabolite id → name map from metabolite elements
    met_id_to_name = {}
    for met in root.findall('.//metabolite'):
        met_id  = met.get('id', '')
        name_el = met.find('name')
        if met_id and name_el is not None:
            met_id_to_name[met_id] = name_el.text.strip()

    sw = DEFAULT_SW
    sf = DEFAULT_SF
    n  = DEFAULT_N

    for simulation in root.findall('.//simulation'):
        met_id_el = simulation.find('metabolite_id')
        if met_id_el is None:
            continue

        met_id = met_id_el.text.strip() if met_id_el.text else ''
        name   = met_id_to_name.get(met_id, f'Metab_{len(results)+1:02d}')

        # Decode ppms, areas, phases
        ppms_el  = simulation.find('ppms')
        areas_el = simulation.find('areas')
        phases_el= simulation.find('phases')

        if ppms_el is None:
            continue

        try:
            ppms  = _decode_array(ppms_el)
            areas = _decode_array(areas_el) if areas_el is not None \
                    else np.ones(len(ppms))
            phases= _decode_array(phases_el) if phases_el is not None \
                    else np.zeros(len(ppms))
        except Exception as ex:
            _log(f"  WARNING: could not decode arrays for {name}: {ex}")
            continue

        peaks = list(zip(ppms.tolist(), areas.tolist(), phases.tolist()))
        fid   = _synthesize_fid(peaks, sw, sf or 123.26, n)

        results.append({
            'fid'   : fid,
            'sw'    : sw,
            'sf'    : sf,
            'n'     : n,
            'name'  : name,
            'peaks' : peaks,
            'source': path,
        })
        _log(f"  Loaded: {name}  {len(peaks)} peaks")

    return results


def _decode_array(elem):
    """Decode a VeSPA base64+zlib encoded numpy array."""
    import base64, zlib

    encoding = elem.get('encoding', '')
    dtype    = elem.get('data_type', 'float64')
    shape_s  = elem.get('shape', '')
    text     = (elem.text or '').strip()

    data = base64.b64decode(text)

    if 'zlib' in encoding:
        data = zlib.decompress(data)

    arr = np.frombuffer(data, dtype=np.dtype(dtype))

    if shape_s:
        try:
            shape = tuple(int(s) for s in shape_s.split(','))
            arr   = arr.reshape(shape)
        except Exception:
            pass

    return arr.ravel()


def _synthesize_fid(peaks, sw, sf, n, lw_hz=2.0, ppm_calib=DEFAULT_PPM_CALIB):
    """
    Synthesize a time-domain FID from a list of (ppm, area, phase) peaks.
    Uses Lorentzian line shapes.
    """
    dt  = 1.0 / sw
    t   = np.arange(n) * dt
    T2  = 1.0 / (np.pi * lw_hz)
    fid = np.zeros(n, dtype=complex)

    for ppm, area, phase in peaks:
        freq_hz = (ppm - ppm_calib) * sf   # Hz offset from center
        fid += area * np.exp(1j * (2 * np.pi * freq_hz * t + np.radians(phase))) \
                    * np.exp(-t / T2)

    return fid


############## Folder reader #############

def read_vespa_folder(folder_path):
    """
    Read all VeSPA .xml files in a folder.

    Parameters
    ----------
    folder_path : str

    Returns
    -------
    list of core struct dicts
    """
    folder_path = os.path.abspath(folder_path)
    xml_files   = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.lower().endswith('.xml') and
        not f.startswith('.') and not f.startswith('._')
    ])

    if not xml_files:
        raise RuntimeError(f"No .xml files found in {folder_path}")

    results = []
    for path in xml_files:
        try:
            results.extend(read_vespa(path))
        except Exception as ex:
            _log(f"  WARNING: could not read {os.path.basename(path)}: {ex}")

    if not results:
        raise RuntimeError(f"No valid VeSPA prior data found in {folder_path}")

    _log(f"\n  Total loaded: {len(results)} metabolites from {len(xml_files)} files")
    return results


############## Quick test #############

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python read_vespa.py <path_to_vespa.xml>")
        sys.exit(1)

    results = read_vespa(sys.argv[1])
    print(f"\nLoaded {len(results)} metabolites")
    for r in results:
        print(f"  {r['name']:15s}  peaks={len(r['peaks'])}  "
              f"n={r['n']}  sw={r['sw']}  sf={r['sf']}")
        for ppm, area, phase in r['peaks'][:3]:
            print(f"    ppm={ppm:.4f}  area={area:.4f}  phase={phase:.2f}")
