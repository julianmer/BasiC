"""
readers/read_inspector.py
Basis Set Converter — INSPECTOR reader

Reads an INSPECTOR .mat basis set file into a list of core struct dicts.

Confirmed structure from Basis_sLASER_26.mat:
    Top-level variable: lcmBasis (NOT 'data')

    lcmBasis.sw_h     - spectral width in Hz         (e.g. 4000)
    lcmBasis.sf       - Larmor frequency in MHz       (e.g. 123.26)
    lcmBasis.ppmCalib - water reference in ppm        (ignored per project convention)
    lcmBasis.data     - object array, shape (n_metabs,)

    Each lcmBasis.data[i] is a 5-element object array:
        [0] : name      (str)           metabolite name e.g. 'Asc'
        [1] : linewidth (float)         simulation linewidth in Hz e.g. 1.5
        [2] : unknown   (float)         0.025 — purpose TBD
        [3] : fid       (complex128)    FID array shape (n,)
        [4] : []        (empty)         unused

Note: ppmCalib is present but ignored — derived when needed.

Core struct returned per metabolite:
    {
        'fid'       : np.ndarray, complex, shape (n,)
        'sw'        : float, spectral width in Hz
        'sf'        : float, Larmor frequency in MHz
        'n'         : int, number of points
        'name'      : str, metabolite name from data[i][0]
        'linewidth' : float, from data[i][1]
        'source'    : str, path to source file
    }

Entry point:
    read_inspector(path) -> list of core struct dicts
"""

import os
import numpy as np
import scipy.io as sio


def _log(*args): pass  # GUI handles display


def read_inspector(path):
    """
    Read an INSPECTOR .mat basis set file.

    Parameters
    ----------
    path : str
        Full path to an INSPECTOR .mat file.

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    raw = _load_mat(path)

    ############## Find lcmBasis struct #############
    basis = _find_struct(raw, path)

    ############## Extract top-level parameters #############
    sw_h = float(np.asarray(basis.get('sw_h', 0)).squeeze())
    sf   = float(np.asarray(basis.get('sf',   0)).squeeze())

    if sw_h == 0:
        raise RuntimeError(
            f"Missing or zero sw_h in {os.path.basename(path)}"
        )
    if sf == 0:
        raise RuntimeError(
            f"Missing or zero sf in {os.path.basename(path)}"
        )

    ############## Try flat arrays first (written by write_inspector.py) #############
    # These survive scipy roundtrip better than nested cell arrays
    if 'fids' in basis and 'names' in basis:
        fids_mat  = np.asarray(basis['fids'])
        names_arr = np.asarray(basis['names']).ravel()
        lws_arr   = np.asarray(basis.get('linewidths',
                        np.ones(len(names_arr)))).ravel()

        results = []
        for i in range(len(names_arr)):
            name = str(names_arr[i]).strip()
            fid  = fids_mat[:, i].astype(complex) if fids_mat.ndim == 2                    else fids_mat.ravel().astype(complex)
            lw   = float(lws_arr[i]) if i < len(lws_arr) else 1.5

            results.append({
                'fid'      : fid,
                'sw'       : sw_h,
                'sf'       : sf,
                'n'        : len(fid),
                'name'     : name,
                'linewidth': lw,
                'source'   : path,
            })
            _log(f"  Loaded: {name}")

        if results:
            _log(f"\n  Total loaded: {len(results)} metabolites")
            return results

    ############## Fallback: parse data cell array #############
    data = basis.get('data')
    if data is None:
        raise RuntimeError(
            f"No 'data' field found in lcmBasis of {os.path.basename(path)}"
        )

    data  = np.asarray(data)
    cells = data.ravel()
    results = []

    for i, cell in enumerate(cells):
        try:
            cell = np.asarray(cell).ravel()
            if len(cell) < 4:
                continue
            name      = str(cell[0]).strip()
            linewidth = float(np.real(cell[1])) if cell[1] is not None else 1.0
            fid       = np.asarray(cell[3]).squeeze().astype(complex)
            if fid.ndim != 1 or len(fid) == 0:
                continue
            results.append({
                'fid'      : fid,
                'sw'       : sw_h,
                'sf'       : sf,
                'n'        : len(fid),
                'name'     : name,
                'linewidth': linewidth,
                'source'   : path,
            })
            _log(f"  Loaded: {name}")
        except Exception as ex:
            _log(f"  WARNING: could not parse data[{i}]: {ex}")

    if not results:
        raise RuntimeError(
            f"No valid metabolites found in {os.path.basename(path)}"
        )

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {sw_h:.1f} Hz  sf = {sf:.4f} MHz  n = {results[0]['n']}")
    return results


############## Helpers #############

def _load_mat(path):
    try:
        return sio.loadmat(path, simplify_cells=True)
    except NotImplementedError:
        try:
            import mat73
            return mat73.loadmat(path)
        except ImportError:
            raise RuntimeError(
                f"HDF5 v7.3 file — install mat73: pip install mat73"
            )
    except Exception as ex:
        raise RuntimeError(f"Could not load {os.path.basename(path)}: {ex}")


def _find_struct(raw, path):
    """
    Find the lcmBasis struct in the loaded .mat file.
    Tries 'lcmBasis' first, then falls back to inspecting field names.
    """
    # Confirmed key from your file
    if 'lcmBasis' in raw:
        return raw['lcmBasis']

    # Some versions may use different names
    for key in ['lcmbasis', 'basis', 'LCMBasis']:
        if key in raw:
            return raw[key]

    # Check if any top-level key has 'data', 'sw_h', 'sf'
    for key, val in raw.items():
        if key.startswith('_'):
            continue
        if isinstance(val, dict):
            if {'data', 'sw_h', 'sf'}.issubset(val.keys()):
                return val

    raise RuntimeError(
        f"Could not find INSPECTOR lcmBasis struct in "
        f"{os.path.basename(path)}. Found keys: "
        f"{[k for k in raw if not k.startswith('_')]}"
    )
