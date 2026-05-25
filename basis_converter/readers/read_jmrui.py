"""
readers/read_jmrui.py
Basis Set Converter — jMRUI reader

Reads a folder of jMRUI .txt basis function files into core struct dicts.

jMRUI .txt file structure (confirmed from FID-A io_writejmrui.m):

    jMRUI Data Textfile

    Filename: <path>

    PointsInDataset: <n>
    DatasetsInFile: 1
    SamplingInterval: <dwelltime_ms>     <- in MILLISECONDS
    ZeroOrderPhase: 0.000000E+00
    BeginTime: 0.000000E+00
    TransmitterFrequency: <txfrq_hz>     <- in HZ
    MagneticField: <Bo>                  <- in Tesla
    TypeOfNucleus: 0.000000E+00
    NameOfPatient: <name>
    DateOfExperiment: <date>
    Spectrometer: <scanner>
    AdditionalInfo: <info>


    Signal and FFT
    sig(real)   sig(imag)   fft(real)   fft(imag)
    Signal 1 out of 1 in file
    <real>  <imag>  <fft_real>  <fft_imag>
    ...

Key units:
    SamplingInterval is in MILLISECONDS → sw = 1000 / SamplingInterval
    TransmitterFrequency is in Hz → sf = TransmitterFrequency / 1e6

Core struct returned per metabolite:
    {
        'fid'    : np.ndarray, complex, shape (n,)
        'sw'     : float, spectral width in Hz
        'sf'     : float, Larmor frequency in MHz
        'n'      : int, number of points
        'name'   : str, from NameOfPatient or filename stem
        'Bo'     : float, field strength in Tesla
        'source' : str, path to source .txt file
    }

Entry points:
    read_jmrui_file(path)       → one core struct dict
    read_jmrui_folder(path)     → list of core struct dicts
"""

import os
import numpy as np


def _log(*args): pass  # GUI handles display

def read_jmrui_file(path, name=None):
    """
    Read one jMRUI .txt basis function file.

    Parameters
    ----------
    path : str
        Full path to a jMRUI .txt file.
    name : str, optional
        Override metabolite name. If None, read from NameOfPatient header.

    Returns
    -------
    core struct dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r') as f:
        lines = f.readlines()

    ############## Parse header #############
    header = {}
    data_start = None

    for i, line in enumerate(lines):
        line_stripped = line.strip()

        # Look for key: value pairs
        if ':' in line_stripped and not line_stripped.startswith('Signal'):
            parts = line_stripped.split(':', 1)
            key = parts[0].strip()
            val = parts[1].strip()
            header[key] = val

        # Data starts after the column header line
        if line_stripped.startswith('Signal') and 'out of' in line_stripped:
            data_start = i + 1
            break

    if data_start is None:
        raise RuntimeError(
            f"Could not find data section in {os.path.basename(path)}"
        )

    ############## Extract header values #############
    n = int(float(header.get('PointsInDataset', 0)))
    if n == 0:
        raise RuntimeError(
            f"PointsInDataset is 0 in {os.path.basename(path)}"
        )

    # SamplingInterval is in MILLISECONDS
    dwell_ms = float(header.get('SamplingInterval', 0))
    if dwell_ms == 0:
        raise RuntimeError(
            f"SamplingInterval is 0 in {os.path.basename(path)}"
        )
    sw = 1000.0 / dwell_ms    # convert ms → Hz

    # TransmitterFrequency is in Hz
    txfrq_hz = float(header.get('TransmitterFrequency', 0))
    sf = txfrq_hz / 1e6       # convert Hz → MHz

    Bo = float(header.get('MagneticField', 0))

    if name is None:
        name = header.get('NameOfPatient', '').strip()
        if not name:
            name = os.path.splitext(os.path.basename(path))[0]

    ############## Parse data #############
    # 4 columns: sig(real)  sig(imag)  fft(real)  fft(imag)
    # We only need columns 0 and 1 (the FID)
    real_vals = []
    imag_vals = []

    for line in lines[data_start:]:
        cols = line.strip().split()
        if len(cols) < 2:
            continue
        try:
            real_vals.append(float(cols[0]))
            imag_vals.append(float(cols[1]))
        except ValueError:
            continue

    if len(real_vals) == 0:
        raise RuntimeError(
            f"No data found after header in {os.path.basename(path)}"
        )

    fid = np.array(real_vals) + 1j * np.array(imag_vals)

    # Trim to n if longer
    if len(fid) > n:
        fid = fid[:n]
    elif len(fid) < n:
        pad = np.zeros(n - len(fid), dtype=complex)
        fid = np.concatenate([fid, pad])

    return {
        'fid'    : fid,
        'sw'     : sw,
        'sf'     : sf,
        'n'      : len(fid),
        'name'   : name,
        'Bo'     : Bo,
        'source' : path,
    }


def read_jmrui_folder(folder_path):
    """
    Read all jMRUI .txt basis function files in a folder.

    Parameters
    ----------
    folder_path : str

    Returns
    -------
    list of core struct dicts
    """
    folder_path = os.path.abspath(folder_path)

    if not os.path.isdir(folder_path):
        raise RuntimeError(f"Folder not found: {folder_path}")

    txt_files = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.lower().endswith('.txt')
    ])

    if not txt_files:
        raise RuntimeError(f"No .txt files found in: {folder_path}")

    # Check first file looks like jMRUI format
    with open(txt_files[0], 'r') as f:
        first_line = f.readline().strip()
    if 'jMRUI' not in first_line:
        raise RuntimeError(
            f"First .txt file does not appear to be jMRUI format "
            f"(expected 'jMRUI Data Textfile', got: '{first_line}')"
        )

    results = []
    failed  = []

    for path in txt_files:
        try:
            core = read_jmrui_file(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(f"No valid jMRUI files loaded from: {folder_path}")

    if failed:
        _log(f"\n  {len(failed)} file(s) skipped: {', '.join(failed)}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {results[0]['sw']:.1f} Hz  "
          f"sf = {results[0]['sf']:.4f} MHz  "
          f"n = {results[0]['n']}")
    return results
