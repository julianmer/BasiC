"""
core/formats.py
Basis Set Converter — format and tool compatibility definitions

Central registry mapping each detected format to:
    - The software tools that can use that format
    - A short display label
    - The file extension(s)

Used by:
    detect.py     — adds 'compatible_tools' to detection result
    screen1_upload.py — displays compatible tools in the result card
    screen2_output.py — groups tools by format compatibility
"""

############## Format → compatible tools mapping #############
#
# Key   : detected format string (must match exactly what detect.py returns)
# Value : list of tool names that can read/use that format

FORMAT_TOOLS = {

    ############## MARSS #############
    'MARSS .mat': [
        'MARSS',
    ],

    'MARSS .raw': [
        'MARSS',
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    ############## LCModel #############
    'LCModel .BASIS': [
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    'LCModel .raw': [
        'LCModel',
        'SPANT',
        'TARQUIN',
    ],

    ############## Osprey #############
    'Osprey .mat': [
        'Osprey',
    ],

    ############## FSL-MRS #############
    'FSL-MRS .json': [
        'FSL-MRS',
    ],

    ############## INSPECTOR #############
    'INSPECTOR .mat': [
        'INSPECTOR',
    ],

    ############## MRSCloud #############
    'MRSCloud .mat': [
        'MRSCloud',
    ],

    ############## FID-A #############
    'FID-A .mat': [
        'FID-A',
    ],

    ############## NIfTI-MRS #############
    'NIfTI-MRS': [
        'FSL-MRS',
        'Osprey',
        'SPANT',
    ],

    ############## jMRUI / NMRScopeB / QUEST #############
    'jMRUI .txt': [
        'jMRUI',
        'QUEST',
        'NMRScopeB',
    ],

    ############## ProFit #############
    'ProFit .mat': [
        'ProFit',
    ],

    ############## PyAMARES / OXSA #############
    'PyAMARES .csv': [
        'PyAMARES',
        'OXSA',
    ],

    ############## VeSPA #############
    'VeSPA .xml': [
        'VeSPA',
    ],
}


def get_compatible_tools(fmt):
    """
    Return the list of compatible tools for a detected format string.
    Returns empty list if format is not in the registry.

    Parameters
    ----------
    fmt : str
        Detected format string e.g. 'LCModel .BASIS'

    Returns
    -------
    list of str
    """
    return FORMAT_TOOLS.get(fmt, [])


def get_tools_display(fmt):
    """
    Return a human-readable string of compatible tools.
    e.g. 'LCModel · SPANT · TARQUIN'

    Parameters
    ----------
    fmt : str

    Returns
    -------
    str
    """
    tools = get_compatible_tools(fmt)
    if not tools:
        return ''
    return '  ·  '.join(tools)
