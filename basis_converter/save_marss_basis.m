function save_marss_basis(fid, metabolite_name, sw, sf, output_folder, varargin)
% SAVE_MARSS_BASIS  Save a Spinach (or any MATLAB) simulation result as a
%                   MARSS-compatible .mat file for visualization in INSPECTOR
%                   or use with the MRS Basis Set Converter.
%
% Usage:
%   save_marss_basis(fid, metabolite_name, sw, sf, output_folder)
%   save_marss_basis(fid, metabolite_name, sw, sf, output_folder, te)
%   save_marss_basis(fid, metabolite_name, sw, sf, output_folder, te, ppmCalib)
%
% Inputs:
%   fid             - complex array, the FID from Spinach or any simulator
%   metabolite_name - string, e.g. 'NAA', 'Cr', 'Cho'
%   sw              - spectral width in Hz (e.g. 4000)
%   sf              - Larmor frequency in MHz (e.g. 123.26 for 3T)
%   output_folder   - path to output folder (created if needed)
%   te              - (optional) echo time in ms (default: 0)
%   ppmCalib        - (optional) water reference ppm (default: 4.675)
%
% Output:
%   Saves <metabolite_name>.mat in output_folder in MARSS format.
%   Readable by INSPECTOR, MARSS, and the MRS Basis Set Converter.
%
% Example:
%   fid = liquid(spin_system, @press_sequence, parameters, 'nmr');
%   save_marss_basis(fid, 'NAA', 4000, 123.26, './marss_basis', 30);
%
% Authors: Kay Igwe

    if nargin < 5
        error('save_marss_basis: at least 5 arguments required.');
    end

    te        = 0;
    ppmCalib  = 4.675;

    if nargin >= 6, te       = varargin{1}; end
    if nargin >= 7, ppmCalib = varargin{2}; end

    % Flatten FID to column vector
    fid = fid(:);
    n   = length(fid);

    %%%%%%%%%%% Create output folder %%%%%%%%%%%
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('Created folder: %s\n', output_folder);
    end

    %%%%%%%%%%% Build MARSS-compatible exptDat struct %%%%%%%%%%%
    % MARSS stores each metabolite in exptDat with these fields:
    %   fid       - complex FID column vector
    %   sw_h      - spectral width in Hz
    %   sf        - Larmor frequency in MHz
    %   nspecC    - number of spectral points
    %   ppmCalib  - water reference frequency in ppm

    exptDat.fid      = fid;
    exptDat.sw_h     = sw;
    exptDat.sf       = sf;
    exptDat.nspecC   = n;
    exptDat.ppmCalib = ppmCalib;

    % Optional fields
    exptDat.te       = te;
    exptDat.name     = metabolite_name;
    exptDat.source   = 'Spinach';
    exptDat.date     = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    %%%%%%%%%%% Save %%%%%%%%%%%
    outpath = fullfile(output_folder, [metabolite_name '.mat']);
    save(outpath, 'exptDat', '-v7.3');

    fprintf('Saved: %s  (n=%d, sw=%.1f Hz, sf=%.4f MHz, te=%.1f ms, ppmCalib=%.3f)\n', ...
            [metabolite_name '.mat'], n, sw, sf, te, ppmCalib);

end
