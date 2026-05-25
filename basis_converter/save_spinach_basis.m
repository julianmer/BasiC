function save_spinach_basis(fid, metabolite_name, sw, sf, output_folder, varargin)
% SAVE_SPINACH_BASIS  Save a Spinach simulation result as a FID-A compatible
%                     .mat file for use with the MRS Basis Set Converter.
%
% Usage:
%   save_spinach_basis(fid, metabolite_name, sw, sf, output_folder)
%   save_spinach_basis(fid, metabolite_name, sw, sf, output_folder, te)
%
% Inputs:
%   fid             - complex array, the FID returned by Spinach liquid()
%                     or equivalent simulation function
%   metabolite_name - string, name of the metabolite (e.g. 'NAA', 'Cr')
%   sw              - spectral width in Hz (e.g. 4000)
%   sf              - Larmor frequency in MHz (e.g. 123.26 for 3T)
%   output_folder   - path to output folder (will be created if needed)
%   te              - (optional) echo time in ms (default: 0)
%
% Output:
%   Saves a .mat file named <metabolite_name>.mat in output_folder
%   in FID-A compatible format, readable by the MRS Basis Set Converter.
%
% Example (after running a Spinach PRESS simulation):
%   % Run Spinach simulation
%   fid = liquid(spin_system, @press_sequence, parameters, 'nmr');
%
%   % Save for converter
%   save_spinach_basis(fid, 'NAA', 4000, 123.26, './spinach_basis', 30);
%
% After saving all metabolites, load the output_folder in the
% MRS Basis Set Converter (Browse folder) to convert to any format.
%
% The converter will detect the files as FID-A format.
%
% Authors: Kay Igwe

    %%%%%%%%%%% Input validation %%%%%%%%%%%
    if nargin < 5
        error('save_spinach_basis: at least 5 arguments required.');
    end

    te = 0;
    if nargin >= 6
        te = varargin{1};
    end

    % Flatten FID to column vector
    fid = fid(:);
    n   = length(fid);

    % ── Create output folder ──────────────────────────────────────────────────
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('Created folder: %s\n', output_folder);
    end

    % ── Build FID-A compatible struct ─────────────────────────────────────────
    % FID-A stores data in a struct with these fields.
    % The converter reads: specs, fids, sz, spectralwidth, txfrq, te, name
    
    out = struct();

    % Time-domain FID
    out.fids  = fid;

    % Frequency-domain spectrum (FFT of FID)
    out.specs = fftshift(fft(fid));

    % Dimensions
    out.sz    = [n, 1, 1, 1];

    % Spectral parameters
    out.spectralwidth = sw;          % Hz
    out.txfrq         = sf * 1e6;   % MHz → Hz
    out.dwelltime     = 1.0 / sw;   % seconds
    out.ppm           = (sf * 1e6) / (sw / n) * ...
                        ((0:n-1) - n/2) / (sf * 1e6);   % approximate ppm axis

    % Sequence parameters
    out.te    = te;                  % ms
    out.tr    = 0;
    out.Bo    = sf / 42.577;         % Tesla
    out.name  = metabolite_name;

    % Flags expected by FID-A readers
    out.flags.averaged    = 1;
    out.flags.addedrcvrs  = 1;
    out.flags.subtracted  = 0;
    out.flags.writtentotext = 0;
    out.flags.gotparams   = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered    = 0;
    out.flags.zeropadded  = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.isISIS      = 0;

    % Source info
    out.sim_source = 'Spinach';
    out.sim_date   = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % ── Save ──────────────────────────────────────────────────────────────────
    outpath = fullfile(output_folder, [metabolite_name '.mat']);
    save(outpath, 'out', '-v7.3');

    fprintf('Saved: %s  (n=%d, sw=%.1f Hz, sf=%.4f MHz, te=%.1f ms)\n', ...
            [metabolite_name '.mat'], n, sw, sf, te);

end
