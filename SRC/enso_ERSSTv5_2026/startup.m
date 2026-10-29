% Resolve the shared library relative to SRC/enso_ERSSTv5_2026.
analysisDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(analysisDir));
addpath(genpath(fullfile(projectDir, 'nlsa')));
clear analysisDir projectDir
beep off
