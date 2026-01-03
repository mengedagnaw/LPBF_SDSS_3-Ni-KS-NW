# the engine, full code:


function rowT = KSNW_core(rootDir, condFolder, condLabel, outDir, tolDeg, phaseBCC, phaseFCC, grainAngDeg, doPNGs, mtexDir)
% KSNW_core
% Robust KS/NW interphase-boundary analysis for MTEX (EBSD .ang)
% Uses ROTATIONS to avoid BCC/FCC symmetry mismatch in orientation comparisons.
%
% OUTPUT:
% rowT table with: Condition, File, PhaseFCC,nFCC, PhaseBCC,nBCC, nGBseg, FracKS_tol, FracNW_tol, tol_deg

  arguments
    rootDir (1,1) string
    condFolder (1,1) string
    condLabel (1,1) string
    outDir (1,1) string
    tolDeg (1,1) double = 5
    phaseBCC (1,1) string = "Ferryt"
    phaseFCC (1,1) string = "Austenite"
    grainAngDeg (1,1) double = 5
    doPNGs (1,1) logical = true
    mtexDir (1,1) string = ""
  end

  angFile = fullfile(rootDir, condFolder, "0.ang");

  fprintf("\n=============================\n");
  fprintf("Condition: %s\n", condLabel);
  fprintf("Folder   : %s\n", fullfile(rootDir, condFolder));
  fprintf("File     : %s\n", angFile);

  if ~exist(angFile,'file')
    warning("File not found. Skipping.");
    rowT = local_empty_row(condLabel, angFile, tolDeg);
    return
  end
  if ~exist(outDir,'dir'); mkdir(outDir); end

  % --- Setup MTEX (optional)
  local_setup_mtex(mtexDir);

  % --- Crystal symmetries (adjust lattice params if you want)
  % IMPORTANT: We will compare OR via ROTATION objects, so mismatch is avoided.
  csBCC = crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', char(phaseBCC));
  csFCC = crystalSymmetry('m-3m', [3.600 3.600 3.600], 'mineral', char(phaseFCC));
  cs = {csBCC, csFCC};

  % --- Load EBSD
  try
    ebsd = EBSD.load(angFile, cs, 'interface','ang', 'convertEuler2SpatialReferenceFrame');
  catch
    ebsd = loadEBSD_ang(angFile, cs, 'convertEuler2SpatialReferenceFrame');
  end

  ebsd = ebsd('indexed');
  if isempty(ebsd)
    warning("No indexed EBSD. Skipping.");
    rowT = local_empty_row(condLabel, angFile, tolDeg);
    return
  end

  % --- Pick phase names robustly (your logs show BCC="Ferryt", FCC="Austenite")
  [bccName, fccName] = local_pick_phases(ebsd, phaseBCC, phaseFCC);

  ebsdBCC = ebsd(char(bccName));
  ebsdFCC = ebsd(char(fccName));
  nBCC = numel(ebsdBCC);
  nFCC = numel(ebsdFCC);

  fprintf("Detected FCC: %s (%d pts)\n", fccName, nFCC);
  fprintf("Detected BCC: %s (%d pts)\n", bccName, nBCC);

  if nBCC == 0 || nFCC == 0
    warning("Missing one phase (BCC or FCC). Skipping.");
    rowT = local_empty_row(condLabel, angFile, tolDeg);
    return
  end

  % --- Grains (two-phase)
  [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', grainAngDeg*degree);
  grains = smooth(grains, 5);

  % --- Interphase boundaries (BCC–FCC)
  gB = grains.boundary(char(bccName), char(fccName));
  nGB = numel(gB);

  if nGB == 0
    warning("No BCC–FCC boundary segments found.");
    rowT = local_empty_row(condLabel, angFile, tolDeg);
    rowT.nFCC = nFCC; rowT.nBCC = nBCC; rowT.nGBseg = 0;
    return
  end

  w = local_segment_weights(gB);

  % --- Measured OR as ROTATION
  gid = gB.grainId;  % Nx2
  oriBCC = grains(gid(:,1)).meanOrientation;
  oriFCC = grains(gid(:,2)).meanOrientation;
  rotMeas = inv(rotation(oriBCC)) .* rotation(oriFCC);

  % --- KS / NW variants as ROTATIONS
  varsKS = local_make_KS_variants(csBCC, csFCC);
  varsNW = local_make_NW_variants(csBCC, csFCC);

  devKS = local_min_dev_deg(rotMeas, varsKS);
  devNW = local_min_dev_deg(rotMeas, varsNW);

  fracKS = sum(w(devKS <= tolDeg)) / sum(w);
  fracNW = sum(w(devNW <= tolDeg)) / sum(w);

  fprintf("K-S: weighted fraction within %.2f deg = %.3f\n", tolDeg, fracKS);
  fprintf("N-W: weighted fraction within %.2f deg = %.3f\n", tolDeg, fracNW);

  % --- PNG exports
  if doPNGs
    local_export_maps_and_hists(condLabel, outDir, grains, gB, devKS, devNW, tolDeg);
  end

  % --- Row output
  rowT = table(string(condLabel), string(angFile), string(fccName), nFCC, string(bccName), nBCC, nGB, ...
               fracKS, fracNW, tolDeg, ...
    'VariableNames', {'Condition','File','PhaseFCC','nFCC','PhaseBCC','nBCC','nGBseg','FracKS_tol','FracNW_tol','tol_deg'});
end

%% ========================= Helper functions =========================

function local_setup_mtex(mtexDir)
  % If user passes mtexDir, use it. Otherwise do nothing (assumes MTEX already initialized).
  if strlength(mtexDir) == 0
    return
  end
  if exist(mtexDir,'dir')
    % remove other MTEX versions if any (best-effort)
    p = string(strsplit(path, pathsep));
    kill = contains(lower(p), "mtex-") | contains(lower(p), filesep+"mtex"+filesep);
    for k = find(kill)
      rmpath(char(p(k)));
    end
    addpath(genpath(char(mtexDir)));
    try
      startup_mtex;
    catch
      % some installs use "startup_mtex.m" but path already OK anyway
    end
  end
end

function rowT = local_empty_row(condLabel, filePath, tolDeg)
  rowT = table(string(condLabel), string(filePath), "", 0, "", 0, 0, NaN, NaN, tolDeg, ...
    'VariableNames', {'Condition','File','PhaseFCC','nFCC','PhaseBCC','nBCC','nGBseg','FracKS_tol','FracNW_tol','tol_deg'});
end

function [bccName, fccName] = local_pick_phases(ebsd, bccWanted, fccWanted)
  bccWanted = string(bccWanted);
  fccWanted = string(fccWanted);

  % Direct hit?
  if ~isempty(ebsd(char(bccWanted))) && ~isempty(ebsd(char(fccWanted)))
    bccName = bccWanted; fccName = fccWanted; return
  end

  names = local_phase_names(ebsd);

  % Aliases
  bccAliases = ["ferryt","ferrite","alpha","α","bcc"];
  fccAliases = ["austenite","austenit","gamma","γ","fcc"];

  bccName = local_first_match(names, bccAliases, bccWanted);
  fccName = local_first_match(names, fccAliases, fccWanted);

  % Safety: if still empty, fall back to wanted strings
  if bccName == "" ; bccName = bccWanted; end
  if fccName == "" ; fccName = fccWanted; end
end

function names = local_phase_names(ebsd)
  names = string.empty;

  % Try mineralList (newer MTEX)
  try
    names = string(ebsd.mineralList);
  catch
  end

  % Try minerals directly (often works even when other APIs fail)
  if isempty(names)
    try
      names = unique(string(ebsd.mineral));
    catch
      names = string.empty;
    end
  end

  names = unique(names);
  names = names(names~="" & lower(names)~="notindexed");
end

function out = local_first_match(names, aliases, fallback)
  out = "";
  ln = lower(names);

  % 1) exact fallback
  if any(names == fallback)
    out = fallback; return
  end

  % 2) substring alias match
  for a = aliases
    hit = names(contains(ln, lower(a)));
    if ~isempty(hit), out = hit(1); return; end
  end

  % 3) if nothing found, return fallback (may still work)
  out = fallback;
end

function w = local_segment_weights(gB)
  try
    w = gB.segLength;
    if isempty(w), w = ones(numel(gB),1); end
  catch
    w = ones(numel(gB),1);
  end
  w = w(:);
end

function vars = local_make_KS_variants(csBCC, csFCC)
  % KS: {110}bcc || {111}fcc and <111>bcc || <110>fcc
  h_b = Miller(1,1,0, csBCC, 'hkl');
  h_f = Miller(1,1,1, csFCC, 'hkl');
  d_b = Miller(1,-1,1, csBCC, 'uvw');
  d_f = Miller(1,-1,0, csFCC, 'uvw');
  R0 = rotation.map(h_b, h_f, d_b, d_f);
  vars = local_expand_variants(R0, csBCC, csFCC);
end

function vars = local_make_NW_variants(csBCC, csFCC)
  % NW: {110}bcc || {111}fcc and <001>bcc || <110>fcc
  h_b = Miller(1,1,0, csBCC, 'hkl');
  h_f = Miller(1,1,1, csFCC, 'hkl');
  d_b = Miller(0,0,1, csBCC, 'uvw');
  d_f = Miller(1,-1,0, csFCC, 'uvw');
  R0 = rotation.map(h_b, h_f, d_b, d_f);
  vars = local_expand_variants(R0, csBCC, csFCC);
end

function Rvars = local_expand_variants(R0, csBCC, csFCC)
  % Use MTEX "symmetries(...)" which is stable; avoids cs.rotations problems.
  sB = symmetries(csBCC);
  sF = symmetries(csFCC);

  Rall = rotation.empty; idx = 0;
  for i = 1:numel(sB)
    for j = 1:numel(sF)
      idx = idx + 1;
      Rall(idx,1) = sB(i) * R0 * inv(sF(j));
    end
  end
  Rvars = local_unique_rotations(Rall, 1e-6);
end

function Runiq = local_unique_rotations(R, tolRad)
  R = R(:);
  keep = true(numel(R),1);
  for i = 1:numel(R)
    if ~keep(i), continue; end
    for j = i+1:numel(R)
      if keep(j) && angle(R(i)*inv(R(j))) < tolRad
        keep(j) = false;
      end
    end
  end
  Runiq = R(keep);
end

function devDeg = local_min_dev_deg(rotMeas, vars)
  rotMeas = rotation(rotMeas);
  vars = vars(:);

  n = numel(rotMeas);
  dev = inf(n,1);
  for k = 1:numel(vars)
    dRot = rotMeas .* inv(vars(k));
    dev  = min(dev, angle(dRot)./degree);
  end
  devDeg = dev;
end

function local_export_maps_and_hists(condLabel, outDir, grains, gB, devKS, devNW, tolDeg)
  % Make figures without toolbars/menus to avoid "axes toolbar" export noise
  baseBoundaryColor = [0.80 0.80 0.80];

  % KS map
  f = figure('Visible','off','Color','w','ToolBar','none','MenuBar','none');
  ax = axes(f); %#ok<NASGU>
  plot(grains.boundary, 'lineColor', baseBoundaryColor, 'lineWidth', 0.5); hold on;
  plot(gB(devKS<=tolDeg), 'lineColor','k', 'lineWidth', 1.5);
  title(sprintf('%s — KS hits (<= %.1f°)', condLabel, tolDeg));
  axis tight off;
  local_export_png(fullfile(outDir, sprintf('%s_map_KS.png', condLabel)));

  % NW map
  f = figure('Visible','off','Color','w','ToolBar','none','MenuBar','none');
  ax = axes(f); %#ok<NASGU>
  plot(grains.boundary, 'lineColor', baseBoundaryColor, 'lineWidth', 0.5); hold on;
  plot(gB(devNW<=tolDeg), 'lineColor','k', 'lineWidth', 1.5);
  title(sprintf('%s — NW hits (<= %.1f°)', condLabel, tolDeg));
  axis tight off;
  local_export_png(fullfile(outDir, sprintf('%s_map_NW.png', condLabel)));

  % KS histogram
  f = figure('Visible','off','Color','w','ToolBar','none','MenuBar','none');
  histogram(devKS, 0:0.25:20);
  xlabel('Deviation to nearest KS variant (deg)'); ylabel('Count');
  title(sprintf('%s — KS deviation histogram', condLabel));
  local_export_png(fullfile(outDir, sprintf('%s_hist_KS.png', condLabel)));

  % NW histogram
  f = figure('Visible','off','Color','w','ToolBar','none','MenuBar','none');
  histogram(devNW, 0:0.25:20);
  xlabel('Deviation to nearest NW variant (deg)'); ylabel('Count');
  title(sprintf('%s — NW deviation histogram', condLabel));
  local_export_png(fullfile(outDir, sprintf('%s_hist_NW.png', condLabel)));
end

function local_export_png(pngPath)
  set(gcf,'Color','w'); drawnow;
  try
    exportgraphics(gcf, pngPath, 'Resolution', 300);
  catch
    print(gcf, pngPath, '-dpng', '-r300');
  end
  close(gcf);
end
