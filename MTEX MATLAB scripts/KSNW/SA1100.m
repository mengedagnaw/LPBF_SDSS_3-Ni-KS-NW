rootDir = "/MATLAB Drive/SDSS 2507_EBSD";
condFolder = "06_SA1100_hex";
outDir = rootDir + "/OR_KS_NW_outputs/SA1100";
tolDeg = 5;

opts = struct(); opts.grainAngDeg = 5; opts.doMaps = true; opts.doHists = true;

T = KSNW_core(rootDir, condFolder, outDir, tolDeg, opts);
disp(T);
