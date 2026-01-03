rootDir = "/MATLAB Drive/SDSS 2507_EBSD";
condFolder = "03_SR450_hex";
outDir = rootDir + "/OR_KS_NW_outputs/SR450";
tolDeg = 5;

opts = struct(); opts.grainAngDeg = 5; opts.doMaps = true; opts.doHists = true;

T = KSNW_core(rootDir, condFolder, outDir, tolDeg, opts);
disp(T);
