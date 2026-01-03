rootDir = "/MATLAB Drive/SDSS 2507_EBSD";
condFolder = "04_SR500_hex";
outDir = rootDir + "/OR_KS_NW_outputs/SR500";
tolDeg = 5;

opts = struct(); opts.grainAngDeg = 5; opts.doMaps = true; opts.doHists = true;

T = KSNW_core(rootDir, condFolder, outDir, tolDeg, opts);
disp(T);
