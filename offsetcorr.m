function imgoc = offsetcorr(img)

% Load offset correction array ('omap') to workspace.
load('OffsetCorrection.mat');

% Crop omap to 128x128 centered square.
omap128 = omap(65:192,65:192);

% Convert image ('img') to single precision.
img_s = single(img);

% Perform offset correction by subtracting omap128 from raw image.
imgoc = img_s - omap128;
end