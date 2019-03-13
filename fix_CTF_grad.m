function [data] = fix_CTF_grad(headerpath, data)
% FIX_CTF_GRAD rebases the CTF coordinate system of the gradiometer positions
% on average subject head position during the experiment (instead of the 
% initial position), which will increase source localisation precision.
% Headposition must have been recorded continously in the CTF headtracking 
% channels. Has (so far) only been tested for CTF275 systems.
%
% Use as
%   data = ft_fix_grad(headerpath, data)
%
% headerpath is the file path of the header of the dataset in data, which
% will be read in using ft_read_header. We purposefully don't accept
% headers which have already been read in, as we need to read them in using
% dewar coordinates and want to avoid wrong coordinates being fed to the
% function.
% data must be 'raw' data, cut into trials.

% do the general setup of the function
ft_defaults

% This function only works for CTF275 MEGs so far.
ft_checkdata(data, 'senstype', 'ctf275', 'datatype', 'raw');

% Check prerequisites.
HLCs = ft_channelselection('HLC*', data);
assert(length(HLCs) == 9, ...
    'All 9 head localizer coil channels are required.');

%% Gather prerequisites
nTrials = size(data.trial, 2);
dewar_coord = ft_read_header(headerpath, 'coordsys', 'dewar');
dewar_coord = dewar_coord.grad;

% find headcoils
CNX_idx = strcmp('HLC0011', data.label);
CNY_idx = strcmp('HLC0012', data.label);
CNZ_idx = strcmp('HLC0013', data.label);
CLX_idx = strcmp('HLC0021', data.label);
CLY_idx = strcmp('HLC0022', data.label);
CLZ_idx = strcmp('HLC0023', data.label);
CRX_idx = strcmp('HLC0031', data.label);
CRY_idx = strcmp('HLC0032', data.label);
CRZ_idx = strcmp('HLC0033', data.label);

%calculate mean x,y,z position of each of the three headcoils during each
%trial
for t = 1:nTrials
coil1(:,t) = [mean(data.trial{1,t}(CNX_idx,:)); ...
              mean(data.trial{1,t}(CNY_idx,:)); ...
              mean(data.trial{1,t}(CNZ_idx,:))];
coil2(:,t) = [mean(data.trial{1,t}(CLX_idx,:)); ...
              mean(data.trial{1,t}(CLY_idx,:)); ...
              mean(data.trial{1,t}(CLZ_idx,:))];
coil3(:,t) = [mean(data.trial{1,t}(CRX_idx,:)); ...
              mean(data.trial{1,t}(CRY_idx,:)); ...
              mean(data.trial{1,t}(CRZ_idx,:))];
end
allC = [coil1;coil2;coil3];
assert(all(allC(:) >-1) && all(allC(:)<1) && range(allC(:))>0.1, ...
    'Headmovement data should be in the original meter unit.');

% Get mean coordinates of each coil over all trials
mN = mean(coil1, 2).*100; % Nasion
mL = mean(coil2, 2).*100; % LPA
mR = mean(coil3, 2).*100; % RPA

% Get transform matrix, apply it to chan- and coilpos and -ori
[h, coordsys] = ft_headcoordinates(mN, mL, mR, [], 'ctf');
newGrad = ft_transform_geometry(h, dewar_coord);
newGrad.coordsys = coordsys;

%% Replace grad field of the data structure with the headmovement-corrected
% grad.
ft_notice(['Replacing channel- and coil-orientation and -position with '...
    'headmotion-corrected ones.']);
data.grad = newGrad;

end