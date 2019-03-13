function [data] = fix_CTF_grad(headerpath, data)
% FIX_CTF_GRAD rebases the CTF coordinate system of the gradiometer positions
% on average subject head position during the experiment (instead of the 
% initial position), which will increase source localisation precision.
% Headposition must have been recorded continously in the CTF headtracking 
% channels. Has (so far) only been tested for CTF275 systems.
%
% For clarity: The chanori- and coiloir-fields have unit vector coordinates
% from chanpos/coilpos (shifted to the origin) towards the orientation of 
% the channel/coil.
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

%% Find unit vectors of the new coordinate system, project old coordinates
% onto them.

% The X-Axis is the unit normal of the vector from the origin (midpoint 
% between LPA and RPA) to the nasion.
origin = mean([mL mR],2);
xNorm = (mN - origin) / norm(mN - origin);

% The Z-Axis is in the direction of the unit normal of the plane spanned by
% the three fiducial coils - for an explanation, see
% http://mathworld.wolfram.com/Point-PlaneDistance.html
% Compute the unit normal from the XY-Plane
zNorm = (cross((mL - mN), (mR - mN))) / ...
    norm(cross((mL - mN), (mR - mN)));

% The Y-Axis is perpendicular to both X- and Z_axis, so we can just use
% their crossproduct. The norm isn't really necessary, but let's keep it
% for consistency. The Y-Axis in CTF is defined as going to the left 
% - so we flip it from the cross product, which goes to the right.
yNorm = - (cross(xNorm, zNorm) / norm(cross(xNorm, zNorm)));

% Apply new coordinate system to channel positions and orientations.
for iChan = 1:length(dewar_coord.chanpos)
    iChanpos = dewar_coord.chanpos(iChan,:)';
    iChanori = dewar_coord.chanori(iChan,:)';
    
    chanPos(iChan,1) = dot(xNorm, (iChanpos - origin));
    chanPos(iChan,2) = dot(yNorm, (iChanpos - origin));
    chanPos(iChan,3) = dot(zNorm, (iChanpos - origin));
        
    chanOri(iChan,1) = dot(xNorm, iChanori);
    chanOri(iChan,2) = dot(yNorm, iChanori);
    chanOri(iChan,3) = dot(zNorm, iChanori);
end

% Apply new coordinate system to coil positions and orientations.
for iCoil = 1:length(dewar_coord.coilpos)
    iCoilpos = dewar_coord.coilpos(iCoil,:)';
    iCoilori = dewar_coord.coilori(iCoil,:)';
    
    coilPos(iCoil,1) = dot(xNorm, (iCoilpos - origin));
    coilPos(iCoil,2) = dot(yNorm, (iCoilpos - origin));
    coilPos(iCoil,3) = dot(zNorm, (iCoilpos - origin));
        
    coilOri(iCoil,1) = dot(xNorm, iCoilori);
    coilOri(iCoil,2) = dot(yNorm, iCoilori);
    coilOri(iCoil,3) = dot(zNorm, iCoilori);
end
    
%% Replace grad field of the data structure with the headmovement-corrected
% grad.
ft_notice(['replacing channel- and coil-orientation and -position with '...
    'headmotion-corrected ones']);
newGrad = data.grad;
newGrad.chanpos = chanPos;
newGrad.chanori = chanOri;
newGrad.coilpos = coilPos;
newGrad.coilori = coilOri;

data.grad = newGrad;

end

