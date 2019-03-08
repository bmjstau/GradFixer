function [data] = ft_fix_grad(cfg, data, header)
% FT_FIX_GRAD rebases the CTF coordinate system of the gradiometer positions
% on average subject head position during the experiment (instead of the 
% initial position), which will increase source localisation precision.
% data must be 'raw'. Headposition must have been recorded continously in
% the CTF headtracking channels.
%
% Use as
%   data = ft_fix_grad(cfg, data, header)
%
% The configuration can contain
% cfg.headerfile
%
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile     = string with the filename
%   cfg.headerfile   = string with the filename

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    elec_original
ft_preamble provenance elec_original
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% This function only works for CTF275 MEGs so far.
ft_checkdata(data, 'senstype', 'ctf275', 'datatype', 'raw');
% Check prerequisites.
assert(isfield(data, 'grad'), 'grad positions are required');
HLCs = ft_channelselection('HLC*', data);
assert(length(HLCs) == 9, ...
    'All 9 head localizer coil channels are required.');

%% Gather prerequisites
nTrials = size(data.trial, 2);
grad = data.grad;
dewar_coord = ft_read_header(subject.megdata, 'coordsys', 'dewar');
dewar_coord = dewar_coord.grad;

% find headcoils
CNX_idx = find(strcmp('HLC0011', data.label));
CNY_idx = find(strcmp('HLC0012', data.label));
CNZ_idx = find(strcmp('HLC0013', data.label));
CLX_idx = find(strcmp('HLC0021', data.label));
CLY_idx = find(strcmp('HLC0022', data.label));
CLZ_idx = find(strcmp('HLC0023', data.label));
CRX_idx = find(strcmp('HLC0031', data.label));
CRY_idx = find(strcmp('HLC0032', data.label));
CRZ_idx = find(strcmp('HLC0033', data.label));

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
% for consistency.
yNorm = cross(xNorm, zNorm) / norm(cross(xNorm, zNorm));

for iChan = 1:length(dewar_coord.chanpos)
    chan_x(iChan) = dot(xNorm, (dewar_coord.chanpos(iChan,:)'- origin));
    chan_y(iChan) = dot(yNorm, (dewar_coord.chanpos(iChan,:)'- origin));
    chan_z(iChan) = dot(zNorm, (dewar_coord.chanpos(iChan,:)'- origin));
end
    
figure
subplot(3,1,1)
scatter3(dewar_coord.chanpos(:,1),dewar_coord.chanpos(:,2),dewar_coord.chanpos(:,3),20, chan_x, 'filled')
colorbar
subplot(3,1,2)
scatter3(dewar_coord.chanpos(:,1),dewar_coord.chanpos(:,2),dewar_coord.chanpos(:,3),20, chan_y, 'filled')
colorbar
subplot(3,1,3)
scatter3(dewar_coord.chanpos(:,1),dewar_coord.chanpos(:,2),dewar_coord.chanpos(:,3),20, chan_z, 'filled')
colorbar

end

