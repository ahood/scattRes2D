function [abserr, relerr, diff_radial_chebfuns] = L2err(obj, f_values, other, g_values, other_theta)
% Computes the L2 error of f-g on B(0,R) using the values of f and g on 
% other_theta.
% Inputs:
%  - obj: a scattResComp2d or scattResComp2d_sparse object
%  - f_values: the values of a function f on mesh given by obj
%  - other: a scattResComp2d or scattResComp2d_sparse object with 
%           other.Rs(end) = obj.Rs(end)
%  - g_values: the values of a function g on mesh given by other
%  - other_theta: mesh of [0,2*pi] with other_theta(1) = 0,
%                 other_theta(end) = 2*pi

% TODO: check the other_theta vector to make sure the user followed
% instructions.

f_values_other_theta = mapToOtherThetaMesh(f_values,obj.Nr  ,other_theta);
g_values_other_theta = mapToOtherThetaMesh(g_values,other.Nr,other_theta);

f_rect = reshape(f_values_other_theta,obj.Nr  ,[]);
g_rect = reshape(g_values_other_theta,other.Nr,[]);

f_radial_chebfuns = cell(1,length(other_theta));
diff_radial_chebfuns = cell(1,length(other_theta));
for ii = 1:length(other_theta)
    f_chebfun = chebfun([flipud(f_rect(:,ii)); f_rect(:,ii)], [-obj.Rs(end), obj.Rs(end)]);
    g_chebfun = chebfun([flipud(g_rect(:,ii)); g_rect(:,ii)], [-other.Rs(end), other.Rs(end)]);
    f_radial_chebfuns{ii} = f_chebfun;
    diff_radial_chebfuns{ii} = f_chebfun - g_chebfun;
end
abserr = L2normOnDisc(diff_radial_chebfuns, other_theta);
relerr = abserr/L2normOnDisc(f_radial_chebfuns, other_theta);
