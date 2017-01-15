function f_values_other_theta = mapToOtherThetaMesh(f_values,Nr,other_theta)
% Maps a vector of values of a function f on a certain mesh onto another
% mesh which is finer in the theta direction (or coarser, but why would we
% do that?). 
% Inputs:
%  - f_values: values of f stacked like [f(r,theta1); f(r,theta2); ...; f(r,thetaN)]
%              where theta1 = 0.
%  - Nr: number of mesh points in the r direction
%  - other_theta: new theta mesh. 

f_values_rect = reshape(f_values,Nr,[]);
Nt = size(f_values_rect,2);

other_theta = other_theta(:); % make it a column vector

% map the values onto different theta mesh
c_rect = fft(f_values_rect.'); % columns should be values on circles
if mod(Nt,2) == 0
    U = [exp(1i*other_theta*(0:Nt/2-1)), cos(other_theta*Nt/2), ...
         exp(1i*other_theta*(-Nt/2+1:-1))]/Nt;
else
    N2 = (Nt-1)/2;
    U = [exp(1i*other_theta*(0:N2)), ...
         exp(1i*other_theta*(-N2:-1))]/Nt;
end
         
f_values_rect_other_theta = (U*c_rect).';
f_values_other_theta = f_values_rect_other_theta(:);
