function theta = saff_sphere(N)
% 'Uniformly' distribute N points on the surface of a sphere. theta(:,1) is
% the measured from the pole(0-pi) and theta(:,2) is ranges over (0..2*pi)
% Function written by John Enright (date unknown).

k=1:N;
h = -1 + 2*(k-1)/(N-1);

theta = zeros(N,2);
theta(:,1) = acos(h);
for m = 2:(N-1)
    theta(m,2) = mod((theta(m-1,2) + (3.6/sqrt(N))*1/sqrt(1-h(m).^2)),2*pi);
end