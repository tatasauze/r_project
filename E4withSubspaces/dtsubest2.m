% desviaciones tipica

function dts = dtsubest2(theta,din,z);

[F, A, V, G] = thd2ech(theta, din);
[theta, din, lab] = arma2thd(F(:,5:8), [], A(:,5:8), [], V, 1);

theta = [theta zeros(size(theta,1),1)];
for k = 1:size(theta,1)
    if theta(k,1) == 1 | theta(k,1) == 0;
        theta(k,2) = 1;
    end
end

std = imod(theta, din, z, 1);
dts = zeros(size(theta,1),1);
dts(find(theta(:,2)~=1)) = std;