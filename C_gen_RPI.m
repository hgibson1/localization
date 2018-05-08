function [C, WW, LL] = C_gen_RPI(W, L, H, Imax, stepsize, fixture, sensor, spot_angle)
%   Generate unsubsampled C matrix with quincunx f-s pairs
    %[fixture,sensor] = config_RPI(H);
    Nf=size(fixture,2); % # of fixtures
    Ns=size(sensor,2); % # of sensors
    
    WW=floor(W/stepsize);
    LL=floor(L/stepsize);
    C=zeros(Nf*Ns,WW*LL);
    ds=stepsize^2;
    for i=1:Nf
        for j=1:Ns
            fp=fixture(:,i);
            sp=sensor(:,j);
            for k=1:WW
                for l=1:LL
                    xp = [k-1/2;l-1/2;0]*stepsize;
                    D1_sq = sum((fp - xp) .^ 2);
                    D2_sq = sum((sp - xp) .^ 2);
                    cos_theta1 = fp(3) / sqrt(D1_sq);
                    theta1 = acos(cos_theta1);
                    cos_theta2_sq = (sp(3)^2) / D2_sq;
                    if spot_angle == -1 % spot_angle controls whether to use the q function or use the spot angle only
                        intensity = q_real(theta1);
                    elseif spot_angle == -2 % Use the q function of Photon LED
                        intensity = q_Photon(theta1);
                    elseif spot_angle == -3 
                        intensity = q_Photon2(theta1);
                    else
                        intensity = (theta1 <= spot_angle * pi / 180);
                    end
                    C((i-1) * Ns + j, (k-1) * LL + l) = Imax * ds * (intensity * cos_theta1 * cos_theta2_sq) / (D1_sq * D2_sq);
                end
            end
            %disp(['Fixture ',num2str(i),', sensor ',num2str(j),' finished!']);
        end
        %disp(['Fixture ',num2str(i),' finished!']);
    end
end

