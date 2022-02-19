clear
% close
clc
%% SETTINGS
N = 1000;
m = 20;
desiredAngle = 50;
interfereAngle = [40 70 20 80];
jammAngle = [60 90 30 10];
patPoints = 180;
dBSNRin = -12:2;
SNRin = 10.^(dBSNRin/10);
d = 3;
% alpha = 0.998;
% gamma = 0.044;
% delta = 0.1;
alpha = 0.998;
gamma = 0.057;
delta = 0.01;
trials = 200;

SINRout = zeros(length(dBSNRin), trials);
steerVec = zeros(m,patPoints);
%% STEERING VECTORS
for i = 1:patPoints
    for k = 1:m
        steerVec(k,i) = exp(1i*pi*cos(deg2rad(i))*(k-1));
    end
end
aTilde = 1/sqrt(2)* ...
         [steerVec(:,desiredAngle).', steerVec(:,desiredAngle)'].';
%% B MATRIX
B = eye(m) - steerVec(:,desiredAngle)*steerVec(:,desiredAngle)'/...
    (steerVec(:,desiredAngle)'*steerVec(:,desiredAngle));
B1 = eye(m) - steerVec(:,desiredAngle)*steerVec(:,desiredAngle)'/...
     (2*steerVec(:,desiredAngle)'*steerVec(:,desiredAngle));
B2 = -steerVec(:,desiredAngle)*steerVec(:,desiredAngle).'/...
     (2*steerVec(:,desiredAngle)'*steerVec(:,desiredAngle));
Btilde = [B1, B2; conj(B2), conj(B1)];
%% AVERAGE OVER TRIALS
for trial = 1:trials    
    fprintf("TRIAL: %g (%g PERCENT DONE)\n",trial,(trial-1)/trials*100)
    %% DEFINITIONS
    r = zeros(m,N);
    rDesired = zeros(m,N);
    rInterferencePlusNoise = zeros(m,N);
    noiseIdx = 0;
    %% DATA GENERATION
    sDesired = sign(rand(N,1)-0.5);
    sInterfere1 = sign(rand(N,1)-0.5);
    sInterfere2 = sign(rand(N,1)-0.5);
    sInterfere3 = sign(rand(N,1)-0.5);
    sInterfere4 = sign(rand(N,1)-0.5);
    sJamm1 = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    sJamm2 = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    sJamm3 = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    sJamm4 = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    
    varDesired = var(sDesired);
    sigmaNu = varDesired./SNRin;
    %% DIFFERENT INPUT SNR
    for noiseVar = sigmaNu
        noiseIdx = noiseIdx + 1;
        %% RECEIVED DATA
        for k=1:m
             rDesired(k,:) = sDesired * ...
                        exp(1i*pi*(k-1)*cos(deg2rad(desiredAngle))); 
             noise = (randn(N,1)+1i*randn(N,1))*sqrt(noiseVar/2);       
             rInterferencePlusNoise(k,:) =  ...
                 sInterfere1 * ...
                    exp(1i*pi*(k-1)* cos(deg2rad(interfereAngle(1)))) + ...
                 sInterfere2 * ...
                    exp(1i*pi*(k-1)*cos(deg2rad(interfereAngle(2)))) + ...
                 sInterfere3 * ...
                    exp(1i*pi*(k-1)*cos(deg2rad(interfereAngle(3)))) + ...
                 sInterfere4 * ...
                    exp(1i*pi*(k-1)*cos(deg2rad(interfereAngle(4)))) + ...
                 sJamm1 * exp(1i*pi*(k-1)*cos(deg2rad(jammAngle(1)))) + ...
                 sJamm2 * exp(1i*pi*(k-1)*cos(deg2rad(jammAngle(2)))) + ...
                 sJamm3 * exp(1i*pi*(k-1)*cos(deg2rad(jammAngle(3)))) + ...
                 sJamm4 * exp(1i*pi*(k-1)*cos(deg2rad(jammAngle(4)))) + ...
                 + noise;
             r(k,:) = rDesired(k,:) + rInterferencePlusNoise(k,:); 
        end
        %% ALGORITHM
        w = [1;zeros(d-1,1)];
        Qinv = delta*eye(d);
        Rhat = zeros(m,m);
        RcHat = zeros(m,m);
        for i = 1:N
            Rhat = Rhat*(i-1)+r(:,i)*r(:,i)';
            RcHat = RcHat*(i-1)+r(:,i)*r(:,i).';
            Rhat = Rhat/i;
            RcHat = RcHat/i;

            P = zeros(m,d);
            P(:,1) = 0.5*Rhat*steerVec(:,desiredAngle) + ...
                     0.5*RcHat*conj(steerVec(:,desiredAngle));
            P(:,1) = P(:,1)/norm(P(:,1));

            for d_prime = 2:d
               P(:,d_prime) = 0.5*Rhat*(B1*P(:,d_prime-1)+B2*conj(P(:,d_prime-1))) + ...
                        0.5*RcHat*conj(B1*P(:,d_prime-1)+B2*conj(P(:,d_prime-1)));
               P(:,d_prime) = P(:,d_prime)/norm(P(:,d_prime));        
            end

            rBar = real((B1*P + B2*conj(P))'*r(:,i));
            y = gamma*real(steerVec(:,desiredAngle)'*r(:,i)) - w'*rBar;
            xTilde = conj(y)*rBar;
            dTilde = gamma*conj(y)*real(steerVec(:,desiredAngle)'*r(:,i))-1;

            kTilde = Qinv*xTilde / (alpha+xTilde'*Qinv*xTilde);
            zetaTilde = dTilde - w'*xTilde;
            Qinv = Qinv/alpha - kTilde*xTilde'*Qinv;
            w = w + kTilde*conj(zetaTilde);
        end
        %% APPROXIMATION OF AUTO-CORRELATION FUNCTION
        Rdesired = zeros(m,m);
        RcDesired = zeros(m,m);
        Rinterference = zeros(m,m);
        RcInterference = zeros(m,m);
        for i = 1:N
            Rdesired = Rdesired + rDesired(:,i)*rDesired(:,i)';
            RcDesired = RcDesired + rDesired(:,i)*rDesired(:,i).';
            Rinterference = Rinterference + ...
                            rInterferencePlusNoise(:,i)* ...
                            rInterferencePlusNoise(:,i)';
            RcInterference = RcInterference + ...
                             rInterferencePlusNoise(:,i)* ...
                             rInterferencePlusNoise(:,i).';
        end
        Rdesired = Rdesired/N;
        RcDesired = RcDesired/N;
        Rinterference = Rinterference/N;
        RcInterference = RcInterference/N;

        Rs = 0.5*[Rdesired, RcDesired; ...
                  conj(RcDesired), conj(Rdesired)];
        Rin = 0.5*[Rinterference, RcInterference; ...
                   conj(RcInterference), conj(Rinterference)];
        %% SINR OUT
        TrTilde = 1/sqrt(2)*[P.', P'].';
        wTilde = gamma*aTilde - Btilde*TrTilde*w;
        SINRout(noiseIdx,trial) = real((wTilde'*Rs*wTilde)/(wTilde'*Rin*wTilde));
    end
    clc
end
SINR = pow2db(mean(SINRout,2));

plot(dBSNRin,SINR,'LineWidth',1,"DisplayName","WLCCM-KS") 
xlabel("Input SNR")
ylabel("Output SINR (dB)")
title("Steady-state SINR versus the input SNR (N = 1000)")
legend show
grid on
hold on