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
SNRin = -3;
sigmaN = db2pow(-SNRin);
d = 3;
alpha = 0.998;
gamma = 0.044;
delta = 0.01;
trials = 200;
%% DIFFERENT NUMBER OF ANTENNAS
SINRout = zeros(trials,1);
steerVec = zeros(m,patPoints);
wTrials = zeros(d,trials);
%% STEERING VECTORS
for i = 1:patPoints
    for k = 1:m
        steerVec(k,i) = exp(1i*pi*cos(deg2rad(i))*(k-1));
    end
end
B = eye(m) - steerVec(:,desiredAngle)*steerVec(:,desiredAngle)'/...
        (steerVec(:,desiredAngle)'*steerVec(:,desiredAngle));
%% AVERAGE OVER TRIALS
for trial = 1:trials 
    fprintf("TRIAL: %g\n",trial)
    %% DEFINITIONS
    r = zeros(m,N);
    rDesired = zeros(m,N);
    rInterferencePlusNoise = zeros(m,N);
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
    %% RECEIVED DATA
    for k=1:m
         rDesired(k,:) = sDesired * ...
                    exp(1i*pi*(k-1)*cos(deg2rad(desiredAngle)));  
         noise = (randn(N,1)+1i*randn(N,1))*sqrt(sigmaN/2);       
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

    for i = 1:N
        Rhat = Rhat*(i-1)+r(:,i)*r(:,i)';
        Rhat = Rhat/i;

        Tr = zeros(m,d);
        Tr(:,1) = Rhat*steerVec(:,desiredAngle);
        Tr(:,1) = Tr(:,1)/norm(Tr(:,1));

        for d_prime = 2:d
           Tr(:,d_prime) = Rhat*B*Tr(:,d_prime-1);
           Tr(:,d_prime) = Tr(:,d_prime)/norm(Tr(:,d_prime));        
        end

        rb = B*r(:,i);
        rBar = Tr'*rb;
        y = gamma*steerVec(:,desiredAngle)'*r(:,i) - w'*rBar;
        xTilde = conj(y)*rBar;
        dTilde = gamma*conj(y)*steerVec(:,desiredAngle)'*r(:,i)-1;

        kTilde = Qinv*xTilde / (alpha+xTilde'*Qinv*xTilde);
        zetaTilde = dTilde - w'*xTilde;
        Qinv = Qinv/alpha - kTilde*xTilde'*Qinv;
        w = w + kTilde*conj(zetaTilde);
    end
    wTrials(:,trial) = w;
    clc
end
%% w
wAvg = mean(wTrials,2);
%% PATTERN
pattern = zeros(180,1);
augSteerVec = zeros(2*m,180);
for angle = 1:180
%     wDim = Tr*wAvg;
    w = gamma*steerVec(:,desiredAngle)-B*Tr*wAvg;
    pattern(angle) = abs(w'*steerVec(:,angle))^2;  
end
pattern = pow2db(pattern/max(pattern));
% pattern = pow2db(pattern);

plot(1:180, pattern,'LineWidth',1,"DisplayName","LCCM-KS") 
xlabel("DOA/Degree")
ylabel("Output Power(dB)")
title("Beampattern of the reduced-rank beamforming scheme")
legend show
grid on
hold on