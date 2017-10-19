%% ME 555 Final Project 
%------------------------

% Permission is NOT granted to duplicate any part of the code (unless for... 
% ... teaching purposes) without the consent of 'all' the members of the....
% ... group. Presenters are graduate students in the Department of....
% ... Mechanical and Nuclear Engineering at The Pennsylvania State....
% ...  University, University Park at the time of this presentation.

% Topics covered: 
% 1. Open loop step response
% 2. LQR full state feedback design
% 3. variation of closed loop poles with change in input weighting matrix R
% 4. simulation of pendulum angle (\theta) variation with variation in R
% 5. study of the effect of R in non-linear range of the system model
% 6. robustness study for 50% decrease and increase in parameters affecting
% ...pre-filter gain, kpf

clear all
close all
clc

disp('Carey Whitehair, caw68@psu.edu')
disp('Krishna Varadarajan, kzv3@psu.edu')
disp('Nicholas Conde, nmc39@psu.edu')
disap('Course Instructor: Dr. Chris Rahn')
disp('ME 555 Final Project Spring 2017: DYNAMIC ANALYSIS OF AN INVERTED PENDULUM USING LQR FULL STATE FEEDBACK CONTROL')

%% Open Loop Simulation Results:
m1o = 0.103; % kg
m2o = 0.785; % kg
mw1 = 0.11; % kg
mw2 = 1; % kg

m1 = m1o + mw1; % m
m2 = m2o + mw2; % m

g = 9.81; % m/s^2
lo = 0.33; % m
lco = 0.071; % m
lt = 0.07; % m
lb = 0.165; %m 
t = 0.014; % m 
lw2 = -(t + lt + lb)/2;
lc = (mw2*lw2 + m2o*lco)/m2;
Jost = 0.0246; 
Joe = Jost + m1*lo^2 + mw2*lw2^2;   
Jst = (Joe - (m1*lo^2));
kf = 0.0013;
kx = 50200;
ka = 2546;
ks=32;


%  system matrices after converting to encoder counts 
A = [0                             1 0                0;
    m2*lc*g/Jst                    0 m1*g*ka/(kx*Jst) 0; 
    0                              0 0                1;
    (Jst - m2*lo*lc)*g*kx/(ka*Jst) 0 -m1*lo*g/Jst     0];

B = ks*kf*(1/Jst)*[0; -lo*ka; 0; Joe*kx/m1]; % F(t) is the input

% only measuring \theta
C11 = 1;
C22 = 0;
C33 = 0;
C44 = 0;

C = [C11 0  0    0;
     0   C22 0    0;
     0   0   C33  0;
     0   0   0   C44];

C = [1 0 0 0];

D = 0;

sys = ss(A,B,C,D); % create the state space model - OPEN LOOP
x0 = [0;0;0;0]; % initial conditions
TF = 0.11; % time for which to run the simulation - similar to experimental plot --- % calcualted by observation by KP

figure(1) % open loop step
t = 0:0.001:8;
u = 400*sign(sin(2*pi*t/8));  
u(:,4002:end)=0;
lsim(sys,u,t,'g')
% [y,t] = step(sys,TF); % open loop system response
% plot(t,35.5*y(:,1)) % we are off by a factor of 35.5 compared to hardware
ylabel('\theta_{c} (counts)','FontSize',14);
xlabel('time (s)','FontSize',14);
title('Open-Loop Step Response');
grid on

%% LQR - full state feedback controller design

% Objective function - integral(\theta^2 + R*F(t)^2)

Clqr = [1 0 0 0]; % \theta as single output

Q = Clqr'*Clqr; % weighting matrix for \theta
R = 1:1:200; % a random value for the weighting on F(t) - b/w 1 and 200 (linear for 1-50: from manual)

% initialize a myriad of matrices
Klqr = zeros(200,4); 
eigLQR = zeros(4,200);
real_z1 = zeros(4,200);
imag_z1 = zeros(4,200);
kpf = zeros(200,1);

for k = 1:200
    Klqr(k,:) = lqr(A,B,Q,R(k)); % lqr gain
    kpf(k,1) = -m1*g/(ka*ks*kf) - kx*((m1*lo+m2*lc)/(m1*ka))*Klqr(k,3)+Klqr(k,1);  % pre-filter gain
    eigLQR(:,k) = eig(A-B*Klqr(k,:)); % closed loop poles
    
    real_z1(:,k) = real(eigLQR(:,k)); % real part
    imag_z1(:,k) = imag(eigLQR(:,k)); % imaginary part
end

figure(2) % close lopp poles Vs R
k = 1:200;
plot(real_z1(:,k),imag_z1(:,k),'r*')
xlabel('real part','FontSize',14);
ylabel('imaginary part','FontSize',14);
title({'variation of closed loop poles with change in R'; 'Observation: R increases, system is slower, lower gain, lower control effort'});
hold on;
annotation('arrow',[0.15 0.9],[0.515 0.515]) % add arrow to show the direction of increment of R
text(-12,0.5,'R increases','FontSize',14)
grid on

%%
% gain values separated according to R values
Klqrnl_low = [Klqr(1,:); Klqr(2,:); Klqr(3,:); Klqr(4,:)]; % R<5
kpflqrn1_low = [kpf(1,1); kpf(2,1); kpf(3,1); kpf(4,1)];

Klqr_perfect = Klqr(5,:); % R=5
kpflqr_perfect = kpf(5,1);

Klqrlinear = [Klqr(10,:); Klqr(15,:); Klqr(20,:); Klqr(25,:); Klqr(30,:); Klqr(35,:); Klqr(40,:); Klqr(45,:); Klqr(50,:)]; % 5< R< 50
kpflqrlinear  = [kpf(10,1); kpf(15,1); kpf(20,1); kpf(25,1); kpf(30,1); kpf(35,1); kpf(40,1); kpf(45,1); kpf(50,1)];

Klqr_nl_high = [Klqr(100,:); Klqr(200,:)]; % R>50
kpflqr_nl_high =[kpf(100,1); kpf(200,1)];

figure(3) % R<5
for k=1:1:4
    subplot(2,2,k)
    t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
    sysclrob1 = ss(A-B*Klqrnl_low(4,:),B,Clqr,[]);
    u = 400*sign(sin(2*pi*t/8));  
    u(:,4002:end)=0;
    lsim(kpflqrn1_low(k,:)*sysclrob1,u,t,'m')
    xlabel('Time (s)','FontSize',14);
    ylabel('Amplitude (Counts)','FontSize',14);
    title(['R = ' num2str(k)],'FontSize',10);
    grid on
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the control input weighting in non-linear range of the system'});

figure(4) % R=5    
t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
sysclrob1 = ss(A-B*Klqr_perfect,B,Clqr,[]);
u = 400*sign(sin(2*pi*t/8));  
u(:,4002:end)=0;
lsim(kpflqr_perfect*sysclrob1,u,t,'g')
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude (Counts)','FontSize',14);
title(['R = ' num2str(5)],'FontSize',10);
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the control input weighting in linear range of the system'});
grid on


figure(5) % 5<R<50
for k=1:1:9
    subplot(3,3,k)
    t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
    sysclrob1 = ss(A-B*Klqrlinear(k,:),B,Clqr,[]);
    u = 400*sign(sin(2*pi*t/8));  
    u(:,4002:end)=0;
    lsim(kpflqrlinear(k,:)*sysclrob1,u,t,'b')
    xlabel('Time (s)','FontSize',14);
    ylabel('Amplitude (Counts)','FontSize',14);
    title(['R = ' num2str(k*5+5)],'FontSize',10);
    grid on
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the control input weighting in the linear range of the system'});

figure(6) % R>50
for k=1:1:2
    subplot(2,1,k)
    t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
    sysclrob1 = ss(A-B*Klqr_nl_high(k,:),B,Clqr,[]);
    u = 400*sign(sin(2*pi*t/8));  
    u(:,4002:end)=0;
    lsim(kpflqr_nl_high(k,:)*sysclrob1,u,t,'r')
    xlabel('Time (s)','FontSize',14);
    ylabel('Amplitude (Counts)','FontSize',14);
    title(['R = ' num2str(k*100)],'FontSize',10);
    grid on
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the control input weighting in non-linear range of the system'});


%% CASE: Varying dwell time for R = 5,25

figure(7) % R=5

tvar = 2*[2,8,16]; % dwell time of 2,8,16 sec

for k=1:3
subplot(2,2,k)
t = 0:0.001:tvar(k); % 2 secs dwell, 400 count step, and goes back to 0
syscl7 = ss(A-B*Klqr(5,:),B,Clqr,[]);
u = 400*sign(sin(2*pi*t/tvar(k)));  
    if k==1
        zero = 2002;
    end
    if k ==2
        zero = 8002;
    end
    if k == 3
        zero = 16002;
    end
u(:,zero:end)=0;
lsim(kpf(5,1)*syscl7,u,t,'r')
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude (Counts)','FontSize',14);
title(['R = 5, dwell time = ' num2str(tvar(k)/2) 's'],'FontSize',10);
grid on
  
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the dwell time for R = 5'});

figure(8) % R=25

tvar = 2*[2,8,16]; % dwell time of 2,8,16 sec

for k=1:3
subplot(2,2,k)
t = 0:0.001:tvar(k); % 4 secs dwell, 400 count step, and goes back to 0
syscl8 = ss(A-B*Klqr(25,:),B,Clqr,[]);
u = 400*sign(sin(2*pi*t/tvar(k))); 
if k==1
        zero = 2002;
    end
    if k ==2
        zero = 8002;
    end
    if k == 3
        zero = 16002;
    end
u(:,zero:end)=0;
lsim(kpf(25,1)*syscl8,u,t,'r')
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude (Counts)','FontSize',14);
title(['R = 25, dwell time = ' num2str(tvar(k)/2) 's'],'FontSize',10);
grid on
   
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the in the dwell time for R = 25'});


%% CASE: Varying step size for R = 25


figure(9)

amp = [100,600];

for k=1:1:2
    subplot(2,1,k)
    t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
    sysclrob1 = ss(A-B*Klqr(25,:),B,Clqr,[]);
    u = amp(k)*sign(sin(2*pi*t/8));  
    u(:,4002:end)=0;
    lsim(kpf(25,1)*sysclrob1,u,t,'r')
    xlabel('Time (s)','FontSize',14);
    ylabel('Amplitude (Counts)','FontSize',14);
    title(['R = 25, step size = ' num2str(amp(k)) 'counts'],'FontSize',10);
    grid on
end
suptitle({'ME 555 Final Project', 'Dynamic analysis of inverted pendulum angle with variation in the step size of the system'});

%% Robustness check: reduce by half the parameters affecting pre-filter gain
% Reduced by 50%

m1onew = 0.103*0.5; % kg
m2onew = 0.785*0.5; % kg
mw1new = 0.11*0.5; % kg
mw2new = 1*0.5; % kg

m1new = m1onew + mw1new; % m
m2new = m2onew + mw2new; % m

lonew = 0.33*0.5; % m
lcnew = (mw2new*lw2 + m2onew*lco)/m2new;

Klqrnew = lqr(A,B,Q,5); % lqr gain, R =5
kpfnew = -m1new*g/(ka*ks*kf) - kx*((m1new*lonew+m2new*lcnew)/(m1new*ka))*Klqrnew(1,3)+Klqrnew(1,1);  % pre-filter gain

figure(10) % R=5    
t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
sysclrob1 = ss(A-B*Klqrnew,B,Clqr,[]);
u = 400*sign(sin(2*pi*t/8));  
u(:,4002:end)=0;
lsim(kpfnew*sysclrob1,u,t,'g')
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude (Counts)','FontSize',14);
title('R = 5, 50% parameter reduction','FontSize',10);


% Increase by 50%

m1onew2 = m1o + m1o*0.5; % kg
m2onew2 = m2o + m2o*0.5; % kg
mw1new2 = mw1 + mw1*0.5; % kg
mw2new2 = mw2 + mw2*0.5; % kg

m1new2 = m1onew2 + mw1new2; % m
m2new2 = m2onew2 + mw2new2; % m

lonew2 = lo + lo*0.5; % m
lcnew2 = (mw2new2*lw2 + m2onew2*lco)/m2new2;

Klqrnew2 = lqr(A,B,Q,5); % lqr gain, R =5
kpfnew2 = -m1new2*g/(ka*ks*kf) - kx*((m1new2*lonew2+m2new2*lcnew2)/(m1new2*ka))*Klqrnew2(1,3)+Klqrnew2(1,1);  % pre-filter gain

figure(11) % R=5    
t = 0:0.001:8; % 4 secs dwell, 400 count step, and goes back to 0
sysclrob2 = ss(A-B*Klqrnew2,B,Clqr,[]);
u = 400*sign(sin(2*pi*t/8));  
u(:,4002:end)=0;
lsim(kpfnew2*sysclrob2,u,t,'g')
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude (Counts)','FontSize',14);
title('R = 5, 50% parameter increase','FontSize',10);

disp('End of simulation')  
% End of code

