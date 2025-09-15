% Author: Sanjay Chaturvedi, Kanpur, India.
% 
% Paper: Unknown Input Observer Design for a Class of Nonlinear Systems: an LMI Approach,
% ACC 2006, Weitian Chen and Mehrdad Saif
% 
% NOTE: This simulation files gives the values which gives correct simulation results. I have solved the LMI using 
% the MATLAB's inbuilt LMI solver. The values given in paper do not give the results as presented in the paper.

% xdot = Ax + Bu + f(x) + Dv, y = cx ; u is input, v is disturbance, f(x) is nonlinear part, y is output
% zdot = Nz + Ly + Mf(xhat), xhat = z - Ey ; 

%%
clear; close; clc;

%% system matrices
A = [-1 -1 0;-1 0 0;0 -1 -1];
B = [0];
C = [1 0 0;0 0 1]; % size = 2x3
D = [-1;0;0];

%% step 1: U and V
U = -D*pinv(C*D); % size = 3x2
V = eye(size((C*D)*pinv(C*D))) - (C*D)*pinv(C*D);   %
gamma = 0.65;

%% lmi in equation 14 and 15
% lmi initialization
setlmis([]);

% define decision variables
P = lmivar(1,[size(A,1) 1]);
Ybar = lmivar(2,[size(U,1) size(V,1)]);
Kbar = lmivar(2,[size(U,1) size(V,1)]);

% lmi 1
lmiterm([-1 1 1 P],1,1);    % P>0

% lmi 2
I = eye(size(A));

lmiterm([2 1 1 P],((I + U*C)*A)',1,'s');
lmiterm([2 1 1 Ybar],1,(V*C*A),'s');
lmiterm([2 1 1 Kbar],-1,C,'s');
lmiterm([2 1 1 1],gamma,I);

lmiterm([2 1 2 P],sqrt(gamma),(I + U*C));
lmiterm([2 1 2 Ybar],sqrt(gamma),(V*C));

lmiterm([2 2 2 1],-I,1);

% final lmi object
LMIs = getlmis;

% solving lmi and lmi feasibility
[Tmin, Xfeas] = feasp(LMIs);

% solutions for the lmi variable, rewritten in matrix form
P = dec2mat(LMIs,Xfeas,P)
Ybar = dec2mat(LMIs,Xfeas,Ybar)
Kbar = dec2mat(LMIs,Xfeas,Kbar)

%% step 3: Y and K
Y = inv(P)*Ybar;
K = inv(P)*Kbar;

%% step 4: E, M, N, G and L | used in simulink file "Nonlinear_UIO.slx"
E = U + Y*V
M = eye(size(E*C)) + E*C
N = M*A - K*C
G = M*B
L = K*(eye(size(C*E)) + C*E) - M*A*E

%% simulation
z_intial = [3;3;3]; % [0.3,0.3,0.3]

%% running simulink model
model = 'Nonlinear_UIO';  
load_system(model);
sim(model);

% opening all scopes
scopes = find_system(model, 'BlockType', 'Scope');
for i = 1:length(scopes)
    set_param(scopes{i}, 'Open', 'on');
end

