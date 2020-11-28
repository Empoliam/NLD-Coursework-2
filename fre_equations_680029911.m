%%
% System of ODEs for the firing rate equations from Montbrio et al 2015
% with sinusoidal input current
% fixed parameters
Jbar=15; eta=-5;  omega=pi; 

%% Set your personal parameter value for EL
% insert your student number xxxyyyzzz (usually starts with a 6)
% rename as fre_equations_xxxyyyzzz.m and call this in your scripts for
% each question
SNumber=680029911;

if SNumber==0
   error(['Enter your 9 digit student number in your local copy of fre_equations.m then rename as instructed\n  ',...
'See instruction sheet or contact James Rankin on Teams or at j.a.rankin@exeter.ac.uk if unsure.']); 
end

rng(SNumber)
Delta=0.95+0.1*rand;
phase=2*pi*rand;


disp('Your personal parameter values:');
disp(['SNumber=',num2str(SNumber)])
format longg
disp(['Delta=',num2str(Delta,'%.6f')])
disp(['phase (phi)=',num2str(phase,'%.6f')])
format short

%%y=[r;v]
%Define the FRQ model, modified from Montbrio et al 2015 equations E1
Gamma=sqrt(Delta);
F=@(t,y,eta,Jbar,a) [((Delta/pi)+(Gamma+a*sin(t*omega+phase)).*y(1,:)/pi+2*y(1,:).*y(2,:)); ...
    (y(2,:).^2+eta+Jbar*y(1,:)-(pi^2)*(y(1,:).^2))];

% to use with your functions:
rhs=@(y,a,t) F(t,y,eta,Jbar,a); 


