clear

syms f f0 mu s2
assume(s2>0)

gForceMF = 1/sqrt(2*pi*s2) * exp(-1/2 * (f- mu)^2 /s2 );

hxx = f0^4 /(36*s2^4) *( ((f-mu)^2 - 3*s2)^2 - 6*s2^2 );
Hxx = f0^4 /(36*s2^4) *( ((f-mu)^2 - 3*s2)^2 - 6*s2^2 ) *gForceMF ;
% Mixed
a =  ( (f-mu)^2 - 3*s2  ) * (11*(f-mu) + 27/40 *f0 * ( (f-mu)^2 -s2)/s2 );
b = 27/20*f0 *(f-mu)^2;
hxy = f0^3/120 /s2^3  * (a- b); hxy = simplify(hxy);
Hxy = hxy*gForceMF;
% Second diagonal
alpha = - f0/s2 * (11^2/20+27*11/400*f0 *(f-mu)/s2 + 27^2/16000*f0^2 *(f-mu)^2/s2^2);
beta =  11* (f-mu)/s2 +27/40*f0/s2^2 *((f-mu)^2 -3*s2);
delta = 11*(f-mu)+27/40*f0/s2 *((f-mu)^2 - s2);
hyy =  f0/20 * (f0/20/s2 *beta *delta + alpha); hyy = simplify(hyy);
Hyy = hyy *gForceMF;

% Third derivatives
Hxxx = 1/3*f0^2*diff(Hxx,s2);
Hxxx = simplify(Hxxx);

Hxxy = 11/20*f0*diff(Hxx,mu) + 27/400*f0^2*diff(Hxx,s2);
Hxxy = simplify(Hxxy);

Hyyy = 11/20*f0*diff(Hyy,mu) + 27/400*f0^2*diff(Hyy,s2);
Hyyy = simplify(Hyyy);

%Fourth derivatives
Hxxxx = 1/3*f0^2*diff(Hxxx,s2);
Hxxxx = simplify(Hxxxx);

Hxxyx = 1/3*f0^2*diff(Hxxy,s2);
Hxxyx = simplify(Hxxyx);

Hxxyy = 11/20*f0*diff(Hxxy,mu) + 27/400*f0^2*diff(Hxxy,s2);
Hxxyy = simplify(Hxxyy);

Hyyyx = 1/3*f0^2*diff(Hyyy,s2);
Hyyyx = simplify(Hyyyx);

Hyyyy = 11/20*f0*diff(Hyyy,mu) + 27/400*f0^2*diff(Hyyy,s2);
Hyyyy = simplify(Hyyyy);


hxxxx = simplify( Hxxxx / gForceMF );
hxxyx = simplify( Hxxyx / gForceMF );
hxxyy = simplify( Hxxyy / gForceMF );
hyyyx = simplify( Hyyyx / gForceMF );
hyyyy = simplify( Hyyyy / gForceMF ); 


% Write the symbolic expression as functions


matlabFunction(hxx,"File","get_hxx","Vars",[f f0 mu s2]);
matlabFunction(hxy,"File","get_hxy","Vars",[f f0 mu s2]);
matlabFunction(hyy,"File","get_hyy","Vars",[f f0 mu s2]);

matlabFunction(hxxxx,"File","get_hxxxx","Vars",[f f0 mu s2]);
matlabFunction(hxxyx,"File","get_hxxyx","Vars",[f f0 mu s2]);
matlabFunction(hxxyy,"File","get_hxxyy","Vars",[f f0 mu s2]);
matlabFunction(hyyyx,"File","get_hyyyx","Vars",[f f0 mu s2]);
matlabFunction(hyyyy,"File","get_hyyyy","Vars",[f f0 mu s2]);






