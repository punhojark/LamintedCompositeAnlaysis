
function Laminate(In,Out)

clc

Out='output.txt';

%Read Data from input.txt
[Approach] = textread('input.txt', '%u', 1, 'headerlines', 2);
[F1t,F2t,F1c,F2c,Fs] = textread('input.txt', '%f %f %f %f %f', 1, 'headerlines', 6);
[Del_C,Del_T] = textread('input.txt', '%f %f', 1, 'headerlines', 10);
[Nxx,Nyy,Nxy,Mxx,Myy,Mxy] = textread('input.txt', '%f %f %f %f %f %f', 1, 'headerlines', 14);
[theta,th,Ef,Em,Gf,Gm,vf,vm,Vf,Alp_f,Alp_m,Beta_1,Beta_2,X] = textread('input.txt', '%f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'headerlines', 21);
Vm = 1-Vf;

%Common Approach
if Approach == 0
    E1 = Ef.*Vf + Em.*Vm;
    E2 = (Ef.*Em)./(Em.*Vf+Ef.*Vm);
    v12 = (vf.*Vf) + (vm.*Vm);
    G12 = (Gf.*Gm)./(Gm.*Vf+Gf.*Vm);

%elasticity Approach    
elseif Approach == 1
    v12 = vf.*Vf + vm.*Vm + (Vf.*Vm.*(vf-vm).*(2.*Ef.*(vm.*vm)+vm.*Ef-Ef+Em-Em.*vf-2.*Em.*vf.*vf))./((2.*(vm.*vm).*Vf-vm+vm.*Vf-1-Vf).*Ef+(2.*(vf.*vf)-Vf.*vf-2.*Vf.*(vf.*vf)+Vf+vf-1).*Em);
    E1 = Ef.*Vf + Em.*Vm - ((2.*Em.*Ef.*Vf.*((vf-vm).*(vf-vm)).*Vm)./(Ef.*(2.*(vm.*vm).*Vf-vm+Vf.*vm-Vf-1)+Em.*(-1-2.*Vf.*(vf.*vf)+vf-Vf.*vf+2.*(vf.*vf)+Vf)));
    
    Im = -4.*vm+3;
    If = -4.*vf+3;
    
    a1 = 3.*Vf.*(Vm.*Vm).*(Gf./Gm-1).*(Gf./Gm+If)+((Gf.*Im./Gm+If.*Im-(Gf.*Im./Gm-If)).*Vf.*Vf.*Vf).*(Vf.*Im.*(Gf./Gm-1)-(Gf.*Im./Gm+1));
    b1 = -3.*Vf.*(Vm.*Vm).*(Gf./Gm-1).*(Gf./Gm+If)+(1/2).*(Im.*Gf./Gm+(Gf./Gm-1).*Vf+1).*((Im-1).*(Gf./Gm+If)-2.*(Gf.*Im./Gm-If).*(Vf.*Vf.*Vf))+(Vf./2).*(Im+1).*(Gf./Gm-1).*(Gf./Gm+If+(Gf.*Im./Gm-If).*(Vf.*Vf.*Vf));
    c1 = 3.*Vf.*(Vm.*Vm).*(Gf./Gm-1).*(Gf./Gm+If)+(Im.*Gf./Gm+(Gf./Gm-1).*Vf+1).*(Gf./Gm+If+(Gf.*Im./Gm-If).*(Vf.*Vf.*Vf)); 
    L = -b1.*b1+4.*a1.*c1;
    q = (-b1+sqrt(L))./(2.*a1);
    G23 = q.*Gm;
    Kf = Ef./(2.*(vf+1).*(-2.*vf+1));
    Km = Em./(2.*(vm+1).*(-2.*vm+1));
    Kstar = (Km.*(Kf+Gm).*Vm + Kf.*(Km+Gm).*Vf)./((Kf+Gm).*Vm+(Km+Gm).*Vf);
    m = 4.*Kstar.*(v12.*v12)./E1+1;
    v23 = (Kstar-m.*G23)./(Kstar+m.*G23);
    
    E2 = 2.*(v23+1).*G23;
    G12 = Gm.*((Gf.*(Vf+1)+Gm.*Vm)./(Gf.*Vm+Gm.*(Vf+1)));

%Semi-empirical Approach
else
    I = ((Ef./Em)-1)./((Ef./Em)+X);
    
    E1 = (Ef.*Vf) + (Em.*Vm);
    E2 = Em.*(1+X.*I.*Vf)./(1-I.*Vf);
    v12 = (vf.*Vf) + (vm.*Vm);
    G12 = Gm.*(1+X.*I.*Vf)./(1-I.*Vf);
    
end


epsMax1t = F1t/E1;
epsMax1c = F1c/E1;
epsMax2t = F2t/E2;
epsMax2c = F2c/E2;
epsMaxs = Fs/G12;


Alpha_1 = (Alp_f.*Ef./E1).*Vf + (Alp_m.*Em./E1).*Vm;
Alpha_2 = (1+vf).*Alp_f.*Vf + (1+vm).*Alp_m.*Vm - Alpha_1.*v12;
%Beta_1 = ((Beta_f*Del_Cf*Vf*Ef + Beta_m*Del_Cm*Vm*Em)*Den_c)/(E1(Del_Cf*Den_f*Vf + Del_Cm*Den_m*Vm));
%Beta_2 = ((Vf(1+vf)*Del_Cf*Beta_f + Vm(1+vm)*Del_m*Beta_m)*(Den_c - Beta_1*v12))/(Vm*Den_m*Del_Cm + Vf*Den_f*Del_f);
Beta_1 = Beta_1;
Beta_2 = Beta_2;


AlpAlp2 = (Alp_m.*Em.*Vm + Alp_f.*Ef.*Vf)./(Em.*Vm + Ef.*Vf) 
AlpAlp1 = (1+vm).*Alp_m.*Vm + Alp_f.*Vf


N = [Nxx ; Nyy ; Nxy];
M = [Mxx ; Myy ; Mxy];

A=0;B=0;D=0;NT=0;NH=0;MT=0;MH=0;

%store h values in a matrix
h=zeros(size(th,1)+1,1);
h(1,1)=sum(th)/2;
for i=2:size(th,1)+1
    h(i)=h(i-1)-th(i-1);
end
h



%Calculate ABD matrix, Qbar
k=1;
disp(['A, B, & D Matrices for [' num2str(theta') '] Laminate:'])
for i=size(th,1):-1:1
    Alpha = [Alpha_1(i);Alpha_2(i);0];
    Beta =  [Beta_1(i);Beta_2(i);0];
    S12 =[1/E1(i,1) -v12(i,1)/E1(i,1) 0; -v12(i,1)/E1(i,1) 1/E2(i,1) 0; 0 0 1/G12(i,1)];
    
    S12
    
    Q12 = inv(S12);
    disp(['Lamina # ' int2str(k) ':' num2str(theta(i,1))])
    eval(['Qbar' int2str(k) '= convert_sig(-theta(' int2str(i) '))*Q12*convert_eps(theta(' int2str(i) '))']);
    Qbar=eval(['Qbar' int2str(k)]);
    
    Sbar = inv(Qbar);
    A = A - Qbar * (h(i+1,1) - h(i,1));                 % Negative summation since i=size(t,1):-1:1
    B = B - 1/2 * Qbar * (h(i+1,1)^2 - h(i,1)^2);
    D = D - 1/3 * Qbar * (h(i+1,1)^3 - h(i,1)^3); 

    NT = NT + Del_T * Qbar * convert_eps(-theta(i)) * Alpha * (h(i,1) - h(i+1,1));
    NH = NH + Del_C * Qbar * convert_eps(-theta(i)) * Beta  * (h(i,1) - h(i+1,1));
    MT = MT + 1/2 * Del_T * Qbar * convert_eps(-theta(i)) * Alpha * (h(i,1)^2 - h(i+1,1)^2);
    MH = MH + 1/2 * Del_C * Qbar * convert_eps(-theta(i)) * Beta  * (h(i,1)^2 - h(i+1,1)^2);
    
    k=k+1;
end

Nbar = N + NT + NH;
Mbar = M + MT + MH;
ABD = [A B;B D];
abd = inv(ABD);
a = abd(1:3,1:3);
b = abd(1:3,4:6);
bT= abd(4:6,1:3);
d = abd(4:6,4:6);

eps0_k = abd * [Nbar ; Mbar];
eps0 = eps0_k(1:3);
k    = eps0_k(4:6);
Alpha_bar = 1/Del_T * ( a * NT + b * MT);
Beta_bar  = 1/Del_C * ( a * NH + b * MH);


%set up distance matrix
z(1,1)=sum(th)/2;
for i=2:size(th,1)
    z(i-1,1)=z(i-1);
    z(i-1,2)=z(i-1)-th(i-1);
    z(i,1) = z(i-1,2);
end
z(size(th,1),2) = -sum(th)/2;




%Calculate strain and stress for each ply
j=size(th,1);
for i=1:size(th,1)
    Alpha_p = convert_eps(-theta(j)) * [Alpha_1(j);Alpha_2(j);0];
    Beta_p =  convert_eps(-theta(j)) * [Beta_1(j);Beta_2(j);0];
    z1 = (z(j,1)+z(j,2))/2;
    
    eval(['eps' int2str(i) '_upper = eps0 + z(' int2str(j) ',1)*k - Alpha_p * Del_T - Beta_p * Del_C;']);
    eval(['eps' int2str(i) '_lower = eps0 + z(' int2str(j) ',2)*k - Alpha_p * Del_T - Beta_p * Del_C;']);
    eval(['eps' int2str(i) '_mid = eps0 + z1*k -Alpha_p * Del_T - Beta_p * Del_C;']);
    
    eval(['sigxy' int2str(i) '_upper = Qbar' int2str(i) ' * eps' int2str(i) '_upper;']);
    eval(['sigxy' int2str(i) '_lower = Qbar' int2str(i) ' * eps' int2str(i) '_lower;']);
    eval(['sigxy' int2str(i) '_mid = Qbar' int2str(i) ' * eps' int2str(i) '_mid;']);
    
    
    j=j-1;
end

%Calculate Values for Tsai-Wu
%used Mises-Hencky Criterion
H1 = 1/F1t-1/F1c;
H11 = 1/(F1t*F1c);
H2 = 1/F2t-1/F2c;
H22 = 1/(F2t*F2c);
H6 = 0;
H66 = 1/(Fs^2);
H12 = (-1/2)*sqrt(1/(F1t*F1c*F2c*F2t));

%Calculate Tsai-Hill Values
j=size(th,1);
for i=1:size(th,1)
    
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %MID PLANE%
    eval(['sig' int2str(i) '_mid = convert_sig(theta(' int2str(j) ')) * sigxy' int2str(i) '_mid;']);
    eval(['sig = sig' int2str(i) '_mid']);
    
    if sig(1) > 0
        F1 = F1t;
    else
        F1 = F1c;
    end
    
    if sig(2) > 0
        F2 = F2t;
    else
        F2 = F2c;
    end
    
    
    %Max Stress/Strain Strength Ratio
    StrRat(i,:) = [F1/sig(1), F2/sig(2), Fs/sig(3)];
    
    v111 = F1/sig(1)
    v222 = F2/sig(2)
    v333 = Fs/sig(3)
    
    %Tsai-Wu Strength Ratio
    
    A = (sig(1)^2)/(F1t*F1c)+(sig(2)^2)/(F2t*F2c)+2*H12*sig(1)*sig(2)+(sig(3)^2)/(Fs^2);
    B = sig(1)*(1/F1t-1/F1c)+sig(2)*(1/F2t-1/F2c);
    C = -1;
    TsaiWuStrength = (-B+sqrt(B^2-4*A*C))/(2*A);
    
    if TsaiWuStrength < 0
        TsaiWuStrength = (-B-sqrt(B^2-4*A*C))/(2*A)
        StrRatTW(i) = TsaiWuStrength;
    else
        StrRatTW(i) = TsaiWuStrength;
    end
    
    
    %Tsai-Hill Strength Ratio
    
    A = (sig(1)/F1t)^2-(sig(1)*sig(2)/F1t^2)+(sig(2)/F2t)^2+(sig(3)/Fs)^2;
    B = 0;
    C = -1;
    TsaiHillStrength = (-B+sqrt(B^2-4*A*C))/(2*A);
    
    if TsaiHillStrength < 0
        TsaiHillStrength = (-B-sqrt(B^2-4*A*C))/(2*A);
        StrRatTH(i) = TsaiHillStrength;
    else
        StrRatTH(i) = TsaiHillStrength;
    end
    
    
    
    %Tsai-Wu Mid Test
    
    Value = H1*sig(1)+H2*sig(2)+H6*sig(3)+H11*(sig(1)^2)+H22*(sig(2)^2)+H66*(sig(3)^2)+2*H12*sig(1)*sig(2);
    if Value < 1 || Value == 0
        MWRes(i,:) = 'Passed';
    else
        MWRes(i,:) = 'Failed';
    end
    
    %Tsai-Hill Mid Test
    eval(['TH' int2str(i) '_mid = (sig(1)/F1)^2 + (sig(2)/F2)^2 - (sig(1)*sig(2)/F1^2) + (sig(3)/Fs)^2';]);
    if eval(['TH' int2str(i) '_mid']) >= 1
        THRes(i,:) = 'Failed';
    else
        THRes(i,:) = 'Passed';
    end   
    
    %Max Stress Mid Test
    if abs(sig(1))>F1 | abs(sig(2))>F2 | abs(sig(3))>Fs
        MSRes(i,:) = 'Failed';
    else
        MSRes(i,:) = 'Passed';
    end
    
    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%
    %UPPER PLANE%
    eval(['sig' int2str(i) '_upper  = convert_sig(theta(' int2str(j) ')) * sigxy' int2str(i) '_upper;']);
    eval(['sig = sig' int2str(i) '_upper']);
    
    
    %Tsai-Wu Upper
    
    Value = H1*sig(1)+H2*sig(2)+H6*sig(3)+H11*(sig(1)^2)+H22*(sig(2)^2)+H66*(sig(3)^2)+2*H12*sig(1)*sig(2);
    if Value < 1 || Value == 0
        MWRes(i,:) = 'Passed';
    else
        MWRes(i,:) = 'Failed';
    end
    
    
    
    if sig(1) > 0
        F1 = F1t;
    else
        F1 = F1c;
    end
    
    if sig(2) > 0
        F2 = F2t;
    else
        F2 = F2c;
    end
    
    
    %Max Stress Test
    if abs(sig(1))>F1 | abs(sig(2))>F2 | abs(sig(3))>Fs
        MSRes(i,:) = 'Failed';
    else
        MSRes(i,:) = 'Passed';
    end
    
   
    %Tsai-Hill Test
    eval(['TH' int2str(i) '_upper = (sig(1)/F1)^2 + (sig(2)/F2)^2 - (sig(1)*sig(2)/F1^2) + (sig(3)/Fs)^2';]);
    if eval(['TH' int2str(i) '_upper']) >= 1
        THRes(i,:) = 'Failed';
    else
        THRes(i,:) = 'Passed';
    end    
    
    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%
    %LOWER PLANE%
    eval(['sig' int2str(i) '_lower  = convert_sig(theta(' int2str(j) ')) * sigxy' int2str(i) '_lower;']);
    eval(['sig = sig' int2str(i) '_lower;']);
    if sig(1) > 0
        F1 = F1t;
    else
        F1 = F1c;
    end
    
    if sig(2) > 0
        F2 = F2t;
    else
        F2 = F2c;
    end

    %Max Stress Test
    if abs(sig(1))>F1 | abs(sig(2))>F2 | abs(sig(3))>Fs
        MSRes(i,:) = 'Failed';
    else
        MSRes(i,:) = 'Passed';
    end    
    
    %Tsai-Hill Test
    eval(['TH' int2str(i) '_lower = (sig(1)/F1)^2 + (sig(2)/F2)^2 - (sig(1)*sig(2)/F1^2) + (sig(3)/Fs)^2;']);
    if eval(['TH' int2str(i) '_upper']) >= 1
        THRes(i,:) = 'Failed';
    else
        THRes(i,:) = 'Passed';
    end
    
    %Tsai-Wu Lower Test
    Value = H1*sig(1)+H2*sig(2)+H6*sig(3)+H11*sig(1)^2+H22*sig(2)^2+H66*sig(3)^2+2*H12*sig(1)*sig(2);
    if Value < 1 || Value == 0
        MWRes(i,:) = 'Passed';
    else
        MWRes(i,:) = 'Failed';
    end
    
    %%%%%%%%%%%%%%%%STRAIN%%%%%%%%
    eval(['epss' int2str(i) '_mid  = convert_eps(theta(' int2str(j) ')) * eps' int2str(i) '_mid;']);
    eval(['epss = epss' int2str(i) '_mid']);
    
    if epss(1) > 0
        S1 = epsMax1t;
    else
        S1 = epsMax1c;
    end
    
    if epss(2) > 0
        S2 = epsMax2t;
    else
        S2 = epsMax2c;
    end
    
    
    %Max Strain Test
    if abs(epss(1))>S1 | abs(epss(2))>S2 | abs(epss(3))>epsMaxs
        SRes(i,:) = 'Failed';
    else
        SRes(i,:) = 'Passed';
    end
    
    eval(['epss' int2str(i) '_upper  = convert_eps(theta(' int2str(j) ')) * eps' int2str(i) '_upper;']);
    eval(['epss = epss' int2str(i) '_upper']);
    
    if epss(1) > 0
        S1 = epsMax1t;
    else
        S1 = epsMax1c;
    end
    
    if epss(2) > 0
        S2 = epsMax2t;
    else
        S2 = epsMax2c;
    end
    
    
    %Max Strain Test
    if abs(epss(1))>S1 | abs(epss(2))>S2 | abs(epss(3))>epsMaxs
        SRes(i,:) = 'Failed';
    else
        SRes(i,:) = 'Passed';
    end
    
    %Max Strain Theory
    eval(['epss' int2str(i) '_lower  = convert_eps(theta(' int2str(j) ')) * eps' int2str(i) '_lower;']);
    eval(['epss = epss' int2str(i) '_upper']);
    
    if epss(1) > 0
        S1 = epsMax1t;
    else
        S1 = epsMax1c;
    end
    
    if epss(2) > 0
        S2 = epsMax2t;
    else
        S2 = epsMax2c;
    end
    
    %Max Strain Test
    if abs(epss(1))>S1 | abs(epss(2))>S2 | abs(epss(3))>epsMaxs
        SRes(i,:) = 'Failed';
    else
        SRes(i,:) = 'Passed';
    end
    
    j=j-1;
end

%Output Results to Output.txt
OUT = fopen(Out, 'w');
fprintf(OUT, '%s\n',datestr(now));
fprintf(OUT, 'Results for [ %s ] Laminate:\n', num2str(theta'));

fprintf(OUT, '\nApproach Chosen:');
if Approach==0
    fprintf(OUT, '\nStrength of Material Approach');
elseif Approach==1
    fprintf(OUT, '\nElasticity Approach');
else
    fprintf(OUT, '\nSemi-empirical Approach\n\n');
end

fprintf(OUT, '\n');

%MAT1 = [eps0,k,Alpha_bar,Beta_bar];
%fprintf(OUT, '\nep0              k               Alpha_bar       Beta_bar\n');
%fprintf(OUT, '%+7.6e \t %+7.6e \t %+7.6e \t %+7.6e\n' ,MAT1');

fprintf(OUT, '\nEffective Moduli of Laminate:\n');
Ex = 1/a(1,1)/sum(th);
Ey = 1/a(2,2)/sum(th);
Gxy = 1/a(3,3)/sum(th);
vxy = -a(2,1)/a(1,1);
fprintf(OUT, 'Ex(psi)         Ey(psi)         Gxy(psi)        vxy\n');
fprintf(OUT, '%+7.6e \t%+7.6e \t%+7.6e \t%+7.6e\n' ,Ex,Ey,Gxy,vxy);

MAT1 = [Nbar,N,NT,NH];
fprintf(OUT,'\nNbar(lb/in)      N(lb/in)        NT(lb/in)       NH(lb/in)\n');
fprintf(OUT,'%+7.6e \t %+7.6e \t %+7.6e \t %+7.6e\n',MAT1');


fprintf(OUT, '\nFailure Test Results:\n');
fprintf(OUT, '\nPly \t Tsai-Hill \t     Max-Stress  \t Max-Strain  \t Tsai-Wu\n');
for i=size(theta,1):-1:1
    fprintf(OUT, ' %d \t %s \t\t %s \t\t %s \t\t %s\n ', theta(size(theta,1)+1-i),THRes(i,:),MSRes(i,:),SRes(i,:),MWRes(i,:));
end



fprintf(OUT, '\nStrength Ratios:\n\n');
fprintf(OUT, 'Ply  Max Stress/Strain   Tsai-Wu             Tsai-Hill\n');
%for i=size(theta,1):-1:1
%    min = intmax('uint64');
%    
%    for j=1:size(StrRat,1)
%        temp = StrRat(i,j);
%        if abs(temp) < min
%            min = abs(StrRat(i,j));
%            Strength = min;
%        end
%    end
%    
%    fprintf(OUT, '%d \t %+7.6e \t\t %+7.6e \t\t %+7.6e \t\t %+7.6e\n', theta(size(theta,1)+1-i), abs(Strength), StrRatTW(i), StrRatTH(i));
%    fprintf(OUT,'\n');
%end
    
    
fprintf(OUT,'\nPly  Upper strain(in/in) Mid sttrain(in/in)  Lower strain(in/in) Upper stress(psi)   Mid stress(psi)     Lower stress(psi)\n\n');
for i=size(th,1):-1:1
    MAT3 = [abs(i-size(th,1)*ones(3,1)-1),eval(['eps' int2str(i) '_upper']), eval(['eps' int2str(i) '_mid']), eval(['eps' int2str(i) '_lower']),eval(['sigxy' int2str(i) '_upper']), eval(['sigxy' int2str(i) '_mid']), eval(['sigxy' int2str(i) '_lower'])];
    fprintf(OUT,'%d \t %+7.6e    \t %+7.6e   \t %+7.6e   \t %+7.6e   \t %+7.6e   \t %+7.6e\n',MAT3');
    fprintf(OUT,'\n');
    MAT3g(3*i-2:3*i,:)=MAT3;
end

fclose(OUT);


%Internal Functions
function convert_sig=convert_sig(angle)

th=angle*pi/180;
m=cos(th);
n=sin(th);

convert_sig=[m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
    
function convert_eps=convert_eps(angle)

th=angle*pi/180;
m=cos(th);
n=sin(th);

convert_eps=[m^2 n^2 m*n; n^2 m^2 -m*n; -2*m*n 2*m*n m^2-n^2];