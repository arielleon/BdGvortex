% BdG vortex:
%
% self consistency equation:
% Delta(rho) = g/(2*Pi*L)*Sum_{n,l,kz,j,j'}
% c_{n,l,j}*d_{n,l,J}*Phi_{l,j}(rho)*Phi_{l+1,J}(rho)
% 
%
% From this version, we calculate free energy.

profile on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical para
invkFa=-1;
Ec=100;
g=8*pi/(invkFa-2/pi*sqrt(Ec));
mu = 0.96;
xi = 1/sqrt(2*0.3496^2*abs(g));

m = 1/2;
hbar = 1;
a = -1;
kmu = sqrt(2*m*mu)/hbar;
% computational para
lmax = 50;   % # of angular momentum modes included
jmax = 150;   % # of Bessel function bases included 
R = 30.0;
nrho = 100;
drho = R/nrho;
Lz = 10; 
kzmax = 15; 
dkz = 2*pi/Lz;
nkzmax = floor(kzmax/dkz); % correspond to Ec = 100
% kzmax = dkz*nkzmax;
% Lz = 2*pi/dkz;
jump = 1;
iteration = 10;
alpha = 0.1;

Energy = zeros(lmax,jmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
% bessel function base Phijl(rho)
% case on: create PhiCell
if jump == 0;
alphamat = zeros(jmax, 2*lmax);
for l = (1:2*lmax+1)
    alphamat(:,l) = besselzero(abs(l-1-lmax), jmax,1);
end
%

%fname = ['alphamatLib.txt'];
%loadalphamat = dlmread(fname);
%alphamat = loadalphamat;

Phimat = zeros(jmax,nrho);
PhiCell = cell(1,2*lmax);
[PhiCell{:}] = deal(Phimat);
for l = 1:(2*lmax+1)
    for j = 1:jmax
        for r = 1:nrho
            Phimat(j,r) = sqrt(2)*besselj(l-lmax-1,alphamat(j,l)*r*drho/R)/(R*besselj(l-lmax,alphamat(j,l)));
        end
    end
   [PhiCell{l}] = Phimat;
end
% case off.

% initiallize Orderpara(rho) and construct OrderparaBig matrix
% case on: vortex delta
%deltanew = 0.3496*transpose([tanh(drho*(1:nrho))]);
% case off.
% case on: bulk superconductor
deltaini = transpose(drho*(1:nrho).*((((1:nrho)*drho).^2+xi^2).^(-0.5)));
deltaini = 0.3496*deltaini;
deltanew = deltaini;
% case off.
figure(8);
plot((1:nrho)*drho, deltanew);
Deltamat = diag(deltanew);
Rhomat = diag([(1:nrho)*drho]);

DeltaCell = cell(1,2*lmax);
[DeltaCell{:}] = deal(zeros(jmax));
for l = 1:(2*lmax)  % DeltaCell has only been given values up to (lmax-1)-th element. 
                    % lmax-th element is so far zeors. 
    [DeltaCell{l}] = (PhiCell{l}*Rhomat*Deltamat*transpose(PhiCell{l+1})).*drho;  % att: choose bulk superfluid here. 
    %[DeltaCell{l}] = 0.3496*eye(jmax);
end
    
% case on: include nkz
 % construct HCell
Hmat = zeros(jmax);
HCell = cell(2*lmax+1,nkzmax);
[HCell{:}] = deal(Hmat);


for nkz = 1:nkzmax
    for l = 1:(2*lmax+1)
        for j = 1:jmax
            % case: include kz
             Hmat(j,j) = (alphamat(j,l))^2/R^2 + ((nkz-1)*dkz)^2 - mu;           
            % case: neglect kz
            % Hmat(j,j) = alphamat(j,l)^2/R^2 - mu;
        end
        [HCell{l,nkz}] = Hmat;
    end
end
end
 % construct lth BdG equation
HamilCell = cell(2*lmax+1,nkzmax);
[HamilCell{:,:}] = deal(zeros(2*jmax));
ECell = cell(2*lmax+1,nkzmax);
[ECell{:,:}] = deal(zeros(j));

cdtemp = zeros(jmax);
cdCell = cell(jmax,2*lmax+1,nkzmax); % index n,l,kz, to store the matrix c_j*d_J
[cdCell{:}] = deal(zeros(jmax));
ddtemp = zeros(jmax);
ddCell = cell(jmax,2*lmax+1,nkzmax); % index n,l,kz, to store the matrix c_j*d_J
[ddCell{:}] = deal(zeros(jmax));
% case off.

EnergyCell = cell(2*lmax+1,nkzmax);
[EnergyCell{:}] = deal(zeros(jmax,1));
EnergynegCell = cell(2*lmax+1,nkzmax);
[EnergynegCell{:}] = deal(zeros(jmax,1));


vector = zeros(2*jmax);
value = zeros(2*jmax);
c = zeros(jmax);
d = zeros(jmax);
egtemp = zeros(nrho/2,1);
egneg = zeros(nrho/2,1);
egcutoff = zeros(nrho/2,1);
egnegcutoff = zeros(nrho/2,1);
cdcutoff = zeros(nrho/2,1);

 %case on: sum over all l
for i = 1:iteration % interation

    Deltamat = diag(deltanew);
    %DeltaCell{1} = (PhiCell{1}*Rhomat*Deltamat*transpose(PhiCell{1})).*drho;
    deltatemp = zeros(nrho,1);
    utemp = zeros(nrho,1);
    vtemp = zeros(nrho,1);
    uvtemp = zeros(nrho,1);
    feg1temp = zeros(nrho,1);
    feg1 = zeros(nrho,1);
    feg2coef = 0;
    feg2 = zeros(nrho,1);
    feg = zeros(nrho,1);
    
    for l = 1:2*lmax  % Sum_nkz % att: 2*lmax - 2 in the case of vortex phase
                                    %      2*lmax - 1 in the case of bulk
                                    %      superconductor.
        [DeltaCell{l}] = (PhiCell{l}*Rhomat*Deltamat*transpose(PhiCell{l+1})).*drho;
        for nkz = 1:nkzmax % Sum_l 
            [HamilCell{l,nkz}] = [HCell{l,nkz}, DeltaCell{l}; transpose(DeltaCell{l}), -HCell{l+1,nkz}];   % att: lmax-th of DeltaCell is zeros.
            [vector, value]=eig(HamilCell{l,nkz});
            c = vector(1:jmax, jmax+1:2*jmax); % att: normalization hasn't been added
            d = vector(jmax+1: 2*jmax, jmax+1:2*jmax);
            egtemp = diag(value);
            eg = egtemp(jmax+1:2*jmax);
            egcutoff = (abs(eg)<Ec);
            eg1 = eg;
            eg (abs(eg)> Ec)= NaN;
            [EnergyCell{l,nkz}] = eg;
            cdcutoff = diag(egcutoff);
            c = c*cdcutoff;
            d = d*cdcutoff;
    
    
            for n = 1:jmax  % Sum_n
                utemp = transpose(PhiCell{l})*c(:,n);
                vtemp = transpose(PhiCell{l+1})*d(:,n);
                uvtemp = utemp.*conj(vtemp);
                deltatemp =  deltatemp - uvtemp;
                
                feg1temp = abs(vtemp).^2*(-2*eg1(n));
                feg1 = feg1 + feg1temp;
                feg2coef = feg2coef + 1/(2*HCell{l,nkz}(n,n));
                % n seems to be relabeld when getting the eigensystems
                % Does it matter when doing the summation? 
      
            end       
        end
        
    end
    feg2 = (feg2coef*2/R^2-Lz*m/(2*hbar^2*a))* abs(deltanew).^2;
    deltatemp = 2*deltatemp*g/(2*pi*Lz);
    deltaold = deltanew;
    deltanew = (1-alpha)*deltatemp+alpha*deltaold;
    feg = feg1+feg2;
    

end
% case off.
t = toc/60



% Calculate free energy from eigenvectors.

figure(8);clf;
subplot(2,2,1)
hold on;
plot((1:nrho)*drho,deltaini,(1:nrho)*drho,deltaold, (1:nrho)*drho, deltatemp);
hold off;
axis([0 R -0.2 1.2])
xlabel('r k_F')
ylabel('\Delta / E_F')

cut = floor(1/2*nrho);
fegnormcoef = feg(cut);
fegnorm = fegnormcoef*ones(nrho,1);
fegvaltemp = Rhomat*feg;
fegval = sum(fegvaltemp(1:cut))*drho
fegnormvaltemp = Rhomat*fegnorm;
fegnormval = sum(fegnormvaltemp(1:cut))*drho;
fegvorval = (fegval - fegnormval)/(mu*kmu*Lz)

%figure(9);
%plot((1:nrho)*drho,densenew);
%figure(9); clf;
subplot(2,2,2)
hold on;
plot((1:nrho)*drho,feg/(mu*kmu*Lz));
plot((1:nrho)*drho,fegnorm/(mu*kmu*Lz));
hold off;
%axis([0 R -0.2 0.2])
xlabel('r k_F')
ylabel('F / (\mu k_\mu  Lz)')



%figure(10); clf;
%plot(1:jmax, flipud(EnergyCell{lmax+1,1}));
%hold on;
%for nkz = 2:nkzmax
%    plot(1:jmax, EnergyCell{41,nkz})
%end
%hold off;
%xlabel 'bessel base index j'; ylabel 'Energy / E_F';


%figure(11); clf;
subplot(2,2,3)
hold on;
for i = -lmax:lmax-1
plot(i, EnergyCell{lmax+1+i,1},'.','color','b','MarkerSize',12)
end
hold off
axis([-lmax lmax-1 -Ec-20 Ec+20])
xlabel('angular momentum number l')
ylabel('E / E_F')

%figure(12); clf;
subplot(2,2,4)
hold on;
for i = -lmax:lmax-1
      plot(i, EnergyCell{lmax+1+i,1},'.','color','b','MarkerSize',12)
end
hold off
line('XData', [-lmax lmax], 'YData', [0.3496 0.3496], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color','m');
axis([-lmax lmax-1 -1 1])
xlabel('angular momentum number l')
ylabel('E / E_F')

figure(9);clf;
hold on;
plot((1:nrho)*drho,feg2/(mu*kmu*Lz),'r');
plot((1:nrho)*drho,feg1/(mu*kmu*Lz),'b');
hold off;
xlabel('r k_F')
ylabel('F / (\mu k_\mu  Lz)')


t=toc/60
beep on;
beep;

%profile viewer
%p = profile('info');
%profsave(p,'profile_results')
    




    






    


       
