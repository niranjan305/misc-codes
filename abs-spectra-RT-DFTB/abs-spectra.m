% THIS PROGRAM PLOTS THE ABSORPTION SPECTRUM OF ANY SYSTEM
% GIVEN THE TIME-DEPENDENT DIPOLE MOMENTS OBTAINED BY
% RUNNING THE REAL-TIME DFTB CALCULATIONS

% First open the files and extract all the information
% related to the dipole moments. 

fileIDx = fopen('muxx.dat','r');
fileIDy = fopen('muyy.dat','r');
fileIDz = fopen('muzz.dat','r');

formatSpecx = '%f %f';
formatSpecy = '%f %f';
formatSpecz = '%f %f';

sizeX = [2 Inf];
sizeY = [2 Inf];
sizeZ = [2 Inf];

X = fscanf(fileIDx,formatSpecx,sizeX);
Y = fscanf(fileIDy,formatSpecy,sizeY);
Z = fscanf(fileIDz,formatSpecz,sizeZ);

% Close the files

fclose(fileIDx);
fclose(fileIDy);
fclose(fileIDz);

% Normalize the three dipole moments around zero and then
% taking the average of the three
xini = X(2,1);
yini = Y(2,1);
zini = Z(2,1);
x = X(2,:)-xini;
y = Y(2,:)-yini;
z = Z(2,:)-zini;
avg_mu = (x+y+z)./3;

% Extract time scale from any one of the files
t=X(1,:);

% Time step in a.u. (This must match the value used while
% running the RT-DFTB calculations
dt = 0.2;

% Since we need to take a Fourier Transform of the time-varying dipole
% moments, we need to apply a damping function to make the dipole
% moment go to zero after some time. pdamp variable is that damping.
% Increasing the damping increases the width of the peak.
% Standard has been decided as 0.008 for pdamp
pdamp = 0.008;

% Find the time period for which dynamics has been done
cnt = length(t);
damp = zeros(1,cnt);

% Applying damping function
for i=1:cnt
    damp(i)=exp(-i*dt*pdamp);
end
avg_mu = avg_mu.*damp;

% Increase time scale to 400fs and add zeroes to the expanded array
expanded = zeros(1,62016);
avg_mu = cat(2,avg_mu,expanded);

% To ensure that 
N=82688;
% Compute the Fourier transform of the average dipole moment
MU = fft(avg_mu,N);
% Extract the imaginary part of the dipole moment
mu_mag = imag(MU);
% Absolute value of the imaginary part of the dipole moment
mu_mag = -(mu_mag);
% Calculate energy scale
E = (0:N-1)/(N*0.2)*2*pi*27.211396641344194;
lamda = 4.13*10^-15*2.99*10^8./E;
lamda = lamda .* 10^9;
% Absorbance is energy * imag part of (fft of dipole)
mu_mag=E.*mu_mag;
% Plot absorbance Vs energy but only first 5000 points
np = 5000;
Ep = E(1:np);
mu_mag = mu_mag(1:np);

% If u want to plot absorption spectrum w.r.t to wavelength,
% replace Ep with lambda in the plot line
lamda = lamda(1:np);

% Plotting the absorption spectrum
figure;
hold on;
box on;
plot(Ep,mu_mag,'r','Linewidth',1.5);
% Set energy axis
xlim([0 10]);
xlabel('Energy (eV)','Fontsize',18,'Fontweight','normal','Fontsize',24);
ylabel('Absorbance (arb. units)','Fontsize',18,'Fontweight','normal','Fontsize',24);
title('Absorption Spectrum','Fontsize',11,'Fontweight','normal','Fontsize',18);
set(gca,'XMinorTick','on','YMinorTick','on','ticklength',2*get(gca,'ticklength'),'LineWidth',1.5,'FontSize',24);
print('abs-spectra','-dpng','-r600');
hold off;
