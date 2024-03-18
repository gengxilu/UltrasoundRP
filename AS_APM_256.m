
clear
%% transducer parameters
c=1.5; %km/s
Lh=12; %mm
f=4.5; %MHz
lamb=1.5/f; %mm
pixel=0.75/12; %mm
N_pixel=round(Lh./pixel); %mm
k=2*pi/lamb;
kx=2*pi/pixel;
z=15; % mm
% Li=14; % mm
% k_lim=pi*(Lh+Li)/lamb/sqrt((Lh+Li)^2/4+z^2);

%% image generation
raw_I=imread(['V.png']);
gray_I=double((raw_I));
I=imresize(gray_I, round([N_pixel N_pixel]));
I(I<0.5)=0;
I(I>0.5)=1;
figure, imagesc(I);

%% Angular Spectrum; zero padding
tic
F_I=fftshift(fft2(1j.*I, 2*N_pixel-1, 2*N_pixel-1));  % range in (-kx/2 : kx/2). Centered
% I2=ifft2(F_I);
% figure, imagesc(abs(F_I));

dkx=kx/(2*N_pixel-2);
PointNum=floor(k/dkx);

% create H, transfer matrix, Centered.
H=ones(2*N_pixel-1, 2*N_pixel-1);
H=([-(N_pixel-1):N_pixel-1].^2).*H;
H=H+H';     % unit: (lamb/2)^2
H_cut=H;
H_cut(H_cut>PointNum^2)=PointNum^2;
H_cut=sqrt(PointNum^2*ones(2*N_pixel-1, 2*N_pixel-1)-H_cut);
H_complex=exp(-1j*z*dkx*H_cut);    % image to hologram: -1; hologram to image: +1;

% create calibration matrix: CM. 
% We need calibration because the ifft was conducted from 0 to 2N-1; But it
% should be conducted from -(N-1) to N-1. There is a linear shift 
CM=[-(N_pixel-1):N_pixel-1].*ones(2*N_pixel-1, 2*N_pixel-1);
CM=CM+CM';
CM=exp(-1j*2*pi*N_pixel/(2*N_pixel-1).*CM);

F_holo=F_I.*H_complex;
holo_recon=ifft2(F_holo).*CM;
toc

figure,imagesc(abs(holo_recon(1:N_pixel,1:N_pixel)));
% figure,imagesc(abs(holo_recon(1:end,1:end)));
figure,imagesc(angle(holo_recon(1:N_pixel,1:N_pixel)));
% save('USC trojan hologram_14mm_8mm_3MHz.mat','holo_recon','N_pixel');

% F_I_recon=F_holo.*conj(H_complex);
% I_recon=ifft2(F_I_recon);
% figure,imagesc(abs(I_recon(1:N_pixel, 1:N_pixel)));

%% downsample
holo_recon_array = holo_recon;
for ii = 1:16
    for jj = 1:16
        holo_recon_array([1:3]+(ii-1)*3,[1:3]+(jj-1)*3) = holo_recon(3*ii-1,3*jj-1);
    end
end
        
figure,imagesc(abs(holo_recon_array(1:N_pixel,1:N_pixel)));
figure,imagesc(angle(holo_recon_array(1:N_pixel,1:N_pixel)));

%% direct integral
holo_recon_real=zeros(size(H));
holo_recon_real(1:N_pixel,1:N_pixel)=holo_recon_array(1:N_pixel,1:N_pixel);
z_lambda=z/pixel;  % unit: pixel
dR=sqrt(H+z_lambda^2); % unit: pixel
g=exp(1j*k*dR*pixel)./(dR.*pixel);           %  propagation function. -1j: forward propagation; +1j: reverse propagation.
% g=exp(1j*k*dR*pixel).*(z./dR)./(dR.*pixel).*(1./(2*pi*dR.*pixel)-1j*k); %  propagation function. -1j: forward propagation; +1j: reverse propagation.
recon_image=conv2(holo_recon_real,g);
figure, imagesc(abs(recon_image(N_pixel+1:2*N_pixel, N_pixel+1:2*N_pixel)));
colorbar

% save('S 16by16.mat','holo_recon_array')
% figuredata=abs(recon_image(N_pixel+1:2*N_pixel, N_pixel+1:2*N_pixel));
% saveas(gcf,'APM theory recon_14mm_8mm_3MHz.png');
% figure, imagesc(abs(recon_image));

% figure, imagesc(angle(recon_image));
