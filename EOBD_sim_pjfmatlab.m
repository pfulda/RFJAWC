% EOBDsim

% load in the phasemap
phase_m = load('interpphasemap.txt');

mm=1e3;
% set up the grid axes
X = linspace(-0.002, 0.002, 401);
Y = X;

[XX, YY] = meshgrid(X, Y);

lambda=1064e-9;
% set an initial beam parameter. Here we do the simple version of the EOBD,
% so the beam is at the waist
gp=FT_init_gauss_param(lambda,1.7,'w0',250e-6,'z',0);

% set a container for mode coefficients
gc=FT_init_gauss_coefficients('HG',6);

% create the initial beam
initbeam=FT_HG_field(gp,0,0,X,Y,[0,0,0]);

% stuff just for plotting the phase map with the beam on top
phase_max=max(max(phase_m));
phase_mnorm=(phase_m/phase_max+0.125)*4;
phase_mrgb = cat(3,phase_mnorm,phase_mnorm,phase_mnorm);
figure()
image(phase_mrgb)
%set(gca,'clim',[min(min(phase_m)) max(max(phase_m))]);
hold all
contour(abs(initbeam).^2)
set(gca,'clim',[min(min(abs(initbeam).^2)) max(max(abs(initbeam).^2))]);
%view(2)
axis square

% apply the phase map to the beam (exaggerate 100x for plot option)
mapscale=1;
%mapscale=100;
defbeam=initbeam.*exp(1i*phase_m*mapscale);

% decompose in the intial beam basis
defbeamcontent=FT_mode_content(gc,gp,defbeam,X,Y,[0,0,0]);

% print the results
FT_print_gauss_coefficients(defbeamcontent,1e-9,1);

% try propagating some distance
gpnew=FT_update_gauss_param(gp,gp.q+1,gp.nr);

% recreate the field from the propagated beam
propbeam=FT_mode_coefficients_to_field(defbeamcontent,gpnew,X,Y,[0,0,0],0);
propbeam_int=abs(propbeam).^2;
propbeam_int=propbeam_int/max(max(propbeam_int));

defbeam_int=abs(defbeam).^2;
defbeam_int=defbeam_int/max(max(defbeam_int));

if mapscale==100
    figure()
    subplot(211)
    plot(X*mm,defbeam_int(200,:))
    title('Deflected beam at deflection location and propagated 1m (100V)')
    hold on 
    plot(X*mm,propbeam_int(200,:))
    ylabel('Norm. Intensity [A.U.]')
    xlabel('Distance from central axis [mm]')
    grid on
    legend('initial','prop')
    subplot(212)
    plot(X*mm,angle(defbeam(200,:))*180/pi)
    hold on
    plot(X*mm,angle(propbeam(200,:))*180/pi)
    grid on
    ylabel('Phase [deg]')
    xlabel('Distance from central axis [mm]')
    legend('initial','prop')
    saveas(gcf,'EOBD_simple_test.pdf')
else 
end
%% now try splitting it up into multiple scatter matrices

%phase map per slice is the total map divided by the number of slices
phase_mslice=phase_m/11;

% some z steps for the slices
zs=linspace(-0.2,0.2,11);

figure()
for n=1:11
    
    % set new q parameter
    gp=FT_init_gauss_param(lambda,1.7,'w0',250e-6,'z',zs(n));
    % if first one, generate initial field
    if n==1
        beam=FT_HG_field(gp,0,0,X,Y,[0,0,0]);
    % otherwise, build field from previous set of mode coefficients and new
    % q parameter
    else
        beam=FT_mode_coefficients_to_field(defgc,gp,X,Y,[0,0,0],0);
    end
    % apply phasemap of this slice
    defbeam=beam.*exp(1i*phase_mslice);
    % calculate new mode coefficients
    defgc=FT_mode_content(gc,gp,defbeam,X,Y,[0,0,0]);
    % plot the phase front of this slice
    plot(X,angle(defbeam(50,:))*180/pi)
    hold on
end

FT_print_gauss_coefficients(defgc,0,1)
