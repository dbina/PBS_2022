function test_ode
%
%
%
clear all
close all

global Ktest data2fit t DAS fit chisq res trace2fit amp

% t = (0:1:100);
% y = -0.3 * exp( -t/5) + 0.5*exp( -t/20 );
% 
% inpars = [2, 20];
% trace2fit = y;
% options = optimset('MaxFunEvals', 1000, 'Display','iter');
% f = lsqnonlin(@expfit, inpars, [], [], options);
% 
% 
% plot(t,y, '.')
% hold on
% plot(t,fit, '-')
% 
% return
% for ii =1:10
%     
%     plot( 0, ii, 'sq', 'MarkerFaceColor', rainbow( ii/10, 0, 1 ) )
%     hold on
%     
% end
% 
% return

% file to analyze:
fid = fopen( 'C:\DataVB\Progs\PBS_sims_2022\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid);   

% get kinetics
[fn, fp] = uigetfile( [fp0, '*.*' ],' load kinetics at selected times' )
D = csvread( [fp, fn]  );
% D = load(  [fp, fn]  );
t = D(2:end,1);
Y = D(2:end,2:end-1);
index = D(1,2:end-1);
b = 1:size( Y, 2);
classes = unique( index);

size( Y )

% form kinetics by subunits, not pigment type
% assumes: rod 1-6 (54 bilins each), then core (72 bilins), 4 carotenoids = 400 pigments
% 
% for n = 1:2
%     min( (n-1)*400+1:(n-1)*400+ 400 )
%     max( (n-1)*400+1:(n-1)*400+ 400 )    
% end
% 
% return
% size( Y )
% Y0 = Y;
% for n = 1:2
%     
%     Y = Y0(:, (n-1)*400+1:(n-1)*400+ 400 );
%     
%     y0 = t * 0; 
%     rodsize = 54;
% 
%     rod = [];
%     core = [];
%     cars = [];
%     for ii = 1:6
% 
%         rod(:,ii) = y0;
%         for jj = 1:rodsize
%         
%             rod(:,ii) = rod(:,ii) + Y( :, (ii-1)*rodsize+jj   ); 
%         
%         end       
%     
%     end 
%     core = sum( Y(:,6*rodsize+1:6*rodsize+72), 2);
%     cars = sum( Y(:,397:400), 2);
%     
%     
%     if n == 1
%         plot( t, rod, 'k-' )
%         hold on
%         plot( t, core, 'b-' )
%         plot( t, cars, 'r-' )    
%     else
%         plot( t, rod,  'k:' )
%         hold on
%         plot( t, core, 'b:' )
%         plot( t, cars, 'r:' )
%     end
%     
%     out = [ t, rod, core,cars];
%     save( ['C:\DataVB\data2analyze\PBS_OCP_complete\subunit_kins', num2str( n),'.txt'], 'out', '-ascii' )
%     out = [];    
% end
% return

% for ii = 1:length( t )
%     
%     bar( b, Y(ii,:) )
%     pause(0.001)
%     
% end

for ii = 1:length( classes)
    
    C(:,ii) = sum( Y(:, index == classes(ii)), 2 );
    
    outC(:,ii) = [classes(ii); C(:,ii)];
    
    leg{ii, 1} = num2str( classes(ii) );
    
end
legU = char( leg  );

outC = [ outC, sum( outC, 2 )  ];
norm_c = outC( 2, end );
outC(2:end, :) = outC(2:end, :) / norm_c;

figure(2)
subplot(1,2,1)
p1 = semilogx( t, C, '-sq' );
legend(p1, legU, 0);
hold on
semilogx( t, sum( C, 2) , '-+' )

subplot(1,2,2)
semilogy( t, C, '-sq' )
hold on
semilogy( t, sum( C, 2) , '-+' )
pause(eps)

outC = [ [999;t], outC];
% save( 'Z:\2021_09\pbs_sims\kins_of_groups_up_up_fix.txt', 'outC', '-ascii' )

% return
% get spectra
% [fnamea, pathnamea] = uigetfile( [ fp0,'\*.*'], 'absorption spectra' );
% [fname, pathname] = uigetfile( [ fp0,'\*.*'], 'emission spectra' );
% pathnamea = 'Z:\2021_06\data2analyze\PBS_OCP_complete\';
pathnamea = 'C:\DataVB\Progs\PBS_sims_2022\';
pathname = pathnamea;
fnamea = 'abs_spectra__exp_alt_current_show.csv';
fname = 'em_spectra__exp_alt_current_show.csv';

a = read_data( [pathnamea, fnamea ] );
f = read_data( [pathname, fname ] );

x = a(1:1:end,1);
a = a(1:1:end,2:end);
f = f(1:1:end,2:end);

% figure
% plot( x, f,'-' )
% 
% return
% 
% psii = exp( -0.5*(( x-1e7/690 )/ 200).^2); % trap
% f = [f, psii];
% 

scales = [1, 16^2;
          2, 16^2;
          3, 13^2; 
          4, 13^2;
          5, 13^2;
          6, 16^2;
          7, 16^2;
          8, 16^2;
          9, 1];
scales(:,2) = scales(:,2) / 16^2;

sindex = 1:size( f,2 );

S = zeros( size(f,1), length(t) );
uu = 1;
for ii = 1:length( classes )
    
    disp( [ 'kinetics for class ', num2str( classes(ii) ) ] )
    a = C(:, ii );
    b = f( :, sindex == classes(ii) );
    S = S + b * a';
%     S = S + b * a' * scales( scales( :, 1 ) == classes(ii), 2 );
 
end

% size(C)
% size(classes)
% disp( classes )
% size( f )
% f = f(:, 1:end-1);
% a = a(:, 1:end-1);
% size( f )
% % S = f * C';
% 
% size(x)
% size(S)

for ii = 1:length(t)
    S(:, ii) = S(:,ii).*(x.^2) * 1e-6;
end
m = max( max ( S));

figure(3)
plot( x, S, '-')
hold on

H = x.*0;
for ii = 1:size( C,2 )
    
    U(ii) = trapz( t, C(:,ii) );
    v = f(:,ii).*(x.^2) * 1e-6;
    H = H + U(ii) * v;
    
end

plot( x, H/max(H) * m, 'sq-')

z = [615, 640, 690];
for ii = 1:length(z)
    
    [k,n] = min( (z(ii) - 1e7./x).^2 )
    disp(1e7./x(n)) 
    kin(:,ii) = S( n, : )'
end

figure(4)
plot( t, kin, 'sq-')
pause(eps)

% analyze the 3D data
% fitting:
t = t';
Snorm = S / max(max(S));
inpars = [0.1, 5, 30, 200 ]; % no OCP
% inpars = [0.1, 5, 15, 200, 1500 ]; % no OCP+trap
% inpars = [0.5, 1, 6, 30, 100, 200]; % OCP
% inpars = [1, 10, 50, 120]; 
% inpars = [0.1, 10,  40, 1300]; % OCP
data2fit = Snorm + randn( size(Snorm,1), size(Snorm,2) )*0.0;
options = optimset('MaxFunEvals', 1000, 'Display','iter');
%[f,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@model2fit, inpars, [], [], options);
f = lsqnonlin(@model2fit, inpars, [], [], options)

disp( fn )

SS = x*0;
% f = [f, 1350];
% for ii = 1:length(f)
%     
%     disp( f(ii) )
%     SS = SS + DAS(:,ii)/f(ii);
%     
% end
% f = f(1:end-1);

figure(5)
plot( x, DAS, '-sq')
hold on
plot( x, SS, 'k-')
pause(eps)
% [1e7./x, DAS]

picker = 1:5:size( data2fit,1 );
figure(6)
subplot(1,4,1)
plot(t( t>-3 & t < 5)' , data2fit(picker,  t>-3 & t < 5)', '.' )
hold on
plot(t( t>-3 & t < 5)' , fit(picker,  t>-3 & t < 5)', '-' )
hold on

subplot(1,4,2)
plot(t( t>-10 & t < 500)' , data2fit(picker,  t>-10 & t < 500)', '.' )
hold on
plot(t( t>-10 & t < 500)' , fit(picker, t>-10 & t < 500)', '-' )

plot(t( t>-10 & t < 500)' , -0.2+fit(picker, t>-10 & t < 500)'-data2fit(picker,  t>-10 & t < 500)', '.' )

pause(eps)


subplot(1,4,3)
semilogy(t( t>-10 & t < 500)' , data2fit(picker,  t>-10 & t < 500)', '.' )
hold on
semilogy(t( t>-10 & t < 500)' , fit(picker, t>-10 & t < 500)', '-' )
pause(eps)

subplot(1,4,4)
plot(t' , Snorm', '.' )
hold on
plot(t' , fit', '-' )
pause(eps)

size(C)

% for ii = 1:size(C,2)
%     
%     inpars = [0.5, 1,100];
%     trace2fit = C(:,ii)';
%     options = optimset('MaxFunEvals', 1000, 'Display','iter');
%     f = lsqnonlin(@expfit, inpars, [], [], options);
%   
%     subplot(1,9,ii)
%     semilogx( t, trace2fit,'.')
%     hold on
%     semilogx( t, fit,'-')
%     pause(eps)
% end

% convert DAS as acquired in [1/cm] x-scale to [nm] scale
dasout = [ [999;x], [[f, 999 ]; DAS] ];
% dasout = [ [999;1e7./x], [f; DAS] ];

namek = fn(1:end-4);
[dasfn, dasfp] = uiputfile( [fp0, 'sim_das_ctrl_3A.txt'] )
if dasfn ~= 0
    save( [dasfp,dasfn], 'dasout', '-ascii')    
    disp( 'DAS saved')
%     out = [ [999;t'], outC ];
%     save( [dasfp,'\', namek,'_classes.txt'], 'out', '-ascii' )

%     out = [ t', fit' ];
%     save( [dasfp,'\sim_kins_of_classes_test_fit_3A.txt'], 'out', '-ascii' )

end


disp( 'end OK' )

%====================================================================
function dy = model( t, y )
% basic model
% 
%
global irf_t0 irf_fwhm f a Ktest

% pulsed excitation
%P = pulse(t, irf_t0, irf_fwhm, f);
% 
% Ktest(4,4) = -P;
% Ktest(1,4) = P;
% dy = Ktest * y;

dy = Ktest * y;

%--------------------------------------------------------------------
function data = read_data( filename )
%
%
%
A = textread( filename, '%s');

for ii =1:size(A,1)-1
   if ii > 3
       
       data(ii-3,:) = str2num( char(A{ii}) );
   else
       disp( A{ii} )
   end
end

%----------------------------------------------------------------------
function v = rainbow(x, xmin, xmax)
%
%
%
% figure(1)

% ii = (1:90)'
%    r = exp( -0.5 * ((ii-15)./13).^2); 
%    g = exp( -0.5 * ((ii-45)./13).^2);
%    b = exp( -0.5 * ((ii-75)./13).^2);
% 
%    plot( ii, [r,g,b], '-')
%    
% figure(2)

ii = (x - xmin) / (xmax - xmin) * 90;
ii = max( [ ii, 0 ] );
    
%for ii = 1:90
    
   r = exp( -0.5 * ((ii-15)/33)^2) ;
   g = exp( -0.5 * ((ii-45)/33)^2);
   b = exp( -0.5 * ((ii-75)/33)^2);
   
   s = r+g+b;
   r = r / s;
   g = g / s;
   b = b / s;   
   
   v = [b,g,r];
   v = v / max(v);
   
  % plot( [0,1], [ii,ii] ,'-', 'Color' , v, 'LineWidth', 3 )
  % hold on
   %end
   
%---------------------------------------------------------------------
function f = model2fit( pars )
%
%
%
global data2fit t DAS fit chisq res

k = [1./pars, 1/1350 ]; % 1/1350
% shift = pars(end-1);
% sigma = pars(end) / 2.355; % from fwhm to std

for ii = 1:length(k) % adds a component in each cycle, i.e. exponential decay (x) irf
%    kins(ii, :) = 0.5 * exp(-k(ii)*t).*(1+erf((t-shift-k(ii)*sigma^2)/sqrt(2)/sigma))*exp(k(ii)*(shift+k(ii)*sigma^2/2));
     kins(ii, :) = exp(-k(ii)*t);
end
% kins(ii+1, :) = ones(1, length(t)); % offset, use for fluorescence
% kins(ii+1,:) = 0.5*exp( -0.5 * ((shift - t)/sigma).^2); % coherent artifact

a = (kins'\data2fit')'; % compute amplitudes
z = a * kins;
% s = data2fit;
% s(s <= 0) = 1;
f = (z - data2fit);% ./ sqrt(s);
fit = z;
DAS = a;
chisq = sum( sum( f.^2)); 
res = f;

% k = [1./pars(1:end-2)];
% % shift = pars(end-1);
% % sigma = pars(end) / 2.355; % from fwhm to std
% 
% for ii = 1:length(k) % adds a component in each cycle, i.e. exponential decay (x) irf
%     kins(ii, :) = 0.5 * exp(-k(ii)*t).*(1+erf((t-shift-k(ii)*sigma^2)/sqrt(2)/sigma))*exp(k(ii)*(shift+k(ii)*sigma^2/2));
% end
% % kins(ii+1, :) = ones(1, length(t)); % offset, use for fluorescence
% kins(ii+1,:) = 0.5*exp( -0.5 * ((shift - t)/sigma).^2); % coherent artifact
% 
% a = (kins'\data2fit')'; % compute amplitudes
% z = a * kins;
% % s = data2fit;
% % s(s <= 0) = 1;
% f = (z - data2fit);% ./ sqrt(s);
% fit = z;
% DAS = a;
% chisq = sum( sum( f.^2)); 
% res = f;

%-----------------------------------------------
function f = expfit( pars )
%
%
%
%global data2fit fit chisq resids t irf x DAS irfpars rchisq offset fixpars freepars
global trace2fit t amp fit
% t is horizontal
t = t(:)';
tau = (1./pars);
TAU = tau' * ones( 1, length(t) ) ; 
T = ones( length(tau), 1 ) * t;  
F = exp( - T ./ TAU ); % generates a matrix in which each row is one decay component
D = ones( length(tau), 1 ) * trace2fit;
a = [F]' \ D';% solves matrix equation to obtain amplitudes of the components for a single data trace  
fit = a' * [F]; % generates the fit of the given trace, step 1; this has n = length(tau) identical rows
% size(fit)
% dsize(trace
f = fit(1,:) - trace2fit;
fit = fit(1,:);

