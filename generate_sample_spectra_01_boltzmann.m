function overlaps_01
% this function creates files of absorption and 
% emission spectra for computation of overlaps
%
% ouput file format:
% 
%
%
% (c) DB 2022-Jan-14

clear all
% close all

global shifttem
% file to analyze:
fid = fopen( 'C:\DataVB\Progs\PBS_sims_2022\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid);   

% define labels for spectroscopic species, NO spaces
splabels = {'APC_a84', 'APC_b84', 'CPC_a84', 'CPC_b84', 'CPC_b155', 'APCbottom_84', 'APCE_190', 'APCbottom_a84', 'car'}
% now asign to them numberical labels
numlabels = [1,2,3,4,5,6,7,8,9];

% comment in spectra names
comm = '_exp_alt_current_show';

% x = (14300: 1 : 18100)' ; % [1/cm]
x = (11500: 1 : 18100)' ; % [1/cm]
% absorption spectra

v(1) = 1e7/651; % a84 APC
v(2) = 1e7/651; % b84 APC
v(3) = 1e7/618; % a84 CPC
v(4) = 1e7/624; % b84 CPC
v(5) = 1e7/594; % b155 CPC
v(6) = 1e7/655; % a84 APCbottom
v(7) = 1e7/655; % APCE_190
v(8) = 1e7/655; % APCbottom_a84_2
v(9) = 13800; % car

shifttem = 2; % [nm]

% assign stokes shifts
dstokes = ones( length( v ),1 ) * 550;
% now, particulars
dstokes(9) = 20000;
dstokes(1:2) = 280;
width = ones( length( v ),1 ) * 630; % before, it was 400
width(9) = 3800;

for ii = 1:length( v )

    % abs spectrum
    a(:,ii) = gpeak( v(ii), width(ii), 1, x );

    % emission spectrum
    f(:,ii) = gpeak( v(ii)-dstokes(ii), width(ii), 1, x );

end
a(:,9) = carotenoidS1( v(9),x );
f(:,9) = carotenoidS1S0f( v(9),x );

sw = 1;
% % use experimental spectra
a(:,1) = exp_spectrum( x, 4, 'abs', sw );
f(:,1) = exp_spectrum( x, 4, 'flu', sw );

a(:,2) = exp_spectrum( x, 4, 'abs', sw );
f(:,2) = exp_spectrum( x, 4, 'flu', sw );

a(:,3) = exp_spectrum( x, 1, 'abs', sw);
f(:,3) = exp_spectrum( x, 1, 'flu', sw );

a(:,4) = exp_spectrum( x, 2, 'abs', sw );
f(:,4) = exp_spectrum( x, 2, 'flu', sw );

a(:,5) = exp_spectrum( x, 3, 'abs', sw );
f(:,5) = exp_spectrum( x, 3, 'flu', sw );

a(:,6) = exp_spectrum( x, 5, 'abs', sw );
f(:,6) = exp_spectrum( x, 5, 'flu', sw );

a(:,7) = exp_spectrum( x, 5, 'abs', sw );
f(:,7) = exp_spectrum( x, 5, 'flu', sw );

a(:,8) = exp_spectrum( x, 5, 'abs', sw );
f(:,8) = exp_spectrum( x, 5, 'flu', sw );

for ii = 1:8
    a(:,ii) = a(:,ii) / max(a(:,ii));
    f(:,ii) = f(:,ii) / max(f(:,ii));
end

% get energies
for ii = 1:8
    
    [u1, u2] = max( f(:,ii) );
    EF(ii,:) = [ii, x( u2 ) ];
    
    [u1, u2] = max( a(:,ii) );
    EA(ii,:) = [ii, x( u2 ) ];
        
end 

Emean = (EF + EA)/2; 
E = Emean;

showE = [ EA(:,2), EF(:,2), Emean(:,2) ];
disp( showE )
disp( 1./showE )

n2 = 4;
K = (72-n2)/n2 * exp( -(E(1,2)-E(6,2))/(0.695*293) ) 

% return

E(9,:) = [9,13800];
range = [550, 750];
q = x > 1e7/(max(range)) & x < 1e7/(min(range));
% now plot the spectra as real fluorescence and absorbance
% in wavelength axis
figure(1)
for ii = 1:length(v)
    
    offset = 0 + (ii-1)*1.1;
    % absorbance
    showa = a(:, ii).* x;  % bandshape to absorbance spectrum, band * freq
    showa = showa / max( showa) ;
    
    % absorbance
    showf = f(:, ii).* x.^3; % bandshape to fluorescence spectrum, band * freq^3
    showf = showf.* x.^2;    % spectrum measured in 1/cm to nm
    showf = showf / max( showf) ;
    plot( 1e7./x(q), showa(q)+offset, 'k-' )
    hold on
    plot( 1e7./x(q), showf(q)+offset, 'r-' )
    plot( [1e7/E(ii, 2), 1e7/E(ii, 2)], [0,1]+offset, 'g:' )
    
    text( min( 1e7./x(q) ), offset+0.15, [ ' ', num2str( ii ),',',...
                                                char( splabels{ii} ),...
                                                ',maxAbs=',num2str( 1e7/v(ii))],...
                                                'Interpreter', 'none' )
    grid on
    set(gca,'xtick',[ 550 : 20 : 750 ])    
    
    ff(:,ii)=showf;  
    aa(:,ii)=showa;    
end

% compute overlaps
% method 1, assume you want to calculate FRET using dipoles
% overlap is computed from spectra normalized to area = 1
% 1. normalize
% 2. overlap, out into a matrix
% 

% figure(2)
for ii = 1:size( a, 2 )
    anorm(:,ii) = a(:,ii) / trapz( x, a(:,ii) );
    fnorm(:,ii) = f(:,ii) / trapz( x, f(:,ii) );     
end
for ii = 1:size( a, 2 )
    for jj = 1:size( a , 2 )
          DF = fnorm(:,ii); % donor fluorescence
          AA = anorm(:,jj); % acceptor absorbance
          J(ii,jj) = trapz( x, AA.*DF );  % overlap
%           subplot( size( a, 2 ),size( a, 2 ), (ii-1)*size( a, 2 ) + jj )
%           plot( x, DF, 'r-' )
%           hold on
%           plot( x, AA, 'b-' )
%           plot( x, anorm(:,ii), ':', 'Color', 0.6*[1,1,1] )
%           text( 20000, 2e-3, [ num2str(ii), '->' ,num2str(jj), '; ', num2str( J(ii,jj) ) ] )
%           pause(eps)
    end
end

% close all
% modify overlaps by Boltzmann weights
J0 = J; % store the original J matrix

T = 293; % [K]
kB = 0.695034800 ; % [1/(cm.K)]

W = ones( 9,9 );
for ii = 1:9
    
    for jj = 1:9
        
        ED = E( jj, 2 ); % energy of donor
        EA = E( ii, 2 ); % energy of acceptor
        
        if ED < EA
            W(jj,ii) = exp( -(EA - ED)/(kB*T) );
            J(jj,ii) = J(ii,jj) * W(jj,ii);
        end
        
    end
    
end

disp( E )
disp( J0 * 1e4)
disp( J * 1e4)
%disp( W )
%disp( (J - J0) ./ J0*100 )
% disp( J0 .* W )
% disp( J(:,9)*1e4 )
% return
%J(1:end-1, end) = 2*J(1:end-1, end);
disp( 'overlaps x 1e4' )
disp( J*1e4 )
fp0 = 'C:\DataVB\Progs\PBS_sims_2022\';
[fn, fp] = uiputfile( [fp0,'new_overlaps_exp_alt_13800_boltz.txt'] )
if fn ~= 0
    csvwrite( [fp,fn], J )
else
    disp( 'overlaps not saved')
end

[fn, fp] = uiputfile( [fp0,'new_overlaps_exp_alt_13800_ctrl.txt'] )
if fn ~= 0
    csvwrite( [fp,fn], J0 )
else
    disp( 'overlaps not saved')
end

% return
% F = ff(:,1);
% A = aa(:,1);
% xnm = 1e7./x(q);
% sortd = sortrows( [xnm, F, A], 1  );
% xnm = sortd(:,1);
% A = sortd(:,3);
% F = sortd(:,2);
% 
% xcm = xnm*1e-7;
% I = 4.5e-18 
% trapz( xcm, A.*F.*xcm.^4 )/trapz( xcm, F )
% Y = 0.23
% t = 1.5 * 1e-9 %[s]
% e = 1.15e-8 %[cm2/mol]
% S = (Y/t)*e
% R = 26.3 % [A]
% n = 1.567
% kappa = 1.13
% 
% C = log(10)*9/(128*pi^5*6.022e23*n^4)
% G = kappa^2 / (R*1e-10)^6;
% kFRET = I*S*C*G
% tfret = 1/kFRET * 1e8

newx = ( 555:1:850)';
for ii = 1:size(ff,2)
    
    
    ff2(:,ii) = interp1( 1e7./x, ff(:,ii), newx );
    aa2(:,ii) = interp1( 1e7./x, aa(:,ii), newx );    
    
end
button = questdlg('Save spectra?',...
'Continue Operation','Yes','No','Help','No');
if strcmp(button,'Yes')
   disp('saving spectra library files...')
   % output the data into files
    absname = ['abs_spectra_',comm,'.csv']
    emname = ['em_spectra_',comm,'.csv']
    pathname = 'C:\DataVB\Progs\PBS_sims_2022\' %fp0
    filep = [ pathname,'\',absname ];
    data2csv( [ newx ,aa2], filep, splabels, numlabels, 'abs' )

    filep = [ pathname,'\',emname ];
    data2csv( [newx,ff2], filep, splabels, numlabels, 'flu' )
 
   disp( '...done')
elseif strcmp(button,'No')
   disp('Canceled file operation')
elseif strcmp(button,'Help')
   disp('Sorry, no help available')
end

disp('end OK')
%===============================================
function f = gpeak( x0, fwhm, a, x )
%
%
%

f = a * exp( -0.5 * ( ( x - x0 ) ./ (fwhm/2.355) ).^2 );

%-----------------------------------------------
function f = data2csv( data, file, splabels, numlabels, header )
%
%
%
fid = fopen( file,'w');
% write header
fprintf( fid, '%s\n', header);
% write 1st line
for ii = 1:max(size( splabels ))+1
    if ii == 1
        textp = 'x,';
    else
        if ii < max(size( splabels ))+1
            textp = [ char( splabels{ii-1} ), ','];
        else
            textp = char( splabels{ii-1} );
        end
    end
    fprintf( fid, '%s', textp);
end
fprintf( fid, '\n');
% write 2nd line
for ii = 1:max(size( splabels ))+1
    if ii == 1
        textp = '[1/cm],'; % units of x-axis
    else
        if ii < max(size( splabels ))+1
            textp = [ num2str( numlabels(ii-1) ), ','];
        else
            textp = num2str( numlabels(ii-1) );
        end
    end
    fprintf( fid, '%s', textp);
end
fprintf( fid, '\n');
for ii = 1:size( data,1)
    
    for jj = 1:size( data,2 )
        
        if jj < max(size( data,2 ))
           textp = [ num2str( data(ii,jj) ), ','];
        else
           textp = num2str( data(ii,jj) ) ;
        end
        fprintf( fid, '%s', textp);
        
    end
    fprintf( fid, '\n');
end
fclose( fid );

disp('...file saved')

%------------------------------------------
function y = carotenoidS1( x00,x )
%
%
%
% ver 2:
a = [0.117, 0.170, 0.763, 0.269];
b = [488.016, 352.427, 1136.202, 1561.849];
x0 = [13917.581, 15025.624, 16439.921, 18126.125] - 14000;
% x0	=[-77.13995474	1244.665093	2465.793407	3591.637644	4835.537891];
% b=[	1349.169935	1349.169935	1349.169935	1349.169935	1349.169935];
% a=[	0.153225397	0.546559673	0.800140134	0.590559293	0.283881307];
%y0=[	0.05631697				];
x0 = x0 + x00;
y = x*0;
for ii = 1:length(x0)
%    y = y + a(ii)*exp( -0.5*( (x-x0(ii)) / (b(ii)/2.355) ).^2);
     y = y + a(ii)*exp( -0.5*( (x-x0(ii)) / b(ii) ).^2); % ver 2
end
y = y ./ x;
y = y / max( y );

%------------------------------------------
function y = carotenoidS1S0f( x00,x )
%
%
%
a = [0.117, 0.170, 0.763, 0.269];
b = [488.016, 352.427, 1136.202, 1561.849];
x0 = [13917.581, 15025.624, 16439.921, 18126.125] - 14000;

% x0	=[-77.13995474	1244.665093	2465.793407	3591.637644	4835.537891];
% b=[	1349.169935	1349.169935	1349.169935	1349.169935	1349.169935];
% a=[	0.153225397	0.546559673	0.800140134	0.590559293	0.283881307];
%y0=[	0.05631697				];
x0 = -x0 + x00;
y = x*0;
for ii = 1:length(x0)
%    y = y + a(ii)*exp( -0.5*( (x-x0(ii)) / (b(ii)/2.355) ).^2);
     y = y + a(ii)*exp( -0.5*( (x-x0(ii)) / b(ii) ).^2); % ver 2    
end
y = y ./ x.^3;
y = y / max( y );

%------------------------------------------------------------
function f = exp_spectrum( vcm, species, abs_em, spctype )
% species = numerical code:
% 1 = CPC a84
% 2 = CPC b84
% 3 = CPC b155
% 4 = APC 660
% 5 = APC 680
% abs_em switches output to absorption/emission spectrum 
% spctype = 0 to output measured/experimental lineshape, or,
%         = 1 the theoretical lineshape  
% spectra are output in wavenumber axis
% ver_ 1.1

global shifttem

if species == 1 % CPC a84
% CPC a84
  if abs_em == 'abs'
     x0 = [572.0904141,	595.7562104,	623.2649391];
     a	=  [34540.55723,	40683.80338,	80209.92383];
     fwhm = [73.4933399,	44.30628051,	34.26792443];
  end

  if abs_em == 'flu'
     x0 = [637.1245623,	665.2971738, 699.911179];
     a =  [0.038663691,	0.012619673,	0.005375097];
     fwhm	= [31.34840537,	41.84466772,	37.50629362];
  end
  
end

if species == 2
% CPC b84	
  if abs_em == 'abs'
     x0 =[	589.2586372,	610.8121206, 629.1168625];
     a	= [27216.01614,	24550.93914,	40802.83227];
     fwhm	= [84.68944833,	33.66579635,	29.2117424];
  end
  if abs_em == 'flu'
     x0	= [641.9609849,	662.3917979,	700.4352219];
     a	= [0.035556482,	0.010583312,	0.006428899];
     fwhm =  [36.51622561,	41.535846,	37.2646874];
  end
end

if species == 3
% CPC b155	
    if abs_em == 'abs'
        x0= [546.5077417,	572.9343465,	600.473332];
        a=	[19565.24511,	42613.22847,	86777.44369];
        fwhm= [60.04713846,	67.50621771,	50.25540478];
    end
    if abs_em == 'flu'
        x0	= [602	624	690]+5;
        a = [	0,	0.028,	0.007];
        fwhm = [15,	40,	100];
    end
end

if species == 4
% APC 660	
    if abs_em == 'abs'
%         x0 = [	595,	630,	654]-1;
%         a = [0.34,	0.46,	0.79];
% %         a(1:2) = a(1:2) * 0.8;
%         fwhm = [40,	40,	25];
%        fwhm = [40,	40,	35];
% alternative in vivo apc
%        x0 =[	602,	630.0020823,	661]-10 + 3 ;
%        a	= [0.6,	0.81,	2.9];
%        fwhm= [	40,	40,	32];       
        x0 =[	602,	630.0020823,	661]-10 + 3 ;
        a	= [0.6,	0.81,	2.9];
        fwhm= [	35,	35,	28];       

    end
    if abs_em == 'flu'
%         x0 =[660, 680,722]+2;
%         a = [0.49, 0.125, 0.09];    
%         fwhm = [25,	30,	50];
        
%        x0 = [	663,	702, 1];
%        a =	[0.85,	0.16, 0];
%        fwhm = [	31,	89.45, 10];

        x0 = [	663,	702, 1];
        a =	[0.85,	0.16, 0];
        fwhm = [	28,	85, 10];
        
    end
end

if species == 5
% APC 680	
    if abs_em == 'abs'
%         x0 = [	595,	630,	654]+3;
%         a = [0.34,	0.46,	0.79];
% %          a(1:2) = a(1:2) * 0.8;
%        fwhm = [40,	40,	25];
% 
       x0 =[	602,	630.0020823,	661]+1  + shifttem + 5;
        a	= [0.6,	0.8,	2.9];
        fwhm= [	35,	35,	28];             
    end
    if abs_em == 'flu'
%         x0 =[678,700,744]-5 ;
%         a = [0.65,0.185,0.1];    
%         fwhm = [30,30,50];
       x0 = [	673.3,	719.5921924, 400]-1 + shifttem-1 + 5;
       a	= [0.70,	0.13, 0];
       fwhm =	[27,	80.2, 11];

    end
end

wnm = 1e7./vcm;
f = wnm * 0;
for ii = 1:3
    
    f = f + a(ii) * exp( -0.5 * ( (wnm - x0(ii) )/(fwhm(ii)/2.355)).^2 );
     
end 
% the parameters are for spectra in nm axis, hence fluorescence must be recomputed to 1/cm:
if abs_em == 'flu' 
    f = f .* wnm.^2;
end
% now recalculate from spectra to lineshapes
if spctype == 1
    
   if abs_em == 'flu' 
      f = f ./ vcm.^3;
   end
   if abs_em == 'abs' 
      f = f ./ vcm;
   end
end

f = f / trapz( vcm, f );