function spec_calc_basic_v1
% given a file of dipole vectors in the format (_ = space):
% indexA_indexB_label(ascii-numerical)_mju_x_y_z[nm]_ex_ey_ez
% where: ex,y,z are components of unit vector in the direction of the dipole
% mju = dipole size in Debye
% x,y,z, are dipole position
% computes dipole couplings and energy transfer rates
% further input: a) file of indexes that are used to address values in the 
% b) file of spectral loverlaps
% ver. 2020-05-14 (c)DB; 
% this procedure has been tested on generated data, using the script 'evaluate_computations.m'  

clear all; 
close all;

global Ktest

% file to analyze:
fid = fopen( 'C:\DataVB\Progs\PBS_sims_2022\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid);    
% fid = fopen(  'C:\DataVB\Progs\PBS_2021\current_path.txt', 'w' );
%           
% for ii = 1:size( fnames, 1 )              
%              
%     dfile = [ fpaths{ii},fnames{ii}] ;
%     disp( [ 'adding file: ', dfile ] )
%     addedfile{ii,1} = dfile;
%     A = textread( dfile ,'%s','delimiter', '\n');   
%     jj = 1;
%     while 1
%         if jj > size(A,1)
%            disp( ['...EOF ', dfile ])
%            break
%         end
%         fprintf( fid, '%s', A{jj} );
%         fprintf( fid, '\n');
%         jj = jj + 1;      
%     end     
%     linecount(ii,1) = jj-1;    
%     
% %     C = read_PDB_coord_02( [IFPath2,'\',dfile] , 1, reslabel );
% end
% fclose(fid);

disp( 'running PBS_w_overlaps.m' )

clear mRe K Vij

% mainpath = 'C:\DataVB\Progs\PBS_2021\tests';
%mainpath = 'C:\Data\2021_04\PBS_2104\data2analyze';
mainpath = fp0;

% load dipole data
[ fname, fpath ] = uigetfile( [mainpath,'\*.*'], 'dipole data' )
if isempty( fname) 
   disp(' no dipoles loaded, exit') 
   return
else
    rawmRe = load( [fpath, fname] );
    index  = rawmRe(:,1);
    chainlabelnum = rawmRe(:,3); % chain label letter, as ASCII number
    mRe = rawmRe(:,4:end); % just the dipole vector data
%     mRe(:,1) = 15;
%     mRe(end-3:end,1) = 2.7; % carotenoid
end
% load index of pigments to assign overlaps
[ fname, fpath ] = uigetfile( [mainpath,'\*.*'] , 'load index' )
if fname == 0 
   disp(' no index file loaded, using generic values ') 
   index = ones( size(mRe,1), 1 );
   return
else
   index = load( [fpath, fname] );
end

% load overlap file
% *** note ovelap is ocded as donor vertical, acceptor horizontal
[ fname, fpath ] = uigetfile( [mainpath,'\new_overlaps_exp_alt.txt'], 'load overlaps' )
if fname == 0 
   disp(' no ovelap file loaded, using generic values ') 
   J = ones( max(index) ) * 2e-4;
   return
else
   J = csvread( [fpath, fname] );
end

%J = J*0.75;
% J(3,1) = 0.6*J(3,1);
% J(3,2) = 0.6*J(3,2);
% J(4,1) = 0.6*J(4,1);
% J(4,2) = 0.6*J(4,2);
% J(5,1) = 0.6*J(5,1);
% J(5,2) = 0.6*J(5,2);


% check consistency of datasets
a = size( mRe ,1 );
b = size( index , 1 );
c = size( J , 1 );
if a ~= b
    disp('inconsistent dipole and index data')
    return
end
if max( index) > max( size( J ) )
    disp('inconsistent index and overlap ')    
    return
end

% % reduce array files to 2 single PBS
% mRe = mRe(1:792, :);
% index = index(1:792); 

% reduce array files to 2 single PBS
% mRe = [mRe(1:792, :); mRe( 1981:1988, : )];
% index = [index(1:792); index( 1981:1988 ) ];

% mRe = [mRe(1:792, :); mRe( 1985:1988, : )];
% index = [index(1:792); index( 1985:1988 ) ];

% mRe = [mRe(1:396, :)];
% index = [index(1:396) ];

% plot3( mRe(:, 2), mRe(:, 3), mRe(:, 4), '.' )
% hold on
% plot3( mRe(22, 2), mRe(22, 3), mRe(22, 4), 'ro' )
% plot3( mRe(397, 2), mRe(397, 3), mRe(397, 4), 'ro' )
% set(gca,'DataAspectRatio', [1 1 1])

% r = mRe( 22 ,2:4 ) - mRe( 397 ,2:4 )
% sqrt( sum( r.^2 ) )
% 
% mRe = [ mRe( 22 , : ); mRe( 397 , : ) ];
% index = [index( 22 ); index( 397 ) ];
% 
% norm( mRe( 1 , 5:7 ) )
% norm( mRe( 2 , 5:7 ) )
% Kij = Kij_local( mRe , 1, index, J )
% kappa_calc_nofile( mRe )
% for ii = 1:size( mRe, 1)
%     
%     v = [mRe(ii, 2:4); mRe(ii, 2:4) + mRe(ii, 5:7) ];
%     plot3( v(:,1), v(:,2),v(:,3), '-' )
%     hold on
% end
% set(gca,'DataAspectRatio', [1 1 1])
% 
% return

% disp( 'index for all residues' )
% disp( index(6*54 + 1) )
% 

% center = mean( mRe(:, 2:4), 1 )
% mRe(:, 2:4) = centercomplex( mRe(:, 2:4) );

% 340.000	1.000	18.034
% 345.000	1.000	58.855
% 364.000	1.000	16.770
% 370.000	1.000	54.265
% 373.000	1.000	19.764
% 381.000	1.000	57.807
% 385.000	1.000	21.357
% 393.000	1.000	59.648
% 398 - 364 : car - bilin
% 398 - 370 : car - bilin, slower
% 399 - 373 : car - bilin
% 399 - 370 : car - bilin, slower
% 397 - 345 : car - bilin
% 397 - 340 : car - bilin, slower
% 400 - 385 : car - bilin
% 400 - 393 : car - bilin, slower

% pivotpt = [44, 43,  30] - center;
% plot3( mRe(:,2), mRe(:,3), mRe(:,4), 'b.' )
% hold on
% plot3( mRe(397:400,2), mRe(397:400,3), mRe(397:400,4), 'r+' )
% plot3( mRe(364,2), mRe(364,3), mRe(364,4), 'rsq' )
% plot3( mRe(370,2), mRe(370,3), mRe(370,4), 'rsq' )

% % 
% % selvec  = [48, 102, 225, 111, 171, 273];
% % plot3( mRe(selvec,2), mRe(selvec,3), mRe(selvec,4), 'go' )
% for ii = 340:400
%     
%     txtcolor = 'k';
%     switch index( ii )
%         case 3
%             txtcolor = ' b ';
%         case 4
%             txtcolor = ' g ';
%         case 5
%             txtcolor = ' r ';
%         case 1
%             txtcolor = ' m ';
%         case 2
%             txtcolor = ' k ';    
%         case 9
%             txtcolor = ' y ';
%     end             
%     plot3( mRe( ii, 2), mRe( ii, 3), mRe( ii, 4), 'sq', 'MarkerFaceColor', txtcolor )
%     hold on
%     text( mRe( ii, 2)+0.1, mRe( ii, 3)+0.1, mRe( ii, 4)+0.1, num2str( ii ), 'Color', 'k' )
%     
% end 
% grid on
% set(gca,'DataAspectRatio', [1 1 1])
% % 
% % 
% return
% set(gca,'DataAspectRatio', [1 1 1])
% hold on
% plot3( pivotpt(1), pivotpt(2), pivotpt(3), 'r+' )
% return

disp( 'overlaps (x 1e4):')
disp(J*1e4)

% assign dipole moments based on index
for ii = 1:length( index )
    mRe(ii,1) = 12;
    switch index(ii)
        case 9
% %              mRe(ii,1) = 5.1;
             mRe(ii,1) = 2.7;
        case 3        
            mRe(ii,1) = 13;            
        case 4        
            mRe(ii,1) = 13;            
        case 5        
            mRe(ii,1) = 13;     
        case 6        
            mRe(ii,1) = 15;            
        case 7        
            mRe(ii,1) = 15;            
        case 8        
            mRe(ii,1) = 15;              
        case 1        
            mRe(ii,1) = 15;            
        case 2        
            mRe(ii,1) = 15;              
            
    end    
end

% mRe(1985:1988, 1) = 5.1;

% h1 = sum( mRe( index == 1, 1 ) );
% n1 = sum( index == 1 );
% 
% h2 = sum( mRe( index == 2, 1 ) );
% n2 = sum( index == 2 );
% 
% h3 = sum( mRe( index == 3, 1 ) );
% n3 = sum( index == 3 );
% 
% h4 = sum( mRe( index == 4, 1 ) );
% n4 = sum( index == 4 );
% 
% h5 = sum( mRe( index == 5, 1 ) );
% n5 = sum( index == 5 );
% 
% h6 = sum( mRe( index == 6, 1 ) );
% n6 = sum( index == 6 );
% 
% h7 = sum( mRe( index == 7, 1 ) );
% n7 = sum( index == 7 );
% 
% h8 = sum( mRe( index == 8, 1 ) );
% n8 = sum( index == 8 );
% 
% % average dipole of core
% (h1+h2+h6+h7+h8) / (n1+n2+n6+n7+n8)
% 
% % average dipole of rod
% (h3+h4+h5) / (n3+n4+n5)
% 
% return

% mRe(340,1) = 15;%1.000	18.034
% mRe(345,1) = 15;%.000	1.000	58.855
% mRe(364,1) = 15;%1.000	18.034
% mRe(370,1) = 15;%.000	1.000	58.855
% mRe(373,1) = 15;%1.000	18.034
% mRe(381,1) = 15;%.000	1.000	58.855
% mRe(385,1) = 15;%1.000	18.034
% mRe(393,1) = 15;%.000	1.000	58.855
% 
% mRe( 340, 1) = 15;
% mRe( 364, 1) = 15;
% mRe( 373, 1) = 15;
% mRe( 385, 1) = 15;
% mRe(397,1) = 0;
% mRe(398,1) = 0;
% mRe(399,1) = 0;
% mRe(400,1) = 0;
% mRe = mRe(1:54, :);
% mRe(:,1) = mRe(:,1) / (1.4^0.5); 
% mRe(398,1) = 2.7;
% mRe(364,1) = 0;
% mRe(370,1) = 0;

% mRe(364,5:7) = rot_vector_local(mRe(364,5:7),[ [1,0,0]; [0,1,0]; [0,0,1] ], [0, 0, -20]) ;

plot3( mRe(:,2), mRe(:,3), mRe(:,4), 'bo' )
hold on
grid on
set(gca,'DataAspectRatio', [1 1 1])
%plot3( mRe(54:end,2), mRe(54:end,3), mRe(54:end,4), 'ro' )
plot3( mRe(50:60,2), mRe(50:60,3), mRe(50:60,4), 'ro' )
return

% mRe2 = mRe;
% mRe2(:,2) = mRe2(:, 2) + 15;
% mRe2(:,3) = mRe2(:, 3) + -10;
% 
% plot3( mRe2(:,2), mRe2(:,3), mRe2(:,4), 'ro' )
% 
% view(2)
% print( [fpath,'test_300dpi.png'], '-dpng', '-r300');
% 
% mRe = [mRe; mRe2];
% index = [index; index];

% plot3( mRe(397,2), mRe(397,3), mRe(397,4), 'r+' )
% plot3( mRe(1:54,2), mRe(1:54,3), mRe(1:54,4), 'g.' )

% startpts = [ (0:5)*54+1 ];
% plot3( mRe(startpts,2), mRe(startpts,3), mRe(startpts,4), 'r+' )

% % rod 1
% W = mRe( 1:54, 2:4);
% [a,m] = axis_by_svd( W );
% 
% W = W - repmat( m , size( W,1),1 );
% 
% for ii = 1:size( W ,1)
%     
%     W( ii,: ) = rot_vector_local( W(ii,:) , [ a';[0,0,1];[0,0,1]] , [45,0,0]) ;
%     
% end
% 
% mRe( 1:54, 2:4) = W + repmat( m , size( W,1),1 );

% plot3( mRe(:,2), mRe(:,3), mRe(:,4), 'r.' )

% % rod 1
% W = mRe( 1:54, 2:4);
% [a,m] = axis_by_svd( W );
% plot3( [m(1),m(1) - 10*a(1)], [ m(2), m(2) - 10*a(2) ] , [ m(3), m(3) - 10*a(3) ],'r-' )

% translate rod 1
% mRe(1:54, 2:4) = blockshift( mRe(1:54, 2:4), 1);
% rotate rod 1
% W = mRe( 1:54, 2:4);
% W = W - repmat( pivotpt , size( W,1),1 );
% for ii = 1:size( W ,1)
%     
%     W( ii,: ) = rot_vector_local( W(ii,:) , [[0,0,1];[0,0,1];[0,0,1]] , [30,0,0]) 
%     
% end
% W = W + repmat( pivotpt , size( W,1),1 );
% 
% mRe( 1:54, 2:4) = W;

% plot3( mRe(:,2), mRe(:,3), mRe(:,4), 'r.' )
% 
% mRe(397,2) = mRe(397,2)+0.05;
% mRe(397,3) = mRe(397,3)+0.05;

% plot3( mRe(397,2), mRe(397,3), mRe(397,4), 'g+' )

% return

% tempout = [index, mRe(:,1)];
% save( 'C:\DataVB\data2analyze\PBS_OCP_complete\temp_out_dips.txt', 'tempout', '-ascii' )
% return
% 
% if kk == 1
% for ii = 1:size( mRe,1)
%     
%     %eay =  rot_vector( ea0y, [1,0,0; 0,1,0;0,0,1] ,ang_a)
%     mRe( ii, 2:4 ) =  rot_vector( mRe( ii, 2:4 ) , [1,0,0; 0,1,0;0,0,1] , [90, 0,0] )    ;
%     mRe( ii, 5:end ) =  rot_vector( mRe( ii, 5:end ) , [1,0,0; 0,1,0;0,0,1] , [90, 0,0] )   ;     
% end
% end

% mRe = mRe( 325:end, :);
% index = index( 325:end);

% disp('kappa')
% kappa_calc_nofile( mRe );
% disp('computing coupling')
% V = Vij_calc_nofile_local( mRe, 1.0 ); %excitonic couplings = off-diagonal part of Hamiltonian


% U = reshape( V, [], 1);
% U = nonzeros( U );
% U = U( U > 50 );
% hist(U)
% 
% [a1,a2,a3] = find( V > 100 )
% % U = sort(U);
% % disp( U( 1:30));
% return
disp('computing FRET rate constants')
Kij = Kij_local( mRe, 1.0, index, J );

% pick the carotenoid-acceptor columns
% CC = Kij( :, end-3:end  );
% TC = 1./CC;
% % disp transer-to-car
% [u,v,w] = find( TC < 50 );
% 
% for ii = 1:length( u )
%     
%     disp( [ '[',num2str( u(ii) ),',',num2str( v(ii) ),']: ',num2str( TC(u(ii),v(ii) ) ) ] )
%     
% end

% return
% disp(Kij(364, 398) )
% return
% 398 - 364 : car - bilin
% 398 - 370 : car - bilin, slower
% 399 - 373 : car - bilin
% 399 - 370 : car - bilin, slower
% 397 - 345 : car - bilin
% 397 - 340 : car - bilin, slower
% 400 - 385 : car - bilin
% 400 - 393 : car - bilin, slower

% kbc = 1/100; % [1/ps]
% Kij(364, 398) = kbc;
% % Kij(370, 398) = kbc/3;
% Kij(373, 399) = kbc;
% % Kij(370, 399) = kbc/3;
% Kij(345, 397) = kbc;
% % Kij(340, 397) = kbc/3;
% Kij(385, 400) = kbc;
% Kij(393, 400) = kbc/3;

% computation test section
% donor = 214
% acceptor = 217
% 
% disp( Kij( donor, acceptor ) )
% disp( index( donor ) )
% disp( index( acceptor ) )
% J( index( donor ), index( acceptor )  )
% vecdonor = mRe( donor, :)
% vecacceptor = mRe( acceptor, :) 

% return

K = Kij;
% rij = K * 0;
% [ fname, fpath ] = uiputfile( [mainpath,'\*.*'] , 'save the matrix of time constant' )
% if fname ~= 0
%     
%     csvwrite( [fpath, fname], 1./K );
%     
%     for ii = 1:size( mRe , 1 )
%         for jj = 1:size( mRe , 1 )
%             rij(jj,ii) = norm( mRe(jj,2:4) - mRe(ii,2:4) );
%         end
%     end
%     csvwrite( [fpath, [fname(1:end-4),'_R.txt'] ], rij );
%     
% disp('...rate constants saved')
% end
% return
% for ii = 1:size( Kij, 1)
%     
%     for jj = 1:size(Kij,2)
%         
%         if jj ~= ii
%            disp( ['k(',num2str(ii),'->',num2str(jj),') = ', num2str(1/Kij(ii,jj)),' ps'] )
%         end
%     end
% end


% return
% 

% % * modify rate matrix
% add acceptor
% index = [index; 10];
% K = [ [ K; zeros( 1, length( index)-1)], zeros( length(index), 1)];
% q = ( index > 5 & index < 9 );
% K(q, end ) = 1/15; %[ps];
% disp( [index, q ] )

% add acceptor to an array
% index = [index; 10; 10; 10; 10; 10];
% K = [ [ K; zeros( 5, length( index)-5)], zeros( length(index), 5)];
% q = ( index > 5 & index < 9 );
% ipos = (1:length( index) )';
% np = 396;
% % class(q)
% for ii = 1:5
%     q2 = q;
%     rcrange1 = (ii-1)*np + 1;
%     rcrange2 = (ii-1)*np + np;
%     q2( ipos < rcrange1 | ipos > rcrange2) = false;
%     K( q2, end-5+ii ) = 1/15; %[ps];
% end

% now turn off all but the last RC
% K( :, end-4:end-1) = 0;

% disp( [ index, q ] )
% % clear diagonal
% for ii = 1:size(K,1)
%     
%     K(ii,ii) = 0;
%     
% end
% % make diagonal meaningful
% edit diagonal
for ii = 1:size(K,1)
    % clear present value
    K(ii,ii) = 0;
    % prepare the S1 decay rate
    k00 = 1/1500;
    % treat special cases
    switch index(ii)
        case 4
            k00 = 1/900;
        case 9
            k00 = 1/3.3;
%               k00 = 0;    
        case 10
            k00 = 0; % trap does not decay
    end
    % make diagonal meaningful
    K(ii,ii) = - (sum( K(ii,:) )+k00) ;
end

% csvwrite( 'C:\Data\2021_08\testK.txt', K(:, end-10:end) );
% return


% K(401,401) = 0;
% K(364, 398) = 0;
% K(370, 398) = 0;

% for ii = size(K,1)-3:size(K,1)
%     
%     K(ii,ii) = - (sum( K(ii,:) )+1/3.3) ;
%     
% end



% save( 'C:\Data\2021_06\data2analyze\PBS_OCP_complete\tempK.txt', 'K', '-ascii' )
% return

% transpose
Ktest = K';
% Kij is ready to use now...
size(  Ktest )

% startpts = [ (0:5)*54+1 ];
% selvec  = [48, 102, 225, 111, 171, 273];

% startpts = 54;
% startpts = [46,47,48];
% startpts = [48,225,114];
% startpts = [60,114,168];  % down down
% startpts = 54+54+54+(1:9);% down down
%startpts = 54+54+72+54+(1:9); % up-up
% startpts = 54+54+72+54+1; % up-up
% startpts = [132,186,240];  % up up
% startpts = ceil( rand(10,1)*396 ) ;

% arrays
startpts = 74;

classes = unique( index);

for ii = 1:length( startpts )
    
    disp( [ 'bilin # ', num2str( startpts(ii) ),' is class ', num2str( index( startpts(ii) ) ) ] )
    
end

% return
% disp( index( [46;47;48] , 1)  )
% disp( ['bilin types to excite: ', num2str( index( startpts ) ) ] )
% return
% for ii = 1:length( t )
%     
%     bar( b, Y(ii,:) )
%     pause(0.001)
%     
% end

clck = clock;
disp( ['ode started ', num2str( clck(1) ),'-',num2str( clck(2) ),'-',num2str( clck(3) ),', ',...
        num2str( clck(4) ),':', num2str( clck(5) ) ] );

for uu = 1:length( startpts)
    
    v0 = zeros( 1, size( Ktest, 1) );
    disp( ['excited pigment no. ', num2str( startpts(uu) ) ])
    v0( startpts(uu) )  = 1;
%     v0( 1:24 )  = 1/24;
    % modeling is done using the ODE solver
    tic
    disp('ode45:')
    init_cond = v0;%[ 1, 0, 0, 0];
    [t,y] = ode45( @model , [0, 300], init_cond);
    disp('- ode finished')
    toc
    %.........
    disp( '- 2nd ode' )
    tic
    init_cond = y(end, :);
    [t2,y2] = ode45( @model , [0, 16.6, 30, 54, 100, 166, 300, 540, 1000, 1660, 3000, 5400, 10000], init_cond);
    toc
    
    n1 = length(t);
    t = [t;t2+300];
    y = [y;y2];
    
    %........
%     figure(1)
%     plot( t, y(:, (1:2:end) ), '-' )
%     hold on
%     plot( t, sum(y,2), 'k-' )
%     hold on
    
    size(y)
    size( index')

%     N = length( t );
    floor( log10( n1 ) );
    q = logspace( 2, log10( n1 ) );
    
    timepick = [(1:1:10), (15:85), round(q), (length(t)-11:1:length(t))]; %[ (1:2:100), (150:50:500), (600:200:1000) ];
    
    % return
    % for ii = 1:length( classes )
    %     
    %     c = index';
    %     Y(: ,ii) = sum( y( timepick, c == classes( ii ) ), 2 );
    %     selA(:, ii) = a(:, classes(ii) == sindex ).* x;;
    %     selF(:, ii) = f(:, classes(ii) == sindex ).* x.^3;
    % end
    % mm = max( max(  selF ) );
    % selF = selF / mm;
    % %selF = selF .* x.^2;
    % plot( t( timepick) , Y, ':')
    % 
    % F = selF * Y';
    % figure(3)
    % plot( 1e7./x, F', '-')
    % 
    % size( t(timepick) )
    % size( F' )
    % size( 1e7./x )
    % 
    % 
    % out = [ [0;t( timepick)], [1e7./x'; F' ] ];
    % 
    % save( 'E:\Work\2021_05\sim.txt', 'out', '-ascii' )
    
    % return
    pt = t( timepick);
    if uu == 1
        py = y( timepick, :);
    else
        py = py + y( timepick, :);
    end
    
    kin_name = 'ctrl_cpc13_apc15_tem15_test_rev_0';
    
    disp( 'saving current data')
    out = [ [999, index', -999 ]; [ t( timepick, :) , y( timepick, :) ], sum( y( timepick, :) , 2) ]; 
    csvwrite( [fpath,kin_name,'_',num2str(startpts(uu)),'.txt'], out  )
    
    tempy =  y( timepick, :);
    
    for ii = 1:length( classes)
    
        C(:,ii) = sum( tempy(:, index == classes(ii)), 2 );
        outC(:,ii) = [ classes(ii); C(:,ii)];
        
    end
    outC = [ [999;t( timepick, :)], outC,  sum( outC, 2 ) ];
    save( [fpath,kin_name,'_by_class_',num2str(startpts(uu)),'.txt'], 'outC','-ascii'  )
    outC = [];
    out = []; 
end % of 'for' over 'uu'

if length( startpts  ) > 1
    disp( ' saving the sum of all excitations: ')
    out = [ [999, index', -999 ]; [pt, py], sum( py, 2) ]; 
    csvwrite( [fpath,kin_name,'_sum.txt'], out  )
end
    
V = Kij;
[ fname, fpath ] = uiputfile( [mainpath,'\*.*'] , 'save the matrix of rate constant' )
if fname ~= 0
    
    csvwrite( [fpath, fname], K );
disp('...rate constants saved')
end

% 
% return
% Vup = max( max( abs( V )) );
% kappa = kappa_calc_nofile( mRe )
% 
% subplot( 1,2,kk)
% figure(1)
% vecplotterA(mRe(:,:), 'b', 0)
% % subplot( 1,2,1 )
% % vecplotterA(carmRe, 'r', 0)
% set(gca,'DataAspectRatio', [1 1 1])
% % axis( [-1.5 1.5 -1.5 1.5 -1.5 1.5 ]  )
% xlabel( 'X' )
% ylabel( 'Y' )
% zlabel( 'Z' )

% V = Vij_calc_nofile_local( mRe, 1.0 ); %excitonic couplings = off-diagonal part of Hamiltonian
% Vcar = V( :, end-3:end );
% disp(' couplings' )
% for ii = 1:length( u )
%     
%     disp( [ '[',num2str( u(ii) ),',',num2str( v(ii) ),']: ',num2str( Vcar( u(ii),v(ii) ) ) ] )
%     
% end

clck = clock;
disp( ['computations completed: ', num2str( clck(1) ),'-',num2str( clck(2) ),'-',num2str( clck(3) ),', ',...
        num2str( clck(4) ),':', num2str( clck(5) ) ] );

disp( 'end OK ')
return


etrace = zeros(50, 5);

for jj = 1:1%size( mRe, 1 )
curpos = jj;
prev = 0;
ii = 1;

etrace(1,jj) = curpos;

while 1    

    % pick a dipole & find the most probable path from it
    u = Kij( curpos,: ) ;% row = donor, column = acceptor
    u( u == 0 ) = -0.001;
    z = sort( u );           
    while 1
 
        [mval, mix] = max( u );    
        if mix == prev
           u( u == max(u) ) = -0.001; 
        else
           break
        end

    end
    
%     disp( mix )
%     disp( u(mix) )
    
    start = mRe( curpos, 2:4 );
    next =  mRe( mix, 2:4 );
    vec = [start; next];
    
    
%     plot3( vec(:,1), vec(:,2), vec(:,3),'k-')
%     text( mean( vec(:,1) ),mean( vec(:,2) ),mean( vec(:,3) ), num2str( 1/u(mix) )  )
%     
%     disp( u )
%     disp( 1./u )    
    plot_direction( start, next, 20 , classes( 1/mval ) );
%     disp( [ num2str(curpos),'->', num2str( mix ),' = ', num2str(1/mval), 'ps'  ] )

    etrace(ii+1,jj) = mix;
    ii = ii + 1;
    
    prev = curpos;
    curpos = mix;

    if ii == 50
       break
    end
       
end

end
print( [fpath,'test_300dpi.png'], '-dpng', '-r300');

etrace = [(1:50)', etrace];
disp(etrace)
return
% 
% % r = 1;
% % for ii = 1:6
% %     
% %     for jj = 1:size( Kij, 2 )
% %         
% %         v = [ mRe( ii, 2:4); mRe( jj,2:4 ) ];
% %         plot3( v(:,1), v(:,2), v(:,3), '-'  )
% %         text( mean(v(:,1) ) , mean(v(:,2) ), mean(v(:,3)), num2str( 1/Kij( ii, jj ) )  )
% %         
% %     end
% % end
% % 
% % return
% % return
% % plot labels
% % for ii = 1:size( mRe, 1)
% %     
% %     text( mRe(:, 2)+0.3, mRe(:, 3)+0.3, mRe(:, 4)+0.3, char( chainlabelnum ) )
% %     
% % end
% 
% veccount = size( mRe, 1); % number of dipole vectors
% % design plot format for couplings or FRET rates etc.
% % try finding & highlighting strongest couplings
% % cutoff(1) = 20;
% % cutoff(2) = 45;
% % cutoff(3) = 1e3;
% 
% cutoff(1) = 0.02;
% cutoff(2) = 5000;
% 
% h = (1:veccount)';
% %wdth = [1,3,5];
% wdth = [3,3,3,3];
% lcol = { [1,0.8,0.8], [1,0.4,0.4], [1,0,0]  };
% 
% disp( '---------------------------------------------------------------')
% disp( 'significant dipole-dipole couplings [1/cm] ')
% disp( '---------------------------------------------------------------')
% for ii = 1:veccount-1
%     
%     u = V( :, ii );
%     
%     for jj = 1:length(cutoff)-1
%         
%         s = abs(u) >= cutoff(jj) & abs(u) < cutoff(jj+1);
%         t = h( s );
%                
%         if ~isempty(t) 
%             
%             for kk = 1:length(t)
%                 if t(kk) > ii
%                    p = [ mRe(ii, 2:4 ); mRe( t(kk), 2:4 ) ];
%                    disp( [ num2str(ii),'-',num2str( t(kk) ),': ', num2str( V( ii, t(kk) ) ) ] )
%                    plot3( p(:,1), p(:,2), p(:,3), 'g-', 'LineWidth', wdth(jj), 'Color', rainbow( abs(V( ii, t(kk) )), cutoff(1), Vup )  )
%                    hold on
%                 end
%             end
%             
%         end
%     end
%     
% end
% 
% % hb = (10:5:180)';
% % W = abs(V( abs(V) >= 10 ));
% %figure(3)
% W = V;
% % [a,b] = hist( reshape( W, [], 1 ), hb );
% [a,b] = hist( reshape( W, [], 1 ));
% disp( fname )
% disp([b(:),a(:)])
% % end
% % subplot( 1,2,2)
% % % figure(2)
% % for ii = cutoff(1):Vup
% %     
% %     plot( [0,1], [ ii,ii], 'Color', rainbow( ii, cutoff(1), Vup ), 'LineWidth', 4 )
% %     hold on

return

% merged
%********************************
DC = 1.1;
X = (400:1:700)';
M = 100;
% vecplotter(mRe, 'r', 0)
set(gca,'DataAspectRatio', [1 1 1])
% axis( [-1.5 1.5 -1.5 1.5 -1.5 1.5 ]  )
xlabel( 'X' )
ylabel( 'Y' )
zlabel( 'Z' )
% print( 'C:\Data\2021_04\vect2.png', '-dpng', '-r300');

eA = X*0;
eCD = eA;
eX = X;
 return
% 
%%### plot setup
figure(2)
a1 = subplot(4,2,1);
set(a1,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a2 = subplot(4,2,2);
set(a2,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a3 = subplot(4,2,3);
set(a3,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a4 = subplot(4,2,4);
set(a4,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');

a5 = subplot(4,2,5);
set(a5,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a6 = subplot(4,2,6);
set(a6,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'add', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a7 = subplot(4,2,7);
set(a7,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');
% title('absorption spectrum');
a8 = subplot(4,2,8);
set(a8,'ButtonDownFcn', 'set(gcf, ''UserData'', 0)', 'NextPlot', 'replacechildren', 'UserData', 1, 'XGrid', 'on');

title('---');

set(gcf, 'UserData', 1);
% rect = [left, bottom, width, height]
set(gcf, 'Units', 'normalized', 'Position', [0.55 0.2 0.40 0.7]);
%###end of plot setup

t = 1;
A0 = X * 0;
CD0 = A0;

wb = waitbar(0,'Computing inhomogeneous broadening, please wait...');
 
while t <= M && (get(gcf, 'UserData'))
             
%% set inhomogeneous broadening 
    if M == 1
       iHBE = SE;
    else
       iHBE = IHBSE(SE, iHBw, 'FWHM');
    end
            
    S = CalcSpcACD_HB2(mRe, iHBE, DC, HB, X, 0);   
%             SLD = CalcSpcLD(mRe, iHBE, DC, HB, X, 0, [0,0,1]);        
    A0 = A0 + S(:, 2);
    CD0 = CD0 + S(:, 3);
%             LD0 = LD0 + SLD(:,2);      
    A = A0 / t;
    scale = 1 ;%max( eA)/max(A);
    CD = CD0 / t;
    A = A * scale;
    CD = CD * scale * 0.01;
%             LD = LD0 / t;

    figure(2)
    get(gcf, 'UserData');
    d = get(gcf, 'Children');
    plot(X, A, 'b-', eX, eA, 'k-', 'Parent', d(7), 'LineWidth', 2);
    plot(X, CD,'b-', eX, eCD, 'k-', 'Parent', d(5), 'LineWidth', 2);
    t = t + 1;
    waitbar( t/M );
    pause(eps)  
end 
% print( 'C:\Data\2021_04\sim2.png', '-dpng', '-r300');

out = [X, A, CD];
save( 'C:\Data\2021_04\sim_spc.txt', 'out', '-ascii' );
close(wb)
disp('end OK')

%=============================================================
function VP = vecplotterA(mRe, co, numlabl)

R = mRe(:, 2:4);
e = mRe(:, 5:7);
d = mRe(:, 1);
N = size(d,1);
sf = 0.1;
LC = co;
MC = ['sq',co];

minR = min(min(R));
maxR = max(max(R));
% plotbox = [1.1 * (minR+0.1), 1.1 * (maxR+0.1), 1.1 * (minR+0.1), 1.1*(maxR+0.1), 1.1 * (minR+0.1), 1.1 * (maxR+0.1)];

for ii = 1:N
    V  = [R(ii,:) - sf * 0.5 * d(ii) * e(ii,:); R(ii,:) + sf * 0.5 * d(ii) * e(ii,:)];       
    EP = R(ii,:) + sf * 0.5 * d(ii) * e(ii,:);           
    C = R(ii,:);
    plot3(C(1,1),C(1,2),C(1,3),'ok', 'MarkerSize', 2, 'MarkerFaceColor', 'k' ) ;
%     axis( plotbox );
    grid on
    hold on
%     plot3(EP(1,1),EP(1,2),EP(1,3), MC,'MarkerFace', LC) ;
    hold on
    plot3(V(:,1),V(:,2),V(:,3), LC,'LineWidth', 2)  ;   
    hold on
%    plot3([0;R(1,1)], [0;R(1,2)], [0;R(1,3)],'k', 'LineWidth', 1)  ;   
%    text([0;0.5*R(1,1)], [0;0.5*R(1,2)], [0;0.5*R(1,3)], num2str( norm(R(1,:)) ) )  ;       
%    hold on   
end

if numlabl ~= 0
    for ii = 1:N
         text(R(ii,1)+0.1,R(ii,2)+0.1,R(ii,3)+0.1, num2str( numlabl(ii) ) )  ;       
         hold on   
    end
end

%----------------------------------------------------------------------
function v = rainbow(x, xmin, xmax)
%
%
%
figure(1)



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
function ofd = Vij_calc_nofile_local(vecdata, DC)
% ver. 2013-20-08
% computes the off-diagonal elements (V(ij)) of matrix Hamiltonian
% based on the coordinate file given in "FName"
% outputs the Vij as a matrix with zeros on the diagonal
% DC = dielectric constant (usu. epsilon)
R = vecdata(:,2:4);
e = vecdata(:,5:7);
m = vecdata(:,1);
N = size(m,1);
% diferent, but equivalent, way of computing the off-diagonal elements
Vij = [];
for ii = 1:N
    for jj = ii+1:N
        rij = R(jj,:) - R(ii,:);
        if norm(rij) ~= 0 
          %   screen = 2.68 * exp(-0.27*norm(rij)) + 0.54;
          Vij(jj,ii) = 5.04 / DC * ( m(ii) * m(jj) * dot(e(ii,:), e(jj,:)) / norm(rij)^3 - 3*m(ii)*m(jj)*dot(e(ii,:),rij)*dot(e(jj,:),rij) / norm(rij)^5 );    
        else
            Vij(jj,ii) = 0;
        end
    end
end 
ofd = [Vij,zeros(N,1)] + [Vij,zeros(N,1)]' ; % off diagonal elements, the final matrix

%---------------------------------------------------------------------
function Kij = Kij_local(vecdata, DC, index, J)
% ver. 2013-20-08
% computes the off-diagonal elements (V(ij)) of matrix Hamiltonian
% based on the coordinate file given in "FName"
% outputs the Vij as a matrix with zeros on the diagonal
% DC = dielectric constant (usu. epsilon)
R = vecdata(:,2:4);
e = vecdata(:,5:7);
m = vecdata(:,1);
N = size(m,1);
% diferent, but equivalent, way of computing the off-diagonal elements
Vij = [];
for ii = 1:N
    for jj = 1:N
        rij = R(jj,:) - R(ii,:);
        if norm(rij) ~= 0 
          %   screen = 2.68 * exp(-0.27*norm(rij)) + 0.54;
           Vij(jj,ii) = 5.04 / DC * ( m(ii) * m(jj) * dot(e(ii,:), e(jj,:)) / norm(rij)^3 - 3*m(ii)*m(jj)*dot(e(ii,:),rij)*dot(e(jj,:),rij) / norm(rij)^5 );
           Kij(jj,ii) = (Vij(jj,ii))^2 * 1.18 * J( index(jj), index(ii) );
           %disp( J( index(jj), index(ii) ) )
%            disp( Vij(jj,ii) * (DC / 5.04) * 1/(m(ii)*m(jj)) * norm(rij)^3 )
           
        else
           Vij(jj,ii) = 0;
           Kij(jj,ii) = 1/1500;           
        end
    end
end 
%Kij = Kij,zeros(N,1)] ; % off diagonal elements, the final matrix

%-------------------------------------------------------------------
function ofd = kappa_calc_nofile(vecdata)
% ver. 2013-20-08
% computes the off-diagonal elements (V(ij)) of matrix Hamiltonian
% based on the coordinate file given in "FName"
% outputs the Vij as a matrix with zeros on the diagonal
% DC = dielectric constant (usu. epsilon)
R = vecdata(:,2:4);
e = vecdata(:,5:7);
m = vecdata(:,1);
N = size(m,1);
% diferent, but equivalent, way of computing the off-diagonal elements
Vij = [];
for ii = 1:N
    for jj = ii+1:N
        rij = R(jj,:) - R(ii,:);
        rij = rij / norm(rij);
        if norm(rij) ~= 0 
            Kij(jj,ii) =  dot(e(ii,:), e(jj,:)) - 3*dot(e(ii,:),rij)*dot(e(jj,:),rij);
            %Kij(jj,ii) = Kij(jj,ii) * 5.04 / 1.34 * m(jj)*m(ii) * 1/norm( R(jj,:) - R(ii,:) )^3;
         else
            Kij(jj,ii) = 0;
        end
    end
end 
ofd = [Kij,zeros(N,1)] + [Kij,zeros(N,1)]' ; % off diagonal elements, the final matrix

%--------------------------------------------------------------------
function f = plot_direction( start, stop, maxwidth, greenperc )
%
%
%
P0 = start;
r = norm( [start;stop] );
d = ( stop - start ) / r;
N = 300;
greenperc = greenperc/100;
for ii = 1:N
   
   %LC = [ ii/N, greenperc/100, 1-ii/N ] ;
   LC = [ greenperc, 0, 1-greenperc ] ;   
   P = P0 + d * (r / N); 
   v = [P0;P];
 %  plot3( v(:,1), v(:,2), v(:,3), '-', 'Color', LC, 'LineWidth', (N+1-ii)/N * maxwidth)    
   plot3( v(:,1), v(:,2), v(:,3), '-', 'Color', LC, 'LineWidth', (N+1-ii)/N * maxwidth)     
   P0 = P; 
   
end

%---------------------------------------------------------------------
function f = classes( x )
%
%
%

support = [ 1e13, 1500, 150, 15, 1.5, 0 ];
w = 1:length( support );

c = [   1,   10,   20,  30,  60, 100];

y = interp1( support, w, x );
f = c(ceil(y));

%--------------------------------------------------------------------
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
function [a,m] = axis_by_svd( P )
%
%
%
m = mean( P, 1 );
h = P - ones( size( P, 1 ), 1 )*m;

[ U , S , V] = svd( h, 0 );
[s,ii] = max( diag(S));
a = V(:,ii);

%---------------------------------------------------------------------
function W = blockshift( V, dv )
% given a set of coordinates of a point cloud (block)
% finds the long axis of the cloud
% and shifts it by normal(radnom) step characterized by std=sstd
% along the axis

[a,m] = axis_by_svd( V );
A = repmat( a', size(V,1),1  );
W = V + dv*A;

%--------------------------------------------------
function newc = centercomplex( C )
%
%
%
center = mean( C, 1 );
newc = C - repmat( center, size( C,1), 1 );

%------------------------------------------------------------
function NV = rot_vector_local(V, U, angles) 
% performs a sequential rotation of given vector about 3 axes: U = [axis 1; axis 2; axis 3]
% angles = [alpha, beta, gamma], in degrees
% ver 1.0., 2014-V-13

d = length(V);
ax1 = U(1,:);
ax2 = U(2,:);
ax3 = U(3,:);

% rotation about axis 1 by alpha degrees
V1 = V * CCW_rotM_local(angles(1), ax1);
% k = ax1 / norm(ax1);
% V1 = V * cos(angles(1)) + cross(k, V)*sin(angles(1)) + k*dot(k, V)*( 1-cos(angles(1)) );
% rotation about axis 2 by beta degrees
V2 = V1 * CCW_rotM_local(angles(2), ax2);
% k = ax2 / norm(ax2);
% V2 = V1 * cos(angles(2)) + cross(k, V1)*sin(angles(2)) + k*dot(k, V1)*( 1-cos(angles(2)) );

% rotation about axis 3 by gamma degrees
V3 = V2 * CCW_rotM_local(angles(3), ax3);
% k = ax3 / norm(ax3);
% V3 = V2 * cos(angles(3)) + cross(k, V2)*sin(angles(3)) + k*dot(k, V2)*( 1-cos(angles(3)) );

NV = V3;

%-----------------------------------------------------------
function R = CCW_rotM_local(deg,u)
% ver. 2013-20-08
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%Input is in degrees.
% deg2rad:
rads = deg / 180 * pi;

R=eye(3);
u=u(:)/norm(u);
x=rads; %abbreviation

for ii=1:3

    v=R(:,ii);
    R(:,ii)=v*cos(x) + cross(u,v)*sin(x) + (u.'*v)*(1-cos(x))*u;
      %Rodrigues' formula     
end 

%--------------------------------------------------------------------