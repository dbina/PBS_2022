function read_coord_1
%
%
%
clear all; close all; 


% file to analyze:
fid = fopen( 'C:\Data\Progs\PBS_2021\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid);     

disp( 'starting tres_reader, selecting the tres data folder')
% get the directory and read data
directory_name = uigetdir( fp0 )
% make the list of files
d = dir( directory_name );
dsize = max( size( {d.name} ) );
str0 = {d.name};

jj = 1;
for ii = 1:dsize
     
   if ~isempty( findstr( str0{ii}, '.xyz' ) )
       str{jj} = str0{ii};
       jj = jj + 1;
   end

end
[ss,v] = listdlg('PromptString','Select a file:',...
                 'SelectionMode','multiple',...
                 'ListString',str);

IFPath2 = directory_name;
uu = 1;
for kk = 1:length( ss ) 

    dfile = str{ss(kk)};
    c = load( [ IFPath2,'\',dfile] );
    x = c(:,1);
    y = c(:,2);
    z = c(:,3);
    chlab = c(:,4);
%     disp( size( x ) );
    length(x)
    
    % plot all atoms as points 
    pcolor = [ kk/length( ss ), 0.2, 1-kk/length( ss )];
    figure(1)   
%     for uu = 1:length( x )

%         plot3( x(uu), y(uu), z(uu), '.', 'Color', pcolor)
        plot3( x, y, z, '.', 'Color', pcolor)
       
        hold on
%         pause(0.1)
%     end
    nat = 94;
    npig = length(x) / nat;
    apig = nat; %length( x ) / npig
    
%     cm(1,:) = [29.26	-26.34	12.23];
%     cm(2,:) = [-17.56	-10.77	17.32];
%     cm(3,:) = [-15.02	-43.93	8.43];
%     cm(4,:) = [-33.52	-20.71	-12.39];
%     cm(5,:) = [15.36	-13.83	-17.49];
%     cm(6,:) = [7.22	-46.12	-8.34];
%     
%     v(1,:) = [0.81	-0.67	-0.69];
%     v(2,:) = [-0.56	-0.11	0.96];
%     v(3,:) = [0.42	-0.2	1.85];
%     v(4,:) = [-1.84	-0.74	0.94];
%     v(5,:) = [1.02	-0.35	-1.14];
%     v(6,:) = [-0.91 0.43  -2.18	];
%     
%     ring1 = [2, 3, 4, 5, 6 ];
%     ring2 = [14, 15, 16, 17, 18 ];
%     ring3 = [23, 24, 25, 26, 27];
%     ring4 = [33, 34, 35, 36, 37 ];
    
    %aroi = [ ring1, ring2, ring3, ring4, 13, 1, 32 ];
    %aroi = [5, 1:42];
    % select atoms belonging to the aromatic part 
    % for computation of the dipoles
%     aroi = [5,10,11,14,15,16,17,22,23,24,25,28,29,30,31,34,35,36,37,38,39,40,41,42,43];
%      aroi = (41:nat-20);
      aroi = (21:40);
    for ii = 1:npig
        
        B{ii} = c( 1 + (ii-1)*apig : ii*apig, : );
        A{ii} = B{ii}( aroi ,:);
        figure(1) 
        plot3( A{ii}(:,1), A{ii}(:,2), A{ii}(:,3), 'b+' )
    end
    
    for ii = 1:npig
        
        figure(1)
        P = A{ii}(:,1:3) ;%B{1}
        
        m = mean( P, 1 );
        h = P - ones( size( P, 1 ), 1 )*m;
        
        [ U , S , V] = svd( h, 0 );
        [s,i] = max( diag(S));
        a = V(:,i);
        
        % other two principal axes
        a2 = V(:,2);
        a3 = V(:,3);        
        
        t = [-5, 5]';
        
        L = (t*0+1) * m + t * a';
        L2 = (t*0+1) * m + t * a2';
        L3 = (t*0+1) * m + t * a3';
       
        plot3( L(:,1), L(:,2), L(:,3),  'k-', 'LineWidth', 3)
        text( m(:,1)+0.3, m(:,2)+0.3, m(:,3)+0.3, [ num2str(uu),'-',mean( A{ii}(:,4), 1 )], 'Color', [ 0, 0.2, 0.2 ]  )
        set(gca,'DataAspectRatio', [1 1 1])
        xlabel('X')
        ylabel('Y')
        zlabel('Z')        
        grid on

        plot3( L2(:,1), L2(:,2), L2(:,3),  'r-', 'LineWidth', 3)
        plot3( L3(:,1), L3(:,2), L3(:,3),  'b-', 'LineWidth', 3)

        %         
        %     disp( a' )
        %     disp( v(ii,:))
        mReC{kk, 1}(ii, :) = [uu, kk, mean( A{ii}(:,4), 1 ), 4, m/10, a' ]; 
        disp( [ num2str( uu ), ': chain ', char( mean( A{ii}(:,4), 1 ) ), ' in file #', num2str(kk),', ', dfile  ]  )
        uu = uu + 1;
        
        normalvecs_cell{kk, 1}(ii,:) = a3';
        in_plane_vec2_cell{ kk, 1}(ii,:) = a2';   
        
    end 
end

% return
mRe = cell2mat( mReC );    
normals = cell2mat( normalvecs_cell );  
plane2 =  cell2mat( in_plane_vec2_cell );  

[dfile,IFPath2] = uiputfile( [IFPath2,'\',dfile(1:end-4),'_mRe_vecdata.txt'] ) 

% outpath = [IFPath2,'\',dfile(1:end-4),'_mRe_vecdata.txt']
if dfile ~= 0
    save( [IFPath2, dfile], 'mRe', '-ascii')
    save( [IFPath2, dfile(1:end-4),'_normals.txt'], 'normals', '-ascii')
    save( [IFPath2, dfile,'_perps.txt'], 'plane2', '-ascii')    
    
else
    disp( 'not saved ')
end
disp( 'end OK')