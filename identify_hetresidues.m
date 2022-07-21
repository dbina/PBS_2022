function read_pdb_1
% reads sequence data from a FASTA file
% uses that information to identify the pigments
% using c-terminal bit of the protein preceeding pigment in the sequence
%

clear all; close all; 

% file to analyze:
fid = fopen( 'C:\DataVB\Progs\PBS_2021\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid);    


% load the data from the dictionary file 
[cterm, index, bindsite] = read_dict_file( [fp0,'dictionary_v0.txt'] )

% load binding sites
[bf, bp] = uigetfile( [fp0, '*.*'], 'load binding sites' )
sites = load( [bp, bf] );

% return
disp( '##############################################################')
disp( 'starting tres_reader, selecting the tres data folder')
% get the directory and read data
directory_name = uigetdir( fp0 )
% make the list of files
d = dir( directory_name );
dsize = max( size( {d.name} ) );
str0 = {d.name};

jj = 1;
% for ii = 1:dsize
%      
%    if ~isempty( findstr( str0{ii}, '*.faa' ) )
%        str{jj} = str0{ii};
%        jj = jj + 1;
%    end
% 
% end
[s,v] = listdlg('PromptString','Select a file:',...
                 'SelectionMode','multiple',...
                 'ListString',str0);
    
IFPath2 = directory_name;

for ii = 1:length( s ) 
    
    CM = [];  
    dfile = str0{s(ii)};

    C = faread( [IFPath2,'\',dfile] );
      
    size(C)
    kk = 1;
    for jj = 1:size(C,1)
        for uu = 1:size( C,2)
           disp( C{jj,uu} )
           
           faline = C{jj,uu};
           k = findstr(faline, 'B');
           if ~isempty( k ) & uu > 2
               
               for nn = 1:length( k )
                   
                   newline{kk,1} = [ num2str(kk),': ', faline( k(nn) ),' ', faline( k(1) - 9 : k(1)-1 ),' ', C{jj,1}(6:end) ];
                   kk = kk + 1;
               end
                            
           end
           
        end
    end
    % write the list to file
    % and prepare the indexes for pigments to assign spectra
    fid = fopen(  [IFPath2,'\',dfile(1:end-4),'_cterm.txt' ], 'w' );
    disp('****')
    for kk = 1:size( newline,1 )
        disp( newline{kk} )
        textp = char( newline{kk} );
        fprintf( fid, '%s', textp);
        fprintf( fid, '\n');
    end
    fclose( fid );
    disp( '...end of file ')

    if size( sites, 1 ) ~= size( newline,1 )
        disp( 'mismatch between sequence and site data ' )
        return
    end
    
    for kk = 1:size( newline,1 )
        found = 0; 
        for nn = 1:size( index, 1) 
            
            u = findstr( newline{kk}, cterm{nn} );
            %disp( [ newline{kk}, '   ', cterm{nn}, '  ', num2str( u ) ] );
            if ~isempty( u ) & bindsite(nn) == sites( kk )

                code(kk,1) = index( nn );
                found = 1;                
                disp( [ newline{kk}, '   ', cterm{nn}, ' at binding site #',...
                        num2str( bindsite(nn) ), ' identified as ', num2str( index( nn ) ) ] );  
                
            end
            
        end
        
        if ~found
           code(kk,1) = 0;            
           disp( [ newline{kk}, '   ', cterm{nn}, ' not identified' ] );       
        end
        
    end
    save( [IFPath2,'\',dfile(1:end-4),'_index.txt' ] , 'code' ,'-ascii')
%     disp( code)
end


disp('end OK' )

%==========================================================
function I = faread( FileName )
%
%
%
a = textread( FileName,'%s','delimiter', '\n');

jj = 1;
kk = 0;
ii = 1;
while 1
    if jj > size(a,1)
       disp('...EOF')
       break
    end     

    if length( a{jj}) ~= 0        
        if a{jj}(1) == '>' 
            kk = kk + 1;
            ii = 1;
        end
        
        if kk > 0
            I{kk,ii} = a{jj};
            ii = ii + 1;           
        end
    end
    jj = jj + 1;
    
end
disp(FileName)
disp([num2str(size(a,1)), ' lines, ', num2str(kk), ' sequences read'])       

%--------------------------------------------------------------------------
function [cterm, index, bindsite] = read_dict_file( dicfile )
%
%
%
a = textread( dicfile,'%s','delimiter', '\n');

jj = 1;
kk = 0;
ii = 1;
while 1
     
    if jj > size(a,1)
       disp('...EOF')
       break
    end     
    textline = char( a{jj} );
    if isempty( strfind( textline(1:4), '%' ) )
        disp( textline  )
        s = findstr( textline, '=' );
        cterm{ii,1} = textline( 1:min(s)-1 );        
        index(ii,1) = str2num( textline( max(s)+1:end ) );        
        bindsite(ii,1) = str2num( textline( s(2)+1:s(3)-1 ) );             
        ii = ii + 1;
    end
    
    jj = jj + 1;
    
end
% bindsite =  index;