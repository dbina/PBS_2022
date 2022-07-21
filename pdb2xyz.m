function read_pdb_1
% reads text data from *.pdb files
% the residues are identified by the variable 'reslabel'
% output as of 2021-04-28 consists of x-y-z coordinate latter label (for
% chain) and number index of atom from original *.pdb file
% ver. 1.1 2021-May-17 (c)DB
% last modified in: GhostXP

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
     
   if ~isempty( findstr( str0{ii}, '.pdb' ) )
       str{jj} = str0{ii};
       jj = jj + 1;
   end

end
[s,v] = listdlg('PromptString','Select a file:',...
                 'SelectionMode','multiple',...
                 'ListString',str);

reslabel = 'CYC';           
IFPath2 = directory_name;

for ii = 1:length( s ) 
    
    CM = [];  
    dfile = str{s(ii)};

    C = read_PDB_coord_02( [IFPath2,'\',dfile] , 1, reslabel );
    
    if ~isempty( C )
        na = length(C(:,1))
        if length( s) == 1
           [dfile,IFPath2] = uiputfile([IFPath2,'\',dfile(1:end-4),'_', reslabel, '.xyz']) 
           if dfile == 0
               break
           else
               save( [IFPath2,'\',dfile], 'C', '-ascii') 
           end
        else
           save( [IFPath2,'\',dfile(1:end-4),'_', reslabel, '.xyz'], 'C', '-ascii')
        end
    else
        disp( ['this file does not contain residue ',reslabel,', nothing saved' ]  )
    end
    
end

disp('end OK' )
%==========================================================
function XYZ = read_PDB_coord_02(FName, startL, reslabel)
% opens and analyses a .PDB file 
% version as of Oct 2015
%
fid = fopen(FName);
disp(['extracting coordinates from file ', FName, ', please wait...'])
a = textread( FName,'%s','delimiter', '\n');
fclose(fid);

XYZ = [];
ii = 1;
jj = 1;
~strcmp('all', reslabel )
if ~strcmp('all', reslabel )
    for ii = startL:size(a,1)-1
        %     disp( a{ii} )
        %if ~isempty( findstr( a{ii}(1:7) , 'ATOM') ) | ~isempty( findstr( a{ii}(1:7), 'HETATM') )
        if (~isempty( findstr( a{ii} , 'HETATM') )|~isempty( findstr( a{ii}(1:7) , 'ATOM') )) & (~isempty( findstr( a{ii}, reslabel ) ) & findstr( a{ii}, reslabel )>16) 
            disp( a{ii} )
            poslab = findstr( a{ii}, reslabel );
            chainlab = a{ii}(poslab+4);
            atomindexinchain = a{ii}(poslab+5:30);
            
            XYZpart1 = a{ii}(31:38);
            XYZpart2 = a{ii}(39:47);
            XYZpart3 = a{ii}(48:56);        
            Coord1 = strread(XYZpart1,'%s');
            Coord2 = strread(XYZpart2,'%s');        
            Coord3 = strread(XYZpart3,'%s');        
            if ~isempty(Coord1)
                XYZ(jj,1) = str2num(char(Coord1));
                XYZ(jj,2) = str2num(char(Coord2));
                XYZ(jj,3) = str2num(char(Coord3));    
                XYZ(jj,4) = double( chainlab);        
                XYZ(jj,5) = str2num( atomindexinchain );        
                
                jj = jj + 1;
            end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end
        
    end
else
    for ii = startL:size(a,1)-1
        %     disp( a{ii} )
        if ~isempty( findstr( a{ii}(1:7) , 'ATOM') ) | ~isempty( findstr( a{ii}(1:7), 'HETATM') )
%         if ~isempty( findstr( a{ii} , 'HETATM') ) & ~isempty( findstr( a{ii}, reslabel ) )    
%               disp( a{ii} )
%             poslab = findstr( a{ii}, reslabel );
            chainlab = a{ii}(22);
            
            XYZpart1 = a{ii}(31:38);
            XYZpart2 = a{ii}(39:47);
            XYZpart3 = a{ii}(48:56);        
            Coord1 = strread(XYZpart1,'%s');
            Coord2 = strread(XYZpart2,'%s');        
            Coord3 = strread(XYZpart3,'%s');        
            if ~isempty(Coord1)
                XYZ(jj,1) = str2num(char(Coord1));
                XYZ(jj,2) = str2num(char(Coord2));
                XYZ(jj,3) = str2num(char(Coord3));    
                XYZ(jj,4) = double( chainlab);        
                
                jj = jj + 1;
            end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end
    end    
end


disp('...done')

%------------------------------------------------------------
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
