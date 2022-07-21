function read_pdb_1
% reads text data from *.pdb files
% and converts them to FASTA file 
% non-aminoacids are also labelled 
% and put to the end of thei respective protein chains
% ver. 1.01 2021-April-28 (c)DB

clear all; close all; 


% file to analyze:
fid = fopen( 'C:\DataVB\Progs\PBS_2021\current_path.txt' );
fp0 = fgetl(fid);
fclose(fid); 

disp( '##############################################################')
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

%    C = read_PDB_coord_02( [IFPath2,'\',dfile] , 1, reslabel );
    C = PDB_to_seq( [IFPath2,'\',dfile] , 1, reslabel );
    
    for jj = 1:size(C,1)
        disp( C{jj} )
    end
    
    disp( '...end of file ')
    
    button = questdlg('Do you want to save the sequence(s) into file?',...
        'Save?','Yes','No','Help','No');
    if strcmp(button,'Yes')
        disp('Creating file')
        [sfile, spath] =  uiputfile( [IFPath2,'\',dfile(1:end-4),'_seq.faa' ] );
        seq2file( C, sfile, spath );
    elseif strcmp(button,'No')
        disp(' ')
        disp(' ')
        disp('Canceled file operation')
    elseif strcmp(button,'Help')
        disp('Sorry, no help available')
    end

end
disp('end OK' )

%==========================================================
function XYZ = read_PDB_coord_02(FName, startL, reslabel)
% opens and analyses a .PDB file 
% version as of Apr 2021
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

%--------------------------------------------------------------------------
function seq = PDB_to_seq(FName, startL, reslabel)
% opens and analyses a .PDB file 
% version as of Oct 2015
%

code3 = ({ 'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO', 'ASP', 'GLU', 'GLY', 'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'CYC', '45D' })';
numat =    [ 5,7,8,8,11,14,8,7,8,9,4,6,7,6,12,8,9,9,11,10,43,94 ];
code1 = ({ 'A',   'V',   'L',   'I',   'F',   'W'  , 'M',   'P',   'D',   'E',   'G',   'S',   'T',   'C',   'Y',   'N',   'Q',   'K',   'R',   'H', 'B', 'X'  })'; 

% lines from a sample PDB file
%ATOM      1  N   SER A   2     303.103 337.644 327.194  1.00158.19      AA   N
%ATOM  28976  N   MET A   1     302.131 337.179 315.160  1.00167.85      BA   N

fid = fopen(FName);
disp(['extracting coordinates from file ', FName, ', please wait...'])
a = textread( FName,'%s','delimiter', '\n');
fclose(fid);

XYZ = [];
ii = 1;
jj = 1;

curchain = '11';
chcount = 1;
c = 1;
r = 0;
acount = 0;
count2match = 1e5;

for ii = startL:size(a,1)-1
        %     disp( a{ii} )
        if ~isempty( findstr( a{ii}(1:7) , 'ATOM') ) | ~isempty( findstr( a{ii}(1:7), 'HETATM') )
%             disp( a{ii} )
            AA = a{ii}(18:20);
            chain = a{ii}(22);            
            chain2 = a{ii}(73:74);          
            
            h = strcmp(curchain , chain2);           
            if ~h
                curchain = chain2;
                curchain1 = chain;
                
                r = r + 1;
                c = 1;
                seq{r,1} = ['>pdb|chain ',curchain1,'_',num2str( chcount )];
                chcount = chcount + 1;
                r = r + 1;
                curAA = 'XXX';                
            end

            k = strcmp(AA , curAA);           
            if ~k
               acount = 0;
               curAA = AA;
               x = strmatch(curAA, code3 );             
%                if ~isempty(x)
                   seq{r,1}(c) = code1{x};
%                else
%                    seq{r,1}(c) = 'X';
%                end
               count2match = numat( x );
               c = c+1;
               if c == 70
                   c = 1;
                   r = r + 1;
               end
           else
               acount = acount  + 1;             
           end
           if k & acount > count2match
               acount = 0;
               curAA = AA;
               x = strmatch(curAA, code3 );             
%                if ~isempty(x)
                   seq{r,1}(c) = code1{x};
%                else
%                    seq{r,1}(c) = 'X';
%                end
               count2match = numat(x);
               c = c+1;
               if c == 70
                   c = 1;
                   r = r + 1;
               end
           end
            
%             poslab = findstr( a{ii}, reslabel );
%             chainlab = a{ii}(poslab+4);
%             atomindexinchain = a{ii}(poslab+5:30);
%             
%             XYZpart1 = a{ii}(31:38);
%             XYZpart2 = a{ii}(39:47);
%             XYZpart3 = a{ii}(48:56);        
%             Coord1 = strread(XYZpart1,'%s');
%             Coord2 = strread(XYZpart2,'%s');        
%             Coord3 = strread(XYZpart3,'%s');        
%             if ~isempty(Coord1)
%                 XYZ(jj,1) = str2num(char(Coord1));
%                 XYZ(jj,2) = str2num(char(Coord2));
%                 XYZ(jj,3) = str2num(char(Coord3));    
%                 XYZ(jj,4) = double( chainlab);        
%                 XYZ(jj,5) = str2num( atomindexinchain );        
%                 
%                 jj = jj + 1;
%             end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end
        
    end



disp('...done')

%------------------------------------------------------------
function f = seq2file( L, fname, fpath )
%
%
%
fid = fopen([fpath,'\',fname],'w');
for ii = 1:size(L,1)
    fprintf(fid, '%s \r\n', L{ii,1});
end

fclose(fid);
f = 0;

%-------------------------------------------
