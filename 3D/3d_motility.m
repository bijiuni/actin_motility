function k=threedactin(varargin)

narginchk(0,1);
if (nargin == 0),
     [filename, pathname] = uigetfile( ...
	       {'*.jpg;*.tif;*.gif;*.png;*.bmp', ...
		'All MATLAB Image Files (*.jpg,*.tif,*.gif,*.png,*.bmp)'; ...
		'*.jpg;*.jpeg', ...
		'TIFF Files (*.tif,*.tiff)'; ...
		'*.gif', ...
		'GIF Files (*.gif)'; ...
		'*.png', ...
		'PNG Files (*.png)'; ...
		'*.bmp', ...
		'Bitmap Files (*.bmp)'; ...
		'*.*', ...
		'All Files (*.*)'}, ...
	       'Select image file');
     if isequal(filename,0) || isequal(pathname,0)
	  return
     else
	  imagename = fullfile(pathname, filename);
     end
elseif nargin == 1,
     imagename = varargin{1};
     [~, file,ext] = fileparts(imagename);
     filename = strcat(file,ext);
end

FileTif=imagename;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
 
FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
for i=1:NumberImages
   FinalImage(:,:,i)=imread(FileTif,'Index',i);
end

testvol = logical(FinalImage);

skel = Skeleton3D(testvol);

% Prompt user for threshold value
prompt={'Enter the shortest brand threshold value'};
def={'0'};
dlgTitle='THRE: user input required';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
if (isempty(char(answer{:})) == 1),
     close(FigName)
     return
else
    thre = str2num(char(answer{:}));
end

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node,link] = Skel2Graph3D(skel,thre);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,0);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)

    wl = wl_new;   
    
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,0);

     wl_new = sum(cellfun('length',{node2.links}));

end;

skel2 = smooth3 (skel2);

figure();
col=[.7 .7 .8];
hiso = patch(isosurface(testvol,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(testvol,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(testvol,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel2,1);
l=size(skel2,2);
h=size(skel2,3);
[x,y,z]=ind2sub([w,l,h],find(skel2(:)));
plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)

points = [x(:),y(:),z(:)];

%find the normals and curvature
[k] = findPointNormals(points,[],[0,0,10],true);

end

function EulerInv =  p_EulerInv(img,LUT)

% Calculate Euler characteristic for each octant and sum up
eulerChar = zeros(size(img,1),1);
% Octant SWU
n = ones(size(img,1),1);
n(img(:,25)==1) = bitor(n(img(:,25)==1),128);
n(img(:,26)==1) = bitor(n(img(:,26)==1),64);
n(img(:,16)==1) = bitor(n(img(:,16)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,22)==1) = bitor(n(img(:,22)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,13)==1) = bitor(n(img(:,13)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SEU
n = ones(size(img,1),1);
n(img(:,27)==1) = bitor(n(img(:,27)==1),128);
n(img(:,24)==1) = bitor(n(img(:,24)==1),64);
n(img(:,18)==1) = bitor(n(img(:,18)==1),32);
n(img(:,15)==1) = bitor(n(img(:,15)==1),16);
n(img(:,26)==1) = bitor(n(img(:,26)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,17)==1) = bitor(n(img(:,17)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NWU
n = ones(size(img,1),1);
n(img(:,19)==1) = bitor(n(img(:,19)==1),128);
n(img(:,22)==1) = bitor(n(img(:,22)==1),64);
n(img(:,10)==1) = bitor(n(img(:,10)==1),32);
n(img(:,13)==1) = bitor(n(img(:,13)==1),16);
n(img(:,20)==1) = bitor(n(img(:,20)==1),8);
n(img(:,23)==1) = bitor(n(img(:,23)==1),4);
n(img(:,11)==1) = bitor(n(img(:,11)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NEU
n = ones(size(img,1),1);
n(img(:,21)==1) = bitor(n(img(:,21)==1),128);
n(img(:,24)==1) = bitor(n(img(:,24)==1),64);
n(img(:,20)==1) = bitor(n(img(:,20)==1),32);
n(img(:,23)==1) = bitor(n(img(:,23)==1),16);
n(img(:,12)==1) = bitor(n(img(:,12)==1),8);
n(img(:,15)==1) = bitor(n(img(:,15)==1),4);
n(img(:,11)==1) = bitor(n(img(:,11)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SWB
n = ones(size(img,1),1);
n(img(:,7)==1) = bitor(n(img(:,7)==1),128);
n(img(:,16)==1) = bitor(n(img(:,16)==1),64);
n(img(:,8)==1) = bitor(n(img(:,8)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,4)==1) = bitor(n(img(:,4)==1),8);
n(img(:,13)==1) = bitor(n(img(:,13)==1),4);
n(img(:,5)==1) = bitor(n(img(:,5)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant SEB
n = ones(size(img,1),1);
n(img(:,9)==1) = bitor(n(img(:,9)==1),128);
n(img(:,8)==1) = bitor(n(img(:,8)==1),64);
n(img(:,18)==1) = bitor(n(img(:,18)==1),32);
n(img(:,17)==1) = bitor(n(img(:,17)==1),16);
n(img(:,6)==1) = bitor(n(img(:,6)==1),8);
n(img(:,5)==1) = bitor(n(img(:,5)==1),4);
n(img(:,15)==1) = bitor(n(img(:,15)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NWB
n = ones(size(img,1),1);
n(img(:,1)==1) = bitor(n(img(:,1)==1),128);
n(img(:,10)==1) = bitor(n(img(:,10)==1),64);
n(img(:,4)==1) = bitor(n(img(:,4)==1),32);
n(img(:,13)==1) = bitor(n(img(:,13)==1),16);
n(img(:,2)==1) = bitor(n(img(:,2)==1),8);
n(img(:,11)==1) = bitor(n(img(:,11)==1),4);
n(img(:,5)==1) = bitor(n(img(:,5)==1),2);
eulerChar = eulerChar + LUT(n);
% Octant NEB
n = ones(size(img,1),1);
n(img(:,3)==1) = bitor(n(img(:,3)==1),128);
n(img(:,2)==1) = bitor(n(img(:,2)==1),64);
n(img(:,12)==1) = bitor(n(img(:,12)==1),32);
n(img(:,11)==1) = bitor(n(img(:,11)==1),16);
n(img(:,6)==1) = bitor(n(img(:,6)==1),8);
n(img(:,5)==1) = bitor(n(img(:,5)==1),4);
n(img(:,15)==1) = bitor(n(img(:,15)==1),2);
eulerChar = eulerChar + LUT(n);

EulerInv(eulerChar==0) = true;
end

function p_is_simple = p_is_simple(N)

% copy neighbors for labeling
n_p = size(N,1);
p_is_simple = ones(1,n_p);

cube = zeros(26,n_p);
cube(1:13,:)=N(:,1:13)';
cube(14:26,:)=N(:,15:27)';

label = 2*ones(1,n_p);

% for all points in the neighborhood
for i=1:26
    
    idx_1 = find(cube(i,:)==1);
    idx_2 = find(p_is_simple);
    idx = intersect(idx_1,idx_2);
    
    if(~isempty(idx))
        
        % start recursion with any octant that contains the point i
        switch( i )
            
            case {1,2,4,5,10,11,13}
                cube(:,idx) = p_oct_label(1, label, cube(:,idx) );
            case {3,6,12,14}
                cube(:,idx) = p_oct_label(2, label, cube(:,idx) );
            case {7,8,15,16}
                cube(:,idx) = p_oct_label(3, label, cube(:,idx) );
            case {9,17}
                cube(:,idx) = p_oct_label(4, label, cube(:,idx) );
            case {18,19,21,22}
                cube(:,idx) = p_oct_label(5, label, cube(:,idx) );
            case {20,23}
                cube(:,idx) = p_oct_label(6, label, cube(:,idx) );
            case {24,25}
                cube(:,idx) = p_oct_label(7, label, cube(:,idx) );
            case 26,
                cube(:,idx) = p_oct_label(8, label, cube(:,idx) );
        end;

        label(idx) = label(idx)+1;
        del_idx = find(label>=4);
        
        if(~isempty(del_idx))
            p_is_simple(del_idx) = 0;
        end;
    end;
end;
end

function cube = p_oct_label(octant, label, cube)

% check if there are points in the octant with value 1
if( octant==1 )
    
    % set points in this octant to current label
    % and recurseive labeling of adjacent octants
    idx_1 = find(cube(1,:) == 1);
    if(~isempty(idx_1))
        cube(1,idx_1) = label(idx_1);
    end;
    
    idx_2 = find(cube(2,:) == 1);
    if(~isempty(idx_2))
        cube(2,idx_2) = label(idx_2);
        cube(:,idx_2) = p_oct_label(2,label(idx_2),cube(:,idx_2));
    end;
    
    idx_4 = find(cube(4,:) == 1);
    if(~isempty(idx_4))
        cube(4,idx_4) = label(idx_4);
        cube(:,idx_4) = p_oct_label(3,label(idx_4),cube(:,idx_4));
    end;
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end;
    
    idx_10 = find(cube(10,:) == 1);
    if(~isempty(idx_10))
        cube(10,idx_10) = label(idx_10);
        cube(:,idx_10) = p_oct_label(5,label(idx_10),cube(:,idx_10));
    end;
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end;
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end;
    
end;

if( octant==2 )
    
    idx_2 = find(cube(2,:) == 1);
    if(~isempty(idx_2))
        cube(2,idx_2) = label(idx_2);
        cube(:,idx_2) = p_oct_label(1,label(idx_2),cube(:,idx_2));
    end;

    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end;

    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end;

    idx_3 = find(cube(3,:) == 1);
    if(~isempty(idx_3))
        cube(3,idx_3) = label(idx_3);
    end;

    idx_6 = find(cube(6,:) == 1);
    if(~isempty(idx_6))
        cube(6,idx_6) = label(idx_6);
        cube(:,idx_6) = p_oct_label(4,label(idx_6),cube(:,idx_6));
    end;
    
    idx_12 = find(cube(12,:) == 1);
    if(~isempty(idx_12))
        cube(12,idx_12) = label(idx_12);
        cube(:,idx_12) = p_oct_label(6,label(idx_12),cube(:,idx_12));
    end;

    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end;

end;

if( octant==3 )
    
    idx_4 = find(cube(4,:) == 1);
    if(~isempty(idx_4))
        cube(4,idx_4) = label(idx_4);
        cube(:,idx_4) = p_oct_label(1,label(idx_4),cube(:,idx_4));
    end;

    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(4,label(idx_5),cube(:,idx_5));
    end;

    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end;

    idx_7 = find(cube(7,:) == 1);
    if(~isempty(idx_7))
        cube(7,idx_7) = label(idx_7);
    end;

    idx_8 = find(cube(8,:) == 1);
    if(~isempty(idx_8))
        cube(8,idx_8) = label(idx_8);
        cube(:,idx_8) = p_oct_label(4,label(idx_8),cube(:,idx_8));
    end;
    
    idx_15 = find(cube(15,:) == 1);
    if(~isempty(idx_15))
        cube(15,idx_15) = label(idx_15);
        cube(:,idx_15) = p_oct_label(7,label(idx_15),cube(:,idx_15));
    end;

    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_13))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end;
    
end;

if( octant==4 )
    
    idx_5 = find(cube(5,:) == 1);
    if(~isempty(idx_5))
        cube(5,idx_5) = label(idx_5);
        cube(:,idx_5) = p_oct_label(1,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(2,label(idx_5),cube(:,idx_5));
        cube(:,idx_5) = p_oct_label(3,label(idx_5),cube(:,idx_5));
    end;

    idx_6 = find(cube(6,:) == 1);
    if(~isempty(idx_6))
        cube(6,idx_6) = label(idx_6);
        cube(:,idx_6) = p_oct_label(2,label(idx_6),cube(:,idx_6));
    end;

    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end;
    
    idx_8 = find(cube(8,:) == 1);
    if(~isempty(idx_8))
        cube(8,idx_8) = label(idx_8);
        cube(:,idx_8) = p_oct_label(3,label(idx_8),cube(:,idx_8));
    end;

    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end;

    idx_9 = find(cube(9,:) == 1);
    if(~isempty(idx_9))
        cube(9,idx_9) = label(idx_9);
    end;

    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(8,label(idx_17),cube(:,idx_17));
    end;

end;

if( octant==5 )
    
    idx_10 = find(cube(10,:) == 1);
    if(~isempty(idx_10))
        cube(10,idx_10) = label(idx_10);
        cube(:,idx_10) = p_oct_label(1,label(idx_10),cube(:,idx_10));
    end;

    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(6,label(idx_11),cube(:,idx_11));
    end;
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(7,label(idx_13),cube(:,idx_13));
    end;

    idx_18 = find(cube(18,:) == 1);
    if(~isempty(idx_18))
        cube(18,idx_18) = label(idx_18);
    end;

    idx_19 = find(cube(19,:) == 1);
    if(~isempty(idx_19))
        cube(19,idx_19) = label(idx_19);
        cube(:,idx_19) = p_oct_label(6,label(idx_19),cube(:,idx_19));
    end;

    idx_21 = find(cube(21,:) == 1);
    if(~isempty(idx_21))
        cube(21,idx_21) = label(idx_21);
        cube(:,idx_21) = p_oct_label(7,label(idx_21),cube(:,idx_21));
    end;

    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end;

end;

if( octant==6 )
    
    idx_11 = find(cube(11,:) == 1);
    if(~isempty(idx_11))
        cube(11,idx_11) = label(idx_11);
        cube(:,idx_11) = p_oct_label(1,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(2,label(idx_11),cube(:,idx_11));
        cube(:,idx_11) = p_oct_label(5,label(idx_11),cube(:,idx_11));
    end;

    idx_12 = find(cube(12,:) == 1);
    if(~isempty(idx_12))
        cube(12,idx_12) = label(idx_12);
        cube(:,idx_12) = p_oct_label(2,label(idx_12),cube(:,idx_12));
    end;

    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(8,label(idx_14),cube(:,idx_14));
    end;
    
    idx_19 = find(cube(19,:) == 1);
    if(~isempty(idx_19))
        cube(19,idx_19) = label(idx_19);
        cube(:,idx_19) = p_oct_label(5,label(idx_19),cube(:,idx_19));
    end;


    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end;
    
    idx_20 = find(cube(20,:) == 1);
    if(~isempty(idx_20))
        cube(20,idx_20) = label(idx_20);
    end;

    idx_23 = find(cube(23,:) == 1);
    if(~isempty(idx_23))
        cube(23,idx_23) = label(idx_23);
        cube(:,idx_23) = p_oct_label(8,label(idx_23),cube(:,idx_23));
    end;
 
end;

if( octant==7 )
    
    idx_13 = find(cube(13,:) == 1);
    if(~isempty(idx_13))
        cube(13,idx_13) = label(idx_13);
        cube(:,idx_13) = p_oct_label(1,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(3,label(idx_13),cube(:,idx_13));
        cube(:,idx_13) = p_oct_label(5,label(idx_13),cube(:,idx_13));
    end;

    idx_15 = find(cube(15,:) == 1);
    if(~isempty(idx_15))
        cube(15,idx_15) = label(idx_15);
        cube(:,idx_15) = p_oct_label(3,label(idx_15),cube(:,idx_15));
    end;

    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(8,label(idx_16),cube(:,idx_16));
    end;

    idx_21 = find(cube(21,:) == 1);
    if(~isempty(idx_21))
        cube(21,idx_21) = label(idx_21);
        cube(:,idx_21) = p_oct_label(5,label(idx_21),cube(:,idx_21));
    end;

    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(8,label(idx_22),cube(:,idx_22));
    end;

    idx_24 = find(cube(24,:) == 1);
    if(~isempty(idx_24))
        cube(24,idx_24) = label(idx_24);
    end;
    
    idx_25 = find(cube(25,:) == 1);
    if(~isempty(idx_25))
        cube(25,idx_25) = label(idx_25);
        cube(:,idx_25) = p_oct_label(8,label(idx_25),cube(:,idx_25));
    end;
end;

if( octant==8 )
    
    idx_14 = find(cube(14,:) == 1);
    if(~isempty(idx_14))
        cube(14,idx_14) = label(idx_14);
        cube(:,idx_14) = p_oct_label(2,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(4,label(idx_14),cube(:,idx_14));
        cube(:,idx_14) = p_oct_label(6,label(idx_14),cube(:,idx_14));
    end;

    idx_16 = find(cube(16,:) == 1);
    if(~isempty(idx_16))
        cube(16,idx_16) = label(idx_16);
        cube(:,idx_16) = p_oct_label(3,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(4,label(idx_16),cube(:,idx_16));
        cube(:,idx_16) = p_oct_label(7,label(idx_16),cube(:,idx_16));
    end;
    
    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(4,label(idx_17),cube(:,idx_17));
    end;
    
    idx_22 = find(cube(22,:) == 1);
    if(~isempty(idx_22))
        cube(22,idx_22) = label(idx_22);
        cube(:,idx_22) = p_oct_label(5,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(6,label(idx_22),cube(:,idx_22));
        cube(:,idx_22) = p_oct_label(7,label(idx_22),cube(:,idx_22));
    end;
    
    idx_17 = find(cube(17,:) == 1);
    if(~isempty(idx_17))
        cube(17,idx_17) = label(idx_17);
        cube(:,idx_17) = p_oct_label(4,label(idx_17),cube(:,idx_17));
    end;
    
    idx_23 = find(cube(23,:) == 1);
    if(~isempty(idx_23))
        cube(23,idx_23) = label(idx_23);
        cube(:,idx_23) = p_oct_label(6,label(idx_23),cube(:,idx_23));
    end;
    
    idx_25 = find(cube(25,:) == 1);
    if(~isempty(idx_25))
        cube(25,idx_25) = label(idx_25);
        cube(:,idx_25) = p_oct_label(7,label(idx_25),cube(:,idx_25));
    end;
    
    idx_26 = find(cube(26,:) == 1);
    if(~isempty(idx_26))
        cube(26,idx_26) = label(idx_26);
    end;
end;

end

function [vox,n_idx,ep] = pk_follow_link(skel,node,k,j,idx,cans,c2n)

vox = [];
n_idx = [];
ep = 0;

% assign start node to first voxel
vox(1) = node(k).idx(j);

i=1;
isdone = false;
while(~isdone) % while no node reached
    i=i+1; % next voxel
    next_cand = c2n(idx);
        cand = cans(next_cand,2);
        if(cand==vox(i-1)) % switch direction
            cand = cans(next_cand,3);
        end;
        if(skel(cand)>1) % node found
            vox(i) = idx;
            vox(i+1) = cand; % first node
            n_idx = skel(cand)-1; % node #
            if(node(n_idx).ep)
                ep=1;
            end;
            isdone = 1;
        else % next voxel
            vox(i) = idx;
            idx = cand;
        end;
end;
end

function nhood = pk_get_nh(img,i)

width = size(img,1);
height = size(img,2);
depth = size(img,3);

[x,y,z]=ind2sub([width height depth],i);

nhood = false(length(i),27);

for xx=1:3
    for yy=1:3
        for zz=1:3
            w=sub2ind([3 3 3],xx,yy,zz);
            idx = sub2ind([width height depth],x+xx-2,y+yy-2,z+zz-2);
            nhood(:,w)=img(idx);
        end;
    end;
end;
end

function nhood = pk_get_nh_idx(img,i)

width = size(img,1);
height = size(img,2);
depth = size(img,3);

[x,y,z]=ind2sub([width height depth],i);

nhood = zeros(length(i),27);

for xx=1:3
    for yy=1:3
        for zz=1:3
            w=sub2ind([3 3 3],xx,yy,zz);
            nhood(:,w) = sub2ind([width height depth],x+xx-2,y+yy-2,z+zz-2);
        end;
    end;
end;
end


function [A,node,link] = Skel2Graph3D(skel,THR)
% SKEL2GRAPH3D Calculate the network graph of a 3D voxel skeleton
%
% [A,node,link] = SKEL2GRAPH3D(skel,THR)
%
% where "skel" is the input 3D binary image, and "THR" is a threshold for 
% the minimum length of branches. A is the adjacency matrix, and node/link
% are structures describing node and link properties
%
% Philip Kollmannsberger (philipk@gmx.net)
%
% For more information, see <a
% href="matlab:web('http://uk.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d')">Skel2Graph3D</a> at the MATLAB File Exchange.

% pad volume with zeros
skel=padarray(skel,[1 1 1]);

% image dimensions
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);

% need this for labeling nodes etc.
skel2 = uint16(skel);

% all foreground voxels
list_canal=find(skel);

% 26-nh of all canal voxels
nh = logical(pk_get_nh(skel,list_canal));

% 26-nh indices of all canal voxels
nhi = pk_get_nh_idx(skel,list_canal);

% # of 26-nb of each skel voxel + 1
sum_nh = sum(logical(nh),2);

% all canal voxels with >2 nb are nodes
nodes = list_canal(sum_nh>3);

% all canal voxels with exactly one nb are end nodes
ep = list_canal(sum_nh==2);

% all canal voxels with exactly 2 nb
cans = list_canal(sum_nh==3);

% Nx3 matrix with the 2 nb of each canal voxel
can_nh_idx = pk_get_nh_idx(skel,cans);
can_nh = pk_get_nh(skel,cans);

% remove center of 3x3 cube
can_nh_idx(:,14)=[];
can_nh(:,14)=[];

% keep only the two existing foreground voxels
can_nb = sort(logical(can_nh).*can_nh_idx,2);

% remove zeros
can_nb(:,1:end-2) = [];

% add neighbours to canalicular voxel list (this might include nodes)
cans = [cans can_nb];

% group clusters of node voxels to nodes
node=[];
link=[];

tmp=false(w,l,h);
tmp(nodes)=1;
cc2=bwconncomp(tmp); % number of unique nodes
num_realnodes = cc2.NumObjects;

% create node structure
for i=1:cc2.NumObjects
    node(i).idx = cc2.PixelIdxList{i};
    node(i).links = [];
    node(i).conn = [];
    [x,y,z]=ind2sub([w l h],node(i).idx);
    node(i).comx = mean(x);
    node(i).comy = mean(y);
    node(i).comz = mean(z);
    node(i).ep = 0;
    
    % assign index to node voxels
    skel2(node(i).idx) = i+1;
end;

tmp=false(w,l,h);
tmp(ep)=1;
cc3=bwconncomp(tmp); % number of unique nodes

% create node structure
for i=1:cc3.NumObjects
    ni = num_realnodes+i;
    node(ni).idx = cc3.PixelIdxList{i};
    node(ni).links = [];
    node(ni).conn = [];
    [x,y,z]=ind2sub([w l h],node(ni).idx);
    node(ni).comx = mean(x);
    node(ni).comy = mean(y);
    node(ni).comz = mean(z);
    node(ni).ep = 1;
    
    % assign index to node voxels
    skel2(node(ni).idx) = ni+1;
end;

l_idx = 1;

c2n=zeros(w*l*h,1);
c2n(cans(:,1))=1:length(cans);

s2n=zeros(w*l*h,1);
s2n(nhi(:,14))=1:length(nhi);

% visit all nodes
for i=1:num_realnodes

    % find all canal vox in nb of all node idx
    link_idx = s2n(node(i).idx);
    
    for j=1:length(link_idx)
        % visit all voxels of this node
        
        % all potential unvisited links emanating from this voxel
        link_cands = nhi(link_idx(j),nh(link_idx(j),:)==1);
        link_cands = link_cands(skel2(link_cands)==1);
        
        for k=1:length(link_cands)
            [vox,n_idx,ep] = pk_follow_link(skel2,node,i,j,link_cands(k),cans,c2n);
            skel2(vox(2:end-1))=0;
            if((ep && length(vox)>THR) || (~ep && i~=n_idx))
                link(l_idx).n1 = i;
                link(l_idx).n2 = n_idx; % node number
                link(l_idx).point = vox;
                node(i).links = [node(i).links, l_idx];
                node(i).conn = [int16(node(i).conn), int16(n_idx)];
                node(n_idx).links = [node(n_idx).links, l_idx];
                node(n_idx).conn = [int16(node(n_idx).conn), int16(i)];
                l_idx = l_idx + 1;
            end;
        end;
    end;
        
end;

% mark all 1-nodes as end points
ep_idx = find(cellfun('length',{node.links})==1);
for i=1:length(ep_idx)
    node(ep_idx(i)).ep = 1;    
end;

% number of nodes
n_nodes = length(node);

% initialize matrix
A = zeros(n_nodes);

% for all nodes, make according entries into matrix for all its links
for i=1:n_nodes
    idx1=find(node(i).conn>0);
    idx2=find(node(i).links>0);
    idx=intersect(idx1,idx2);
    for j=1:length(idx) % for all its links
        if(i==link(node(i).links(idx(j))).n1) % if we are the starting point
            A(i,link(node(i).links(idx(j))).n2)=length(link(node(i).links(idx(j))).point);
            A(link(node(i).links(idx(j))).n2,i)=length(link(node(i).links(idx(j))).point);
        end;
        if(i==link(node(i).links(idx(j))).n2) % if we are the end point
            A(i,link(node(i).links(idx(j))).n1)=length(link(node(i).links(idx(j))).point);
            A(link(node(i).links(idx(j))).n1,i)=length(link(node(i).links(idx(j))).point);
        end;
    end;
end;

% convert to sparse
A = sparse(A);

% transform all voxel and position indices back to non-padded coordinates
for i=1:length(node)
    [x,y,z] = ind2sub([w,l,h],node(i).idx);
    node(i).idx = sub2ind([w-2,l-2,h-2],x-1,y-1,z-1);
    node(i).comx = node(i).comx - 1;
    node(i).comy = node(i).comy - 1;
    node(i).comz = node(i).comz - 1;
end;

% transform all link voxel indices back to non-padded coordinates
for i=1:length(link)
    [x,y,z] = ind2sub([w,l,h],link(i).point);
    link(i).point = sub2ind([w-2,l-2,h-2],x-1,y-1,z-1);
end;
end

function skel = Graph2Skel3D(node,link,w,l,h)

% create binary image
skel = false(w,l,h);

% for all nodes
for i=1:length(node)
    if(~isempty(node(i).links)) % if node has links
        skel(node(i).idx)=true; % node voxels
        a = [link(node(i).links(node(i).links>0)).point];
        if(~isempty(a))
            skel(a)=1; % edge voxels
        end;
    end;
end;

end

function LUT = FillEulerLUT

LUT(1)  =  1;
LUT(3)  = -1;
LUT(5)  = -1;
LUT(7)  =  1;
LUT(9)  = -3;
LUT(11) = -1;
LUT(13) = -1;
LUT(15) =  1;
LUT(17) = -1;
LUT(19) =  1;
LUT(21) =  1;
LUT(23) = -1;
LUT(25) =  3;
LUT(27) =  1;
LUT(29) =  1;
LUT(31) = -1;
LUT(33) = -3;
LUT(35) = -1;
LUT(37) =  3;
LUT(39) =  1;
LUT(41) =  1;
LUT(43) = -1;
LUT(45) =  3;
LUT(47) =  1;
LUT(49) = -1;
LUT(51) =  1;

LUT(53) =  1;
LUT(55) = -1;
LUT(57) =  3;
LUT(59) =  1;
LUT(61) =  1;
LUT(63) = -1;
LUT(65) = -3;
LUT(67) =  3;
LUT(69) = -1;
LUT(71) =  1;
LUT(73) =  1;
LUT(75) =  3;
LUT(77) = -1;
LUT(79) =  1;
LUT(81) = -1;
LUT(83) =  1;
LUT(85) =  1;
LUT(87) = -1;
LUT(89) =  3;
LUT(91) =  1;
LUT(93) =  1;
LUT(95) = -1;
LUT(97) =  1;
LUT(99) =  3;
LUT(101) =  3;
LUT(103) =  1;

LUT(105) =  5;
LUT(107) =  3;
LUT(109) =  3;
LUT(111) =  1;
LUT(113) = -1;
LUT(115) =  1;
LUT(117) =  1;
LUT(119) = -1;
LUT(121) =  3;
LUT(123) =  1;
LUT(125) =  1;
LUT(127) = -1;
LUT(129) = -7;
LUT(131) = -1;
LUT(133) = -1;
LUT(135) =  1;
LUT(137) = -3;
LUT(139) = -1;
LUT(141) = -1;
LUT(143) =  1;
LUT(145) = -1;
LUT(147) =  1;
LUT(149) =  1;
LUT(151) = -1;
LUT(153) =  3;
LUT(155) =  1;

LUT(157) =  1;
LUT(159) = -1;
LUT(161) = -3;
LUT(163) = -1;
LUT(165) =  3;
LUT(167) =  1;
LUT(169) =  1;
LUT(171) = -1;
LUT(173) =  3;
LUT(175) =  1;
LUT(177) = -1;
LUT(179) =  1;
LUT(181) =  1;
LUT(183) = -1;
LUT(185) =  3;
LUT(187) =  1;
LUT(189) =  1;
LUT(191) = -1;
LUT(193) = -3;
LUT(195) =  3;
LUT(197) = -1;
LUT(199) =  1;
LUT(201) =  1;
LUT(203) =  3;
LUT(205) = -1;
LUT(207) =  1;

LUT(209) = -1;
LUT(211) =  1;
LUT(213) =  1;
LUT(215) = -1;
LUT(217) =  3;
LUT(219) =  1;
LUT(221) =  1;
LUT(223) = -1;
LUT(225) =  1;
LUT(227) =  3;
LUT(229) =  3;
LUT(231) =  1;
LUT(233) =  5;
LUT(235) =  3;
LUT(237) =  3;
LUT(239) =  1;
LUT(241) = -1;
LUT(243) =  1;
LUT(245) =  1;
LUT(247) = -1;
LUT(249) =  3;
LUT(251) =  1;
LUT(253) =  1;
LUT(255) = -1;
end

function skel = Skeleton3D(img,spare)
% SKELETON3D Calculate the 3D skeleton of an arbitrary binary volume using parallel medial axis thinning.
%
% skel = SKELETON3D(img) returns the skeleton of the binary volume 'img'
% skel = SKELETON3D(img,mask) preserves foreground voxels in 'mask'
%
% MATLAB vectorized implementation of the algorithm by Lee, Kashyap and Chu
% "Building skeleton models via 3-D medial surface/axis thinning algorithms."
% Computer Vision, Graphics, and Image Processing, 56(6):462–478, 1994.
%
% Inspired by the ITK implementation by Hanno Homann
% http://hdl.handle.net/1926/1292
% and the Fiji/ImageJ plugin by Ignacio Arganda-Carreras
% http://fiji.sc/wiki/index.php/Skeletonize3D
%
% Philip Kollmannsberger (philipk@gmx.net)
%
% For more information, see <a
% href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d')">Skeleton3D</a> at the MATLAB File Exchange.

% pad volume with zeros to avoid edge effects
skel=padarray(img,[1 1 1]);

if(nargin==2)
    spare=padarray(spare,[1 1 1]);
end;

% fill lookup table
eulerLUT = FillEulerLUT;

width = size(skel,1);
height = size(skel,2);
depth = size(skel,3);

unchangedBorders = 0;

while( unchangedBorders < 6 )  % loop until no change for all six border types
    unchangedBorders = 0;
    for currentBorder=1:6 % loop over all 6 directions
        cands=zeros(width,height,depth);
        switch currentBorder
            case 4,
                x=2:size(skel,1); % identify border voxels as candidates
                cands(x,:,:)=skel(x,:,:) - skel(x-1,:,:);
            case 3,
                x=1:size(skel,1)-1;
                cands(x,:,:)=skel(x,:,:) - skel(x+1,:,:);
            case 1,
                y=2:size(skel,2);
                cands(:,y,:)=skel(:,y,:) - skel(:,y-1,:);
            case 2,
                y=1:size(skel,2)-1;
                cands(:,y,:)=skel(:,y,:) - skel(:,y+1,:);
            case 6,
                z=2:size(skel,3);
                cands(:,:,z)=skel(:,:,z) - skel(:,:,z-1);
            case 5,
                z=1:size(skel,3)-1;
                cands(:,:,z)=skel(:,:,z) - skel(:,:,z+1);
        end;
        
        % if excluded voxels were passed, remove them from candidates
        if(nargin==2)
            cands = cands.*~spare;
        end;
        
        % make sure all candidates are indeed foreground voxels
        cands = intersect(find(cands(:)==1),find(skel(:)==1));
        
        noChange = true;
                    
        if(~isempty(cands))
            % get subscript indices of candidates
            [x,y,z]=ind2sub([width height depth],cands);
            
            % get 26-neighbourhood of candidates in volume
            nhood = logical(pk_get_nh(skel,cands));
            
            % remove all endpoints (exactly one nb) from list
            di1 = find(sum(nhood,2)==2);
            nhood(di1,:)=[];
            cands(di1)=[];
            x(di1)=[];
            y(di1)=[];
            z(di1)=[];
            
            % remove all non-Euler-invariant points from list
            di2 = find(~p_EulerInv(nhood, eulerLUT'));
            nhood(di2,:)=[];
            cands(di2)=[];
            x(di2)=[];
            y(di2)=[];
            z(di2)=[];
            
            % remove all non-simple points from list
            di3 = find(~p_is_simple(nhood));
            nhood(di3,:)=[];
            cands(di3)=[];
            x(di3)=[];
            y(di3)=[];
            z(di3)=[];
            
            
            % if any candidates left: divide into 8 independent subvolumes
            if(~isempty(x))
                x1 = find(mod(x,2));
                x2 = find(~mod(x,2));
                y1 = find(mod(y,2));
                y2 = find(~mod(y,2));
                z1 = find(mod(z,2));
                z2 = find(~mod(z,2));
                ilst(1).l = intersect(x1,intersect(y1,z1));
                ilst(2).l = intersect(x2,intersect(y1,z1));
                ilst(3).l = intersect(x1,intersect(y2,z1));
                ilst(4).l = intersect(x2,intersect(y2,z1));
                ilst(5).l = intersect(x1,intersect(y1,z2));
                ilst(6).l = intersect(x2,intersect(y1,z2));
                ilst(7).l = intersect(x1,intersect(y2,z2));
                ilst(8).l = intersect(x2,intersect(y2,z2));
                
                idx = [];
                
                % do parallel re-checking for all points in each subvolume
                for i = 1:8                    
                    if(~isempty(ilst(i).l))
                        idx = ilst(i).l;
                        li = sub2ind([width height depth],x(idx),y(idx),z(idx));
                        skel(li)=0; % remove points
                        nh = logical(pk_get_nh(skel,li));
                        di_rc = find(~p_is_simple(nh));
                        if(~isempty(di_rc)) % if topology changed: revert
                            skel(li(di_rc))=1;
                        else
                            noChange = false; % at least one voxel removed
                        end;
                    end;
                end;
            end;
        end;
        
        if( noChange )
            unchangedBorders = unchangedBorders + 1;
        end;
        
    end;
end;

% get rid of padded zeros
skel = skel(2:end-1,2:end-1,2:end-1);


end


function [curvature] = findPointNormals(points, numNeighbours, viewPoint, dirLargest) %% check inputs
validateattributes(points, {'numeric'},{'ncols',3});

if(nargin < 2)
    numNeighbours = [];
end
if(isempty(numNeighbours))
    numNeighbours = 9;
else
    validateattributes(numNeighbours, {'numeric'},{'scalar','positive'});
    if(numNeighbours > 100)
        warning(['%i neighbouring points will be used in plane'...
            ' estimation, expect long run times, large ram usage and'...
            ' poor results near edges'],numNeighbours);
    end
end

if(nargin < 3)
    viewPoint = [];
end
if(isempty(viewPoint))
    viewPoint = [0,0,0];
else
    validateattributes(viewPoint, {'numeric'},{'size',[1,3]});
end

if(nargin < 4)
    dirLargest = [];
end
if(isempty(dirLargest))
    dirLargest = true;
else
    validateattributes(dirLargest, {'logical'},{'scalar'});
end

%% setup

%ensure inputs of correct type
points = double(points);
viewPoint = double(viewPoint);

%create kdtree
kdtreeobj = KDTreeSearcher(points,'distance','euclidean');

%get nearest neighbours
n = knnsearch(kdtreeobj,points,'k',(numNeighbours+1));

%remove self
n = n(:,2:end);

%find difference in position from neighbouring points
p = repmat(points(:,1:3),numNeighbours,1) - points(n(:),1:3);
p = reshape(p, size(points,1),numNeighbours,3);

%calculate values for covariance matrix
C = zeros(size(points,1),6);
C(:,1) = sum(p(:,:,1).*p(:,:,1),2);
C(:,2) = sum(p(:,:,1).*p(:,:,2),2);
C(:,3) = sum(p(:,:,1).*p(:,:,3),2);
C(:,4) = sum(p(:,:,2).*p(:,:,2),2);
C(:,5) = sum(p(:,:,2).*p(:,:,3),2);
C(:,6) = sum(p(:,:,3).*p(:,:,3),2);
C = C ./ numNeighbours;

%% normals and curvature calculation

curvature = zeros(size(points,1),1);
for i = 1:(size(points,1))
    
    %form covariance matrix
    Cmat = [C(i,1) C(i,2) C(i,3);...
        C(i,2) C(i,4) C(i,5);...
        C(i,3) C(i,5) C(i,6)];  
    
    %get eigen values and vectors
    [v,d] = eig(Cmat);
    d = diag(d);
    [lambda,k] = min(d);
    
    
    %store curvature
    curvature(i) = lambda / sum(d);
end


end