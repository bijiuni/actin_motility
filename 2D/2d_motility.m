function k=LineCurvature2D(Vertices,Lines,varargin)

error(nargchk(0,1,nargin));
if (nargin == 0),
     [filename, pathname] = uigetfile( ...
	       {'*.jpg;*.tif;*.gif;*.png;*.bmp', ...
		'All MATLAB Image Files (*.jpg,*.tif,*.gif,*.png,*.bmp)'; ...
		'*.jpg;*.jpeg', ...
		'JPEG Files (*.jpg,*.jpeg)'; ...
		'*.tif;*.tiff', ...
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
     if isequal(filename,0) | isequal(pathname,0)
	  return
     else
	  imagename = fullfile(pathname, filename);
     end
elseif nargin == 1,
     imagename = varargin{1};
     [path, file,ext] = fileparts(imagename);
     filename = strcat(file,ext);
end

BW = imread(imagename);
image(BW)
FigName = ['IMAGE: ' filename];

BW2 = imcomplement(BW);
imshow (BW2,[],'InitialMagnification',100);
OriginButton = questdlg('Select a point with same x value with the starting point but under it', ...
			'DIGITIZE: user input required', ...
			'OK','Cancel','OK');
switch OriginButton,
     case 'OK',
	  drawnow
	  [Xopixels,Yopixels] = ginput(1);

     case 'Cancel',
	  close(FigName)
	  return
end

Xopixels=round(Xopixels);
Yopixels=round(Yopixels);


[rows,colums] = find (BW2(1:Yopixels,Xopixels)==0,1,'last');
r = rows;
c = Xopixels;


contour = bwtraceboundary(BW,[r c],'W',8,Inf,'counterclockwise');
hold on;
plot(contour(:,2),contour(:,1),'g','LineWidth',2);

SecondButton = questdlg(...
	 'Select a point with same x value with ending point and very near it', ...
			'DIGITIZE: user input required', ...
			'OK','Cancel','OK');

    
        
switch SecondButton,
     case 'OK',
	  drawnow
	  [Xend,Yend] = ginput(1);
     Xend=round(Xend);
     Yend=round(Yend);
     
      bool=1;
          for j=1:size(contour,1) 
              if contour(j,2)==Xend
                  
                  save(bool)=j;
                  bool = bool +1;
              
              end
          end    

     case 'Cancel',
	  close(FigName)
	  return
      
end
  
[minabsolute, la2]=min (abs(contour (save,1)-Yend));

final = save(la2);

        
middle (:,1) = contour (:,2);
middle (:,2) = contour (:,1);

plot(middle(1:final,1),middle(1:final,2),'r','LineWidth',2);


Vertices = middle (1:4:final,1:2);
Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];

% Get left and right neighbor of each points
Na=zeros(size(Vertices,1),1); Nb=zeros(size(Vertices,1),1);
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);

% Check for end of line points, without a left or right neighbor
checkNa=Na==0; checkNb=Nb==0;
Naa=Na; Nbb=Nb;
Naa(checkNa)=find(checkNa); Nbb(checkNb)=find(checkNb);

% If no left neighbor use two right neighbors, and the same for right... 
Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)=Naa(Naa(checkNb));

% Correct for sampeling differences
Ta=-sqrt(sum((Vertices-Vertices(Na,:)).^2,2));
Tb=sqrt(sum((Vertices-Vertices(Nb,:)).^2,2)); 

% If no left neighbor use two right neighbors, and the same for right... 
Ta(checkNa)=-Ta(checkNa); Tb(checkNb)=-Tb(checkNb);

% Fit a polygons to the vertices 
% x=a(3)*t^2 + a(2)*t + a(1) 
% y=b(3)*t^2 + b(2)*t + b(1) 
% we know the x,y of every vertice and set t=0 for the vertices, and
% t=Ta for left vertices, and t=Tb for right vertices,  
x = [Vertices(Na,1) Vertices(:,1) Vertices(Nb,1)];
y = [Vertices(Na,2) Vertices(:,2) Vertices(Nb,2)];
M = [ones(size(Tb)) -Ta Ta.^2 ones(size(Tb)) zeros(size(Tb)) zeros(size(Tb)) ones(size(Tb)) -Tb Tb.^2];
invM=inverse3(M);
a(:,1)=invM(:,1,1).*x(:,1)+invM(:,2,1).*x(:,2)+invM(:,3,1).*x(:,3);
a(:,2)=invM(:,1,2).*x(:,1)+invM(:,2,2).*x(:,2)+invM(:,3,2).*x(:,3);
a(:,3)=invM(:,1,3).*x(:,1)+invM(:,2,3).*x(:,2)+invM(:,3,3).*x(:,3);
b(:,1)=invM(:,1,1).*y(:,1)+invM(:,2,1).*y(:,2)+invM(:,3,1).*y(:,3);
b(:,2)=invM(:,1,2).*y(:,1)+invM(:,2,2).*y(:,2)+invM(:,3,2).*y(:,3);
b(:,3)=invM(:,1,3).*y(:,1)+invM(:,2,3).*y(:,2)+invM(:,3,3).*y(:,3);

% Calculate the curvature from the fitted polygon
k = 2*(a(:,2).*b(:,3)-a(:,3).*b(:,2)) ./ ((a(:,2).^2+b(:,2).^2).^(3/2));


end



function  Minv  = inverse3(M)
% This function does inv(M) , but then for an array of 3x3 matrices
adjM(:,1,1)=  M(:,5).*M(:,9)-M(:,8).*M(:,6);
adjM(:,1,2)=  -(M(:,4).*M(:,9)-M(:,7).*M(:,6));
adjM(:,1,3)=  M(:,4).*M(:,8)-M(:,7).*M(:,5);
adjM(:,2,1)=  -(M(:,2).*M(:,9)-M(:,8).*M(:,3));
adjM(:,2,2)=  M(:,1).*M(:,9)-M(:,7).*M(:,3);
adjM(:,2,3)=  -(M(:,1).*M(:,8)-M(:,7).*M(:,2));
adjM(:,3,1)=  M(:,2).*M(:,6)-M(:,5).*M(:,3);
adjM(:,3,2)=  -(M(:,1).*M(:,6)-M(:,4).*M(:,3));
adjM(:,3,3)=  M(:,1).*M(:,5)-M(:,4).*M(:,2);
detM=M(:,1).*M(:,5).*M(:,9)-M(:,1).*M(:,8).*M(:,6)-M(:,4).*M(:,2).*M(:,9)+M(:,4).*M(:,8).*M(:,3)+M(:,7).*M(:,2).*M(:,6)-M(:,7).*M(:,5).*M(:,3);
Minv=bsxfun(@rdivide,adjM,detM);
end
