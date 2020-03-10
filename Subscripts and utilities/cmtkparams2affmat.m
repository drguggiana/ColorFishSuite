function affmat = cmtkparams2affmat(cmtk)

%cmtk is a matrix with the output parameters of cmtk
%1)translation
tx = cmtk(1,1);
ty = cmtk(1,2);
tz = cmtk(1,3);
%2) rotation
rx = cmtk(2,1);
ry = cmtk(2,2);
rz = cmtk(2,3);
%3) scale
sx = cmtk(3,1);
sy = cmtk(3,2);
sz = cmtk(3,3);
%4) shear
shx = cmtk(4,1);
shy = cmtk(4,2);
shz = cmtk(4,3);
%5) center
cx = cmtk(5,1);
cy = cmtk(5,2);
cz = cmtk(5,3);

%turn the angles from degrees to radians
alpha = deg2rad(rx);
theta = deg2rad(ry);
phi = deg2rad(rz);
  
cos0 = cos(alpha);
sin0 = sin(alpha);
cos1 = cos(theta);
sin1 = sin(theta);
cos2 = cos(phi);
sin2 = sin(phi);

sin0xsin1 = sin0 * sin1;
cos0xsin1 = cos0 * sin1;
  
% rval = zeros(4,4);
% diag(rval)<-1
rval = diag(ones(4,1));
  % nb in R matrices are indexed m[row,col]
  % whereas in C looks like T indexed them
  % m[col-1][row-1]
  % Regexps to transform these forms:
  % \[\d\]\[\d\]
  % (\2+1,\1+1)
  % 
  rval(1,1) =  cos1*cos2;
  rval(2,1) = -cos1*sin2;
  rval(3,1) = -sin1;
  rval(1,2) =  (sin0xsin1 * cos2 + cos0*sin2);
  rval(2,2) = (-sin0xsin1*sin2 + cos0*cos2);
  rval(3,2)=  sin0*cos1;
  rval(1,3) =  (cos0xsin1*cos2 - sin0*sin2);
  rval(2,3) = (-cos0xsin1*sin2 - sin0*cos2);
  rval(3,3) =  cos0*cos1;
  
  
  % generate scales and shears according to CMTK >=v.2.4.0 / svn r5050
%     scaleShear=matrix(0,4,4)
%     diag(scaleShear)=c(sx,sy,sz,1)
scaleShear = diag([sx,sy,sz,1]);
scaleShear(1,2) = shx;
scaleShear(1,3) = shy;
scaleShear(2,3) = shz;
    

%   cM = rval(1:3,1:3)*[cx;cy;cz];
%   rval(1,4) = tx - cM(1) + cx;
%   rval(2,4) = ty - cM(2) + cy;
%   rval(3,4) = tz - cM(3) + cz;

    % NB matrix multiplication must be in opposite order from C original
    rval = rval*scaleShear;
%     rval = scaleShear*rval;
    
     % transform rotation center
  cM = [  cx*rval(1,1) + cy*rval(1,2) + cz*rval(1,3),...
           cx*rval(2,1) + cy*rval(2,2) + cz*rval(2,3),...
           cx*rval(3,1) + cy*rval(3,2) + cz*rval(3,3)  ];

%   cM = rval(1:3,1:3)*[cx;cy;cz];
  
  % set translations
  rval(1,4) = tx - cM(1) + cx;
  rval(2,4) = ty - cM(2) + cy;
  rval(3,4) = tz - cM(3) + cz;
  
  affmat = rval';