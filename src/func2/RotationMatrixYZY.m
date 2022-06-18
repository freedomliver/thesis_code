function [Rotation] = RotationMatrixYZY(alpha,beta,gamma)
%{
Rotation Matrix
Ry*Rz*Ry
y1-z′-y2″ (intrinsic rotations) or y2-z-y1 (extrinsic rotations) from wiki

alpha, beta, gamma: Euler angles
Rotation: rotation matrix

2021/08/25 qifengfeng
%}

s1=sin(alpha);
s2=sin(beta);
s3=sin(gamma);
c1=cos(alpha);
c2=cos(beta);
c3=cos(gamma);

Rotation=[c1*c2*c3-s1*s3,-1*c1*s2,c3*s1+c1*c2*s3;
          c3*s2,c2,s2*s3;
          -1*c1*s3-c2*c3*s1,s1*s2,c1*c3-c2*s1*s3];

end

