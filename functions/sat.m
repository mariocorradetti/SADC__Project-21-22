function sat(Att,r,scal)
% % Function that creates the shape of s/c 
% % rotate body to inertial
A=0.5*Att*[0.3405 0.1 0.1]';
B=0.5*Att*[0.3405 -0.1 0.1]';
C=0.5*Att*[0.3405 -0.1 -0.1]';
D=0.5*Att*[0.3405 0.1 -0.1]';
E=0.5*Att*[-0.3405 0.1 0.1]';
F=0.5*Att*[-0.3405 -0.1 0.1]';
G=0.5*Att*[-0.3405 -0.1 -0.1]';
H=0.5*Att*[-0.3405 0.1 -0.1]';
% 
X=scal.*[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=scal.*[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=scal.*[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
Co=[99/255,99/255,99/255]; 
Xs = [B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Ys = [B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Zs = [B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
% 
xcg_b=r(1);
ycg_b=r(2);
zcg_b=r(3);

Xbody= xcg_b+X;
Ybody= ycg_b+Y;
Zbody= zcg_b+Z;
Xsp= xcg_b+Xs;
Ysp= ycg_b+Ys;
Zsp= zcg_b+Zs;

fill3(Xbody,Ybody,Zbody,Co,'FaceAlpha',1)
hold on
Co=[20/255,80/255,160/255]; 
fill3(Xsp,Ysp,Zsp,Co,'FaceAlpha',1)
axis equal

%% Solar panel

cg1 = Att*[0.3405/2-0.01/2 , 0.1/2+0.3405/2,0]';

xy1= xcg_b+cg1(1)*scal;
yy1= ycg_b+cg1(2)*scal;
zy1= zcg_b+cg1(3)*scal;

A=0.5*Att*[0.01 0.3405 0.1]';
B=0.5*Att*[-0.01 0.3405 0.1]';
C=0.5*Att*[-0.01 0.3405 -0.1]';
D=0.5*Att*[0.01 0.3405 -0.1]';
E=0.5*Att*[0.01 -0.3405 0.1]';
F=0.5*Att*[-0.01 -0.3405 0.1]';
G=0.5*Att*[-0.01 -0.3405 -0.1]';
H=0.5*Att*[0.01 -0.3405 -0.1]';

X=scal.*[A(1) B(1) C(1) D(1); D(1) C(1) G(1) H(1); G(1) H(1) E(1) F(1); E(1) F(1) B(1) A(1); B(1) C(1) G(1) F(1); A(1) D(1) H(1) E(1)]';
Y=scal.*[A(2) B(2) C(2) D(2); D(2) C(2) G(2) H(2); G(2) H(2) E(2) F(2); E(2) F(2) B(2) A(2); B(2) C(2) G(2) F(2); A(2) D(2) H(2) E(2)]';
Z=scal.*[A(3) B(3) C(3) D(3); D(3) C(3) G(3) H(3); G(3) H(3) E(3) F(3); E(3) F(3) B(3) A(3); B(3) C(3) G(3) F(3); A(3) D(3) H(3) E(3)]';
Co=[20/255,80/255,160/255]; 
Xy1=xy1+X;
Yy1=yy1+Y;
Zy1=zy1+Z;
%-----------------------------------------------------------------% 
cg2 = Att*[+0.3405/2-0.01/2 , -0.1/2-0.3405/2,0]';

xy2= xcg_b+cg2(1)*scal;
yy2= ycg_b+cg2(2)*scal;
zy2= zcg_b+cg2(3)*scal;

Xy2=xy2+X;
Yy2=yy2+Y;
Zy2=zy2+Z;
%-----------------------------------------------------------------%

fill3(Xy2,Yy2,Zy2,Co,'FaceAlpha',1)
fill3(Xy1,Yy1,Zy1,Co,'FaceAlpha',1)

hold on

line(r(1)+[0 Att(1,1)]./3,r(2)+[0 Att(2,1)]./3,r(3)+[0 Att(3,1)]./3,'Color','k','Linewidth',2)
line(r(1)+[0 Att(1,2)]./3,r(2)+[0 Att(2,2)]./3,r(3)+[0 Att(3,2)]./3,'Color','k','Linewidth',2)
line(r(1)+[0 Att(1,3)]./3,r(2)+[0 Att(2,3)]./3,r(3)+[0 Att(3,3)]./3,'Color','k','Linewidth',2)
% text(r(1)+Att(1,1)/3+0.05,0,0,'X')
% text(0,r(1)+Att(1,1)/3+0.05,0,'Y')
% text(0,0,r(1)+Att(1,1)/3+0.05,'Z')

xlabel('x')
ylabel('y')
zlabel('z')
grid on

% [x,y,z] = cylinder([0 0.1],10);
% x = (x+r(1));
% y = (y+r(2));
% z =-z/2 - 0.3405.*scal/2;
% z = (z+r(3));
% cyl_base = Att*[x(1,:)',y(1,:)',z(1,:)']';
% cyl_top = Att*[x(2,:)',y(2,:)',z(2,:)']';
% x =[cyl_base(1,:),cyl_top(1,:)];
% y =[cyl_base(2,:),cyl_top(2,:)];
% z =[cyl_base(3,:),cyl_top(3,:)];
% fill3(z,x,y,Co)
end