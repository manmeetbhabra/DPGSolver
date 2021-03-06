format shorte;
clc, clf, hold on;

d = 3;

view([10 40 30])
daspect([1 1 1])

%% std

nt = 15;
np = 8;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZsq = sqrt(2)*[sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 sqrt(2); 0 0 0];
% XYZsq = [-1.0 -1.0/sqrt(3.0) -1.0/sqrt(6.0);
% 1.0 -1.0/sqrt(3.0) -1.0/sqrt(6.0);
% 0.0 0.0  -1.0/sqrt(6.0)+sqrt(2/3);
% 0.0 2.0/sqrt(3.0) -1.0/sqrt(6.0)]; % same result



LHS = [2*XYZsq ones(d+1,1)];
RHS = sum(XYZsq.^2,2);
abcd = LHS\RHS;
r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));
% [abcd(1:d)' r]

subplot(1,2,1); 
hold on;
view([10 80 30])
daspect([1 1 1])

XYZp = XYZsq;
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');
XYZp = XYZp([3 1 4 2],:);
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');

XYZ = zeros((np+1)*(nt+1),d);
XYZ(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZ(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZ(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ro');


XYZ = 1/3*[sum(XYZsq([1 2 3],:));
           sum(XYZsq([1 2 4],:));
           sum(XYZsq([1 3 4],:));
           sum(XYZsq([2 3 4],:))];

LHS = [2*XYZ ones(d+1,1)];
RHS = sum(XYZ.^2,2);
abcd = LHS\RHS;
r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));
% [abcd(1:d)' r]

% Computed from below (largest contained sphere):
abcd(1:d) = [ 0  0 -1.093897997411785e-001];
r = 2.988584907226846e-001; 
rIn = r;


nt = 150;
np = 80;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZ = zeros((np+1)*(nt+1),d);
XYZ(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZ(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZ(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ko');


% Smallest possible enclosing sphere
abcd(1:d) = [0 0 -1/sqrt(6)];
r = 2*sqrt(3)/3;
rOut = r;

nt = 15;
np = 8;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZ = zeros((np+1)*(nt+1),d);
XYZ(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZ(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZ(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'go');

% disp('12T:')
% [rIn rOut rOut/rIn]


%% eq

nt = 15;
np = 8;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZeq = [-1.0 -1.0/sqrt(3.0) -1.0/sqrt(6.0);
1.0 -1.0/sqrt(3.0) -1.0/sqrt(6.0);
0.0 2.0/sqrt(3.0) -1.0/sqrt(6.0);
0.0 0.0  3.0/sqrt(6.0)];

LHS = [2*XYZeq ones(d+1,1)];
RHS = sum(XYZeq.^2,2);
abcd = LHS\RHS;
r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));
rOut = r;
% [abcd(1:d)' r]


subplot(1,2,2);
hold on;
view([10 80 30])
daspect([1 1 1])

XYZp = XYZeq;
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');
XYZp = XYZp([3 1 4 2],:);
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');

XYZ = zeros((np+1)*(nt+1),d);
XYZ(:,1) = r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZ(:,2) = r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZ(:,3) = r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ro');





XYZ = 1/3*[sum(XYZeq([1 2 3],:));
           sum(XYZeq([1 2 4],:));
           sum(XYZeq([1 3 4],:));
           sum(XYZeq([2 3 4],:))];

LHS = [2*XYZ ones(d+1,1)];
RHS = sum(XYZ.^2,2);
abcd = LHS\RHS;
r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));
rIn = r;
% [abcd(1:d)' r]


nt = 150;
np = 80;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZ = zeros((np+1)*(nt+1),d);
XYZ(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZ(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZ(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ko');

% disp('Eq T:')
% [rIn rOut rOut/rIn]


% error('exiting');

%% smallest contained sphere

%% 2D
% d = 2;
% Nf = d+1;
% 
% t_s = [0 0 0];
% t_t = [1 5 9]*pi/6;
% 
% nr = zeros(Nf,d);
% 
% for f = 1:Nf
%     nr(f,1) =  cos(t_s(f))*cos(t_t(f));
% 	nr(f,2) =  cos(t_s(f))*sin(t_t(f));
% end
% 
% XYZeq = [-1.0 -1.0/sqrt(3.0);
%           1.0 -1.0/sqrt(3.0);
%           0.0 2.0/sqrt(3.0)];
%       
% X = XYZeq(:,1); Y = XYZeq(:,2);
%       
% LHS = [X(3)-X(2) 0 0 -nr(1,1) -1 0;
%        Y(3)-Y(2) 0 0 -nr(1,2) 0 -1;
%        0 X(3)-X(1) 0 -nr(2,1) -1 0;
%        0 Y(3)-Y(1) 0 -nr(2,2) 0 -1;
%        0 0 X(2)-X(1) -nr(3,1) -1 0;
%        0 0 Y(2)-Y(1) -nr(3,2) 0 -1];
% 
% RHS = -[X(2) Y(2) X(1) Y(1) X(1) Y(1)]';
% 
% % LHS\RHS
% 
% LHS = XYZeq([1 3],:);
% RHS = ones(2,1);
% n = (LHS\RHS); n = n/norm(n,2);

%% 3D

XYZ = XYZsq
VOL = 8/6

XYZ = 1/2*[2*XYZsq(1,:);
           XYZsq(1,:)+XYZsq(2,:);
           XYZsq(1,:)+XYZsq(3,:);
           XYZsq(1,:)+XYZsq(4,:)]
VOL = 1/6

XYZ = 1/2*[XYZsq(1,:)+XYZsq(3,:);
           XYZsq(3,:)+XYZsq(4,:);
           XYZsq(2,:)+XYZsq(4,:);
           XYZsq(2,:)+XYZsq(3,:)]
VOL = 1/6
       error('exit');

clf;
hold on;
view([10 -80 30])
daspect([1 1 1])

d = 3;
Nf = d+1;

XYZ = [  8.7500e-01  8.7500e-01  1.0000e+00 
 8.7500e-01  9.3750e-01  9.3750e-01 
 7.5000e-01  9.3750e-01  9.3750e-01 
 8.1250e-01  9.3750e-01  1.0000e+00];

XYZp = XYZ;
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');
XYZp = XYZp([3 1 4 2],:);
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-bo');

% XYZp = [0 0 0; 1 0 0; 1 1 0; 0 1 0;
%         0 0 1; 1 0 1; 1 1 1; 0 1 1];
% plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-rs');
xlabel('x');
ylabel('y');
zlabel('z');


FNodeInds = [[2 3 4]; [1 3 4]; [1 2 4]; [1 2 3]];

nr = zeros(Nf,d);
d_p = zeros(Nf,1);

XYZcT = sum(XYZ,1)/4;

for f = 1:Nf
   XYZF =  XYZ(FNodeInds(f,:),:);
   
   Vec1 = XYZF(1,:)-XYZF(3,:);
   Vec2 = XYZF(2,:)-XYZF(3,:);
   
   n       = cross(Vec1,Vec2);
   nr(f,:) = n/norm(n,2);
   d_p(f)  = dot(nr(f,:),XYZF(3,:));
   
   % Ensure that normals point in the outward direction
%    XYZc = sum(XYZF,1)/3;
%    r = norm(XYZcT-XYZc,2);
%    
%    XYZc = XYZc + nr(f,:)*1e2*eps;

   if (nr(f,:)*XYZ(f,:)'-d_p(f) > 0.0)
%    if (norm(XYZcT-XYZc,2)-r < 0.0)
%        f
       nr(f,:) = -nr(f,:);
       d_p(f)  = -d_p(f);
   end
   
  
   XF = 1/3*sum(XYZ(FNodeInds(f,:),1));
   YF = 1/3*sum(XYZ(FNodeInds(f,:),2));
   ZF = 1/3*sum(XYZ(FNodeInds(f,:),3));
   quiver3(XF,YF,ZF,nr(f,1),nr(f,2),nr(f,3),0.02)
end
% nr
% error('exiting');


% Note: For a plane a x + b y + c z = d:
%       Distance point to plane = r = |a x0 + b y0 + c z0 - d|/norm([a b c],2)
%       [a b c] == n (normalized, pointing towards the center):
%       r = a x0 + b y0 + c z0 - d
%       => Solve [nr -ones(4,1)]*[a b c r]' = d_p
%       If coordinates not needed, use Gaussian elimination to simplify?
%       See Wikipedia.
% Much simpler!
LHS = [nr ones(Nf,1)];
RHS = d_p;
tmp = LHS\RHS;
abcd(1:d) = tmp(1:3)';
rIn = tmp(4);



nt = 15*20;
np = 8*20;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

r = rIn;

XYZp = zeros((np+1)*(nt+1),d);
XYZp(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZp(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZp(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'ko');




Inde = [ 1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

Ne = 6;

lenE = zeros(Ne,1);

for e = 1:Ne
    lenE(e) = norm(XYZ(Inde(e,1),:)-XYZ(Inde(e,2),:),2);
end
% lenE
[~,IndlenEmax] = max(lenE);

XYZc = 1/2*(XYZ(2,:)+XYZ(4,:));

for i = 1:Nf
    r = norm(XYZ(i,:)-XYZc);
    [i r];
end

% LHS = [2*XYZ([2 3 4],1:2) ones(3,1)]
% RHS = sum(XYZ([2 3 4],1:2).^2,2)
% 




% case 1
abcd = 1/2*(XYZ(Inde(IndlenEmax,1),:)+XYZ(Inde(IndlenEmax,2),:));
r = 0.5*norm(XYZ(Inde(IndlenEmax,1),:)-XYZ(Inde(IndlenEmax,2),:),2);

nt = 15*2;
np = 8*2;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZp = zeros((np+1)*(nt+1),d);
XYZp(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZp(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZp(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);

% error('exiting');

done = 1;

for i = 1:Nf
    check = norm(XYZ(i,:)-abcd(1:d),2)-r;
    [i norm(XYZ(i,:)-abcd(1:d),2) r norm(XYZ(i,:)-abcd(1:d),2)-r];
    if (check > 0.0)
        done = 0;
        break;
    end
end

if (done)
    [0 r/rIn r rIn]
    plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'ko');
    error('found it.');
end


% case 2
nt = 15*2;
np = 8*2;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;
XYZp = zeros((np+1)*(nt+1),d);

abcF = zeros(Nf,d);
rF   = zeros(Nf,1);
for f = 1:Nf
    XYZF = XYZ(FNodeInds(f,:),:);


    abc_p = nr(f,:);

    LHS = [2*XYZF ones(d,1); abc_p 0];
    RHS = [sum(XYZF.^2,2); d_p(f)];

    abcd = LHS\RHS;
    r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));
    
    abcF(f,:) = abcd(1:d)';
    rF(f) = r;
end

[~,IndrF] = sort(rF);

for f = IndrF'
    
    abcd(1:d) = abcF(f,:)';
    r = rF(f);

%     [f-1 r abcd(1:d)']
%     [f-1 r]
    
    XYZp(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
    XYZp(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
    XYZp(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
%     plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'ro');
    
    done = 1;
    for i = 1:Nf
        check = norm(XYZ(i,:)-abcd(1:d)',2)-r;
%         [i check > 1e2*eps]
        if (check > 1e2*eps)
            done = 0;
            break;
        end
    end

    if (done)
        [1 r/rIn r rIn]
        plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'ro');
        error('found it.');
    end
end




% error('exiting');


LHS = [2*XYZ ones(d+1,1)];
RHS = sum(XYZ.^2,2);
abcd = LHS\RHS;
r = sqrt(abcd(d+1)+sum(abcd(1:d).^2));

nt = 15*2;
np = 8*2;
t = 0:2*pi/nt:2*pi;
p = 0:pi/np:pi;

XYZp = zeros((np+1)*(nt+1),d);
XYZp(:,1) = abcd(1)+r*reshape(cos(t)'*sin(p),(np+1)*(nt+1),1);
XYZp(:,2) = abcd(2)+r*reshape(sin(t)'*sin(p),(np+1)*(nt+1),1) ;
XYZp(:,3) = abcd(3)+r*reshape(ones(size(t))'*cos(p),(np+1)*(nt+1),1);
plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'go');


[3 r/rIn r rIn sqrt(3/8)*max(lenE)]


%%
1;
    

%% PYR base in split TET
error('exiting');
clc;

XYZ = [9.571669482429456e-01     4.217612826262750e-01     6.557406991565868e-01
     4.853756487228412e-01     9.157355251890671e-01     3.571167857418955e-02
     8.002804688888001e-01     7.922073295595544e-01     8.491293058687771e-01
     1.418863386272153e-01     9.594924263929030e-01     9.339932477575505e-01]
% XYZ = [1/4*[2*[XYZsq(1,:)+XYZsq(2,:); XYZsq(1,:)+XYZsq(3,:); XYZsq(1,:)+XYZsq(4,:)]; sum(XYZsq,1)]]

centroid = 1/4*(sum(XYZ,1));

XYZmid_e = 1/2*[XYZ(1,:)+XYZ(2,:);
                XYZ(1,:)+XYZ(3,:);
                XYZ(1,:)+XYZ(4,:);
                XYZ(2,:)+XYZ(3,:);
                XYZ(2,:)+XYZ(4,:);
                XYZ(3,:)+XYZ(4,:)];
            
for e = 1:6
    [e norm(centroid-XYZmid_e(e,:),2)]
end

clf;
hold on;
view([10 -80 30])
daspect([1 1 1])

plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-bo');
XYZ = XYZ([3 1 4 2],:);
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-bo');

plot3(centroid(1),centroid(2),centroid(3),'rd');
for e = 1:6
    XYZp = [XYZmid_e(e,:); centroid];
    plot3(XYZp(:,1),XYZp(:,2),XYZp(:,3),'-rs');
end

