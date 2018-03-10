close all
clc
format rat
n=6;            % Order of the graph

B= [[0, 1, 2, 2, 3, 4], [1, 0, 1, 1, 2, 3], [2, 1, 0, 2, 1, 2], [2, 1, 2, 0, 1, 2], [3, 2, 1, 1, 0, 1], [4, 3, 2, 2, 1, 0]];



D= reshape(B,[n,n]);
if D==D'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end
l=size(D,1);
e=[];
s=[];
t=[];
k=[];
for i=1:l
    for j=i+1:l
        if D(i,j)==1
            e=[e;[i,j]];
        end
    end
end

ress = zeros(l,l);
rest = zeros(l,l);
resk = zeros(l,l);
for i=1:size(e,1)
    ss = [];
    tt = [];
    kk = [];
    for ii=1:l
        if D(e(i,1),ii)<D(e(i,2),ii)
           ss=[ss ii];
        end
        if D(e(i,1),ii)>D(e(i,2),ii)
           tt=[tt ii];
         end
        if D(e(i,1),ii)==D(e(i,2),ii)
           kk=[kk ii];
        end
    end
    ress(i, 1:length(ss)) = ss;
    rest(i, 1:length(tt)) = tt;
    resk(i, 1:length(kk)) = kk;
%     s=[s;ss];
%     t=[t;tt];
%     k=[k;kk]; 
end
res1=[e ress];
res2=[e rest];
res3=[e resk];
[sm, ~] = size(ss);
U=[];
% resU(i, 1:length(ss)) = ss;
for i=1:sm
    U = [U; length(ss(i,:))];
end
Us = (ress~=0);
U=sum(Us');
U=U';
eU=[e U];

[tm, ~] = size(t);
V=[];
for i=1:tm
    V = [V; length(t(i,:))];
end
Vt = (rest~=0);
V=sum(Vt');
V= V';
eV=[e V];
[km, ~] = size(k);
W=[];
for i=1:km
W = [W; length(k(i,:))];
end
Wk = (resk~=0);
W=sum(Wk');
W=W';
eW=[e W];

% saki=[U V W];
% saki=sortrows(saki,[-1,-2]);
siki=[e U V W];
saki=sortrows(siki,[-3,-4]);

www=sum(D);
wwx=www';
wwy=[];
for ij=1:l
    wws=[ij wwx(ij,:)];                % For Wiener-like partition
    wwy=[wwy;wws];
end


W1 = U+(1/2).*W;
W2 = V+(1/2).*W;
w = W1.*W2;
RS=sum(w);                     % Revised Szeged index
%RS1=rats(sum(w))
x = U.*V;
Szeged = sum(x);               % Szeged index
y = U+V;
PI = sum(y);                   % PI index

edges=size(saki,1);
j=ones(edges,1);
a=y-2*j;
aa=sqrt(a./x);
abc2=sum(aa);                  % ABC2 index
b=2*sqrt(x);
bb=b./y;
ga2=sum(bb);                    % GA2 index

wiener=(0.5)*(sum(sum(D)));    % Wiener index
M = max(D);                    % eccentriciy vector
d=[];
for i=1:l
    r=D(i,:);
    dd=nnz(r==1);
    d=[d;dd];
end
ecc=sum(M);                   % Total eccentricity index
ecv=[M' d];
ex = M'.*d;
tec=sum(ex);                  % Eccentric connectivity index
G=[];
for i=1:edges
    GG=[e(i,:) M(e(i,1)) M(e(i,2))];
    G=[G;GG];
end
S=G(:,3);
T=G(:,4);
c1= S.*T;
c2=S+T;
c3=c2-2*j;
c5=sqrt(c3./c1);
abc5=sum(c5);                  % ABC5 index

cc1=2*sqrt(c1);
cc2=cc1./c2;
ga4=sum(cc2);                 % GA5 index

diam=max(max(D));
A=[];
for r=2:diam
    D(D==r)=0;
    A=D;
end
D= reshape(B,[n,n]);
k=sum(A);
MTI=sum((k)*(D+A));       % Schultz Molecular TI

ww=[];
xx=[];
yy=[];
zz=[];
for i=1:n
    for j=1:n
        if i~=j
            DegD=(((k(i))+k(j))*D(i,j));
            Gut=(((k(i))*k(j))*D(i,j));
            HA=(((k(i))+k(j))./D(i,j));
            HM=(((k(i))*k(j))./D(i,j));
            ww=[ww;DegD];
            xx=[xx;Gut];
            yy=[yy;HA];
            zz=[zz;HM];
        end
    end
end
DegDis=(0.5)*sum(ww);                                          % Degree Distance TI
Gutman=(0.5)*sum(xx);                                          % Gutman TI
HararyA=(0.5)*sum(yy);                                         % Additive Harary TI
HararyM=(0.5)*sum(zz);                                         % Multiplicative Harary TI

% fprintf('n =%d\n',n');
fprintf('The Wiener Index is =%d\n',wiener');

fprintf('The Szeged Index is =%d\n',Szeged');

fprintf('The PI Index is =%d\n',PI');

fprintf('The Revised Szeged Index is %4.4f\n',RS');

fprintf('The eccentric connectivity index is =%d\n',tec');

fprintf('The total eccentricity index is =%d\n',ecc');

fprintf('The Second ABC index is %4.4f\n',abc2');

fprintf('The Second GA index is %4.4f\n',ga2');

fprintf('The fifth ABC index is %4.4f\n',abc5');

fprintf('The fourth GA index is %4.4f\n',ga4');

fprintf('The Schultz Molecular Index is =%d\n',MTI');

fprintf('The Degree Distance Index is =%d\n',DegDis');

fprintf('The Gutman Index is =%d\n',Gutman');

fprintf('The Additive Harary Index is %4.4f\n',HararyA');

fprintf('The Multiplicative Harary Index is %4.4f\n',HararyM');
