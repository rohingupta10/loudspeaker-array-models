N=input('Enter Order of Soundfield\n');
r=input('Enter radius of Sphere where point source modeled\n');
k=input('Enter Wavenumber of sound\n');
m=625;
theta=zeros(1,m);
phi=zeros(1,m);
s=zeros(m,1);
lambda=2;
thetals=0.55;
phils=0.62;
c=1;
for i=1:m
    if mod(i,25)==0
        phi(i)=c*2*pi/25;
        c=0;
    else
        phi(i)=c*2*pi/25;
    end
    c=c+1;
    theta(i)=pi*ceil(i/25)/25;
end
sphericalharmonics=(N+1)*(N+1);
H=zeros(sphericalharmonics,m);
alpha=sh2(N,thetals,phils);

for i = 1:m
    Y=sh2(N,theta(i),phi(i));
    %disp(max(Y));
    H(1:(N+1)^2 ,i)=conj(Y);
    j=1;
    order=0;
    while j<sphericalharmonics
        terms=0;
        while terms<2*order+1
            H(terms+j,i)=4*pi*sqrt(-1)*k*H(terms+j,i)*besselhs(order,k*r);
            if i==1
                alpha(terms+j)=4*pi*sqrt(-1)*k*alpha(terms+j)*besselhs(order,k*r);
                %disp(besselhs(order,k*r));
            end
            terms=terms+1;
        end
        j=j+2*order+1;
        order=order+1;
    end
end
alpha=real(alpha);
H=real(H);
epsilon=0.001;
cvx_begin
    variable s(m)
    minimize( norm( H * s - alpha, 2 ) + lambda * norm(s,1));
    subject to
        norm( H * s - alpha, 2 ) <= epsilon
cvx_end