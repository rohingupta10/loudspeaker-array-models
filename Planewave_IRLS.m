N=input('Enter Order of Soundfield\n');
r=input('Enter radius of Sphere where point source modeled\n');
m=625;
theta_mike=[69, 90, 111, 90, 32, 55, 90, 125, 148, 125, 90, 55, 21, 58, 121, 159, 69, 90, 111, 90, 32, 55, 90, 125, 148, 125, 90, 55, 21, 58, 122, 159];
phi_mike=[0, 32, 0, 328, 0, 45, 69, 45, 0, 315, 291, 315, 91, 90, 90, 89, 180, 212, 180, 148, 180, 225, 249, 225, 180, 135, 111, 135, 269, 270, 270, 271];
theta=zeros(1,m);
phi=zeros(1,m);
s=zeros(m,1);
alpha=zeros((N+1)^2,1);
thetals=0.55;
phils=0.62;
thetamike=zeros(5,32);
phimike=zeros(5,32);
pressureth=zeros(160,1);
pressurecal=zeros(160,1);


for i=1:5
    x0=-0.3*cos(2*pi*i/5);
    y0=-0.3*sin(2*pi*i/5);
    z0=0;
    for j=1:32
         x1=0.04*sin(theta_mike(j))*cos(phi_mike(j));
         y1=0.04*sin(theta_mike(j))*sin(phi_mike(j));
         z1=0.04*sin(phi_mike(j));
         x2=x1-x0;
         y2=y1-y0;
         z2=z1-z0;
         elevation=acos(z2/sqrt(x2^2 +y2^2 +z2^2));
         azimuth=atan(y2/z2);
         thetamike(i,j)=elevation;
         phimike(i,j)=azimuth;
    end
end
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

frequency=zeros(181,0);
error=zeros(181,0);
c=1;

for freq = 200:5:1000
    
    for i=1:m
        s(i)=rand()+1;
    end
    
    k=2*pi*freq/340;
    sphericalharmonics=(N+1)*(N+1);
    H=zeros(sphericalharmonics,m);
    alpha=sh2(N,thetals,phils);
    for i = 1:m
        Y=sh2(N,theta(i),phi(i));
        j=1;
        order=0;
        H(1:(N+1)^2 ,i)=conj(Y);
        j=1;
        order=0;
        while j<sphericalharmonics
            terms=0;
            while terms<2*order+1
                H(terms+j,i)=4*pi*(sqrt(-1)^order)*H(terms+j,i);
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
    %H=real(H);
    %alpha=real(alpha);
    iter=100;
    p=0.8;
    w=zeros(1,m);  
    Q=zeros(m,m);
    for i = 1:iter
        for j =1:m
            w(j)=(abs(s(j)))^(p-2);
            Q(j,j)=1/w(j);
        end
        s=Q*(H.')*pinv(H*Q*(H.'))*alpha;
    end
    
    %disp(s);
    alpha=H*s;
    for i=1:5
        for j=1:32
            pressureth((i-1)*32 +j)=0;
            for z =1:m
                y=[sin(theta(i))*cos(phi(i)),sin(theta(i))*sin(phi(i)),cos(theta(i))];
                x=[0.3*sin(thetamike(i,j))*cos(phimike(i,j)),0.3*sin(thetamike(i,j))*sin(phimike(i,j)),0.3*cos(thetamike(i,j))];
                dot=y*(x');
                pressureth((i-1)*32+j)=pressureth((i-1)*32+j)+ s(z)*exp(sqrt(-1)*k*dot);
            end
            
            Y=sh2(N,thetamike(i,j),phimike(i,j));
            
            for a=0:N
                for z=-a:a
                    pressurecal((i-1)*32+j)=pressurecal((i-1)*32+j)+ Y(a^2+a+1+z)*besseljs(order,k*0.3)*alpha(a^2+a+1+z);
            
                end
            end
        end
    end
    
    pressureth=abs(pressureth);
    pressurecal=abs(pressurecal);
    num=sum((pressureth-pressurecal).*(pressureth-pressurecal));
    den=sum(pressurecal.*pressurecal);
    
    if freq==220
        disp(pressureth);
        disp(pressurecal);
    end
    error(c)=10*log10(num/den);
    frequency(c)=freq;
    c=c+1;
end

plot(frequency,error);