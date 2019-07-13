M=input('Enter Order of Soundfield\n');
azimuth_source0=pi/20;
elevation_source=0.62;
n=9;
elevation_speaker=[pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4];
azimuth_speaker=pi/180*[90;70;45;-15;-30;-45;-60;-75;-90];
b=zeros((M+1)^2,1);
D=zeros((M+1)^2,n);
gain=zeros(n,1);

%S = sum(D.*D).^0.5;
%D = D./repmat(S,(M+1)^2,1);
error=zeros(21,1);
theta1=zeros(21,1);
c=1;

for z =-10:10
    azimuth_source=z*azimuth_source0;
    b=1*sh2(M,elevation_source,azimuth_source);
    b=real(b);
    gain=zeros(n,1);
    epsilon=0.001;
    R=b;
    iter=1;
    ind=zeros(n,1);
    
    for i=1:n
        D(: ,i)=sh2(M,elevation_speaker(i),azimuth_speaker(i));
        % gain(i)=rand()+1;
    end
    
    for i =1:n
        if b'*D(:,i)<=0
            continue;
        end
        D(:,i)=D(:,i)*(b'*D(:,i));
    end
    
    D=real(D);
    
    for i =1:n
       D(:,i)=D(:,i)/norm(D(:,i),2);
    end
    
    while norm(R) > epsilon && iter<=4
        innerprod=0;
        iter=iter+1;
        index=-1;
        mat=D';
        for j =1:n
            if mat(j,:)*R > innerprod &&ind(j)==0
                innerprod=mat(j,:)*R;
                index=j;
            end
        end
        
        disp(innerprod);
        if index ==-1
            continue;
        end
        gain(index)=innerprod;
        ind(index)=1;
        R=R-innerprod*D(:,index);
        %disp(gain);
    end
    %disp(gain)
    sum=0;
    for i=1:n
        sum=sum+(gain(i)*0.2821)^2;
    end
    %for i =1:n
    %    gain(i)=gain(i)*0.2821/sqrt(sum);
    %end
    for i =1:n
        error(c)= error(c)+ gain(i)*gain(i);
    end
    error(c)=20*log10(error(c));
    theta1(c)=azimuth_source*180/pi;
    c=c+1;
    %disp(gain);
end

plot(theta1,error);
ylim([-5 5]);