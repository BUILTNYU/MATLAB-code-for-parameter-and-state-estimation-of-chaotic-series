%%COPY AND PASTE THE FOLLOWING CONTENTS INTO AN M-FILE IN MATLAB

function [coeff,x,L,k]=chaoticLS(s,MSA,quad,tol,kmax,repeat)
%estimation of parameters lambda and state var x, given noisy s
%for logistic map
%quad = 0 if logistic map with 1 param, quad = 2 if quadratic map, quad = 1
%if logistic map with 2 param
%MSA = 0 if alpha = beta = 0.5 (Nakamura et al’s method), else we use MSA
%repeat is for if this is a first run or a repeat run; if repeat=[] this is
%a first run, else if it's size Coeff then it's a repeat run
%s is Sx1, where S is the number of samples 
numS=size(s,1);
k=1;
x(:,1)=s;
 
if size(repeat,1)==0
%construct data set based on logistic map function
    if quad==2
        xvar=zeros(numS-1,3);
        yvar=zeros(numS-1,1);
        for i=2:numS
            xvar(i,1)=1;
            xvar(i,2)=x(i-1,1);
            xvar(i,3)=x(i-1,1)^2;
            yvar(i,1)=x(i,1);
        end
    elseif quad==1
        xvar=zeros(numS-1,2);
        yvar=zeros(numS-1,1);
        for i=2:numS
            xvar(i,1)=x(i-1,1);
            xvar(i,2)=x(i-1,1)^2;
            yvar(i,1)=x(i,1);
        end
    else
        xvar=zeros(numS-1,1);
        yvar=zeros(numS-1,1);
        for i=2:numS
            xvar(i,1)=x(i-1,1)*(1-x(i-1,1));
            yvar(i,1)=x(i,1);
        end
        
    end
 
    coeff(:,1)=(xvar'*xvar)^-1*xvar'*yvar;
    
else
    coeff(:,1)=repeat;
end
%solve LS for Ikeda using lsqcurvefit
%coeff=[mu;a;b]
%xvar=s(1:999,:);
%yvar=s(2:1000,:);
%coeff=lsqcurvefit(@Ikedafun,[0;0;0],xvar,yvar);
L(1)=0;
for i=1:numS-1
    if quad==2
        L(1)=L(1)+.5*(x(i+1,1)-(coeff(1,1)+coeff(2,1)*x(i,1)+coeff(3,1)*x(i,1)^2))^2;
    elseif quad==1
        L(1)=L(1)+.5*(x(i+1,1)-(coeff(1,1)*x(i,1)+coeff(2,1)*x(i,1)^2))^2;
    else
        L(1)=L(1)+.5*(x(i+1,1)-(coeff(1,1)*x(i,1)*(1-x(i,1))))^2;
    end
end
 
stop1=0;
while and(stop1==0,k<=kmax)
 
    %gradient descent state estimation
    %indeterminism of the sample data
    L1=L(k);
%    for i=1:numS-1
%        if quad==1
            %L1=L1+.5*(s(i+1,1)-(coeff(1,k)+coeff(2,k)*s(i,1)+coeff(3,k)*s(i,1)^2))^2;
%            L1=L1+.5*(x(i+1,k)-(coeff(1,k)+coeff(2,k)*x(i,k)+coeff(3,k)*x(i,k)^2))^2;
%        else
            %L1=L1+.5*(s(i+1,1)-(coeff(1,k)*s(i,1)+coeff(2,k)*s(i,1)^2))^2;
%            L1=L1+.5*(x(i+1,k)-(coeff(1,k)*x(i,k)+coeff(2,k)*x(i,k)^2))^2;
%        end
%    end
    tau=1;
    delta=0.1;
    tempx1=x(:,k);
    stop=0;
    tempx2=tempx1;
    while and(stop==0,tau<=500)
        for i=1:numS
            if quad==2
                if i==1
                    tempx2(i,1)=tempx1(i,1)-delta*(-(coeff(2,k)+2*coeff(3,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)-coeff(2,k)*tempx1(i,1)-coeff(3,k)*tempx1(i,1)^2));
                elseif i<numS
                    tempx2(i,1)=tempx1(i,1)-delta*((tempx1(i,1)-coeff(1,k)-coeff(2,k)*tempx1(i-1,1)-coeff(3,k)*tempx1(i-1,1)^2)-(coeff(2,k)+2*coeff(3,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)-coeff(2,k)*tempx1(i,1)-coeff(3,k)*tempx1(i,1)^2));
                else
                    tempx2(i,1)=tempx1(i,1)-delta*(tempx1(i,1)-coeff(1,k)-coeff(2,k)*tempx1(i-1,1)-coeff(3,k)*tempx1(i-1,1)^2);
                end
            elseif quad==1
                if i==1
                    tempx2(i,1)=tempx1(i,1)-delta*(-(coeff(1,k)+2*coeff(2,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)*tempx1(i,1)-coeff(2,k)*tempx1(i,1)^2));
                elseif i<numS
                    tempx2(i,1)=tempx1(i,1)-delta*((tempx1(i,1)-coeff(1,k)*tempx1(i-1,1)-coeff(2,k)*tempx1(i-1,1)^2)-(coeff(1,k)+2*coeff(2,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)*tempx1(i,1)-coeff(2,k)*tempx1(i,1)^2));
                else
                    tempx2(i,1)=tempx1(i,1)-delta*(tempx1(i,1)-coeff(1,k)*tempx1(i-1,1)-coeff(2,k)*tempx1(i-1,1)^2);
                end
            else
                if i==1
                    tempx2(i,1)=tempx1(i,1)-delta*(-(coeff(1,k)-2*coeff(1,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)*tempx1(i,1)+coeff(1,k)*tempx1(i,1)^2));
                elseif i<numS
                    tempx2(i,1)=tempx1(i,1)-delta*((tempx1(i,1)-coeff(1,k)*tempx1(i-1,1)+coeff(1,k)*tempx1(i-1,1)^2)-(coeff(1,k)-2*coeff(1,k)*tempx1(i,1))*(tempx1(i+1,1)-coeff(1,k)*tempx1(i,1)+coeff(1,k)*tempx1(i,1)^2));
                else
                    tempx2(i,1)=tempx1(i,1)-delta*(tempx1(i,1)-coeff(1,k)*tempx1(i-1,1)+coeff(1,k)*tempx1(i-1,1)^2);
                end
            end                
        end
        L2=0;
        for i=1:numS-1
            if quad==2
                L2=L2+.5*(tempx2(i+1,1)-(coeff(1,k)+coeff(2,k)*tempx2(i,1)+coeff(3,k)*tempx2(i,1)^2))^2;
            elseif quad==1
                L2=L2+.5*(tempx2(i+1,1)-(coeff(1,k)*tempx2(i,1)+coeff(2,k)*tempx2(i,1)^2))^2;
            else
                L2=L2+.5*(tempx2(i+1,1)-(coeff(1,k)*tempx2(i,1)-coeff(1,k)*tempx2(i,1)^2))^2;
            end
        end
        if abs(L1-L2)<.00001
            stop=1;
        else
            tempx1=tempx2;
            L1=L2;
            tau=tau+1;
        end
    end
 
    %state update
    if MSA==0
        x(:,k+1)=.5*x(:,k)+.5*tempx2;
    else
        x(:,k+1)=(k/(k+1))*x(:,k)+(1/(k+1))*tempx2;
    end
    
    %parameter update
    if quad==2
        for i=2:numS
            xvar(i,1)=1;
            xvar(i,2)=x(i-1,k+1);
            xvar(i,3)=x(i-1,k+1)^2;
            yvar(i,1)=x(i,k+1);
        end
    elseif quad==1
        for i=2:numS
            xvar(i,1)=x(i-1,k+1);
            xvar(i,2)=x(i-1,k+1)^2;
            yvar(i,1)=x(i,k+1);
        end
    else
        for i=2:numS
            xvar(i,1)=x(i-1,k+1)*(1-x(i-1,k+1));
            yvar(i,1)=x(i,k+1);
        end
    end
    tempcoeff=(xvar'*xvar)^-1*xvar'*yvar;
    if MSA==0
        coeff(:,k+1)=.5*coeff(:,k)+.5*tempcoeff;
    else
        coeff(:,k+1)=(k/(k+1))*coeff(:,k)+(1/(k+1))*tempcoeff;
    end
    
    %stop criterion
    L(k+1)=0;
    for i=1:numS-1
        if quad==2
            L(k+1)=L(k+1)+.5*(x(i+1,k+1)-(coeff(1,k+1)+coeff(2,k+1)*x(i,k+1)+coeff(3,k+1)*x(i,k+1)^2))^2;
        elseif quad==1
            L(k+1)=L(k+1)+.5*(x(i+1,k+1)-(coeff(1,k+1)*x(i,k+1)+coeff(2,k+1)*x(i,k+1)^2))^2;
        else
            L(k+1)=L(k+1)+.5*(x(i+1,k+1)-(coeff(1,k+1)*x(i,k+1)-coeff(1,k+1)*x(i,k+1)^2))^2;
        end
    end
    if abs(L(k+1)-L(k))<tol
        stop1=1;
    else
        k=k+1;
    end
end