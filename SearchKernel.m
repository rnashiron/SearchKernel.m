%特性関数の入力 f(v12,v13,v23,vN)
%example
f(1,1,0,1)
f(1,1,1,1)
f(6,0,8,20)

f(-4,-3,-2,10)
f(-4,-3,5,10)
f(-4,2,5,10)

f(2,6,18,10)
f(2,15, 18,10)
f(12,15,18,10)


function f(v12,v13,v23,vN)
%---------特性関数を任意に設定してください---------
v(3)=v12; %v(12)
v(4)=v13; %v(13)
v(5)=v23; %v(23)
v(6)=vN; %v(N)
%------------------------------------------------

v(1)=0; %v(0)=0,fixed
v(2)=0; %v(1)=v(2)=v(3)=0,fixed
v

p=1000; %分割の量
[a,b]=meshgrid(0:vN/p:vN);
C12=max(v(4)-vN+b,-a)-max(v(5)-vN+a,-b); %s12-s21
C13=max(v(3)-a-b ,-a)-max(v(5)-vN+a,-vN+a+b); %s13-s31
C23=max(v(3)-b-a ,-b)-max(v(4)-vN+b,-vN+a+b); %s23-s32


D12=zeros(p+1,p+1);
D13=zeros(p+1,p+1);
D23=zeros(p+1,p+1);

for i=1:p+1
    for j=1:p+1
        if i+j>p+2
            a(i,j)=0;
            b(i,j)=0;
            C12(i,j)=NaN;
            C13(i,j)=NaN;
            C23(i,j)=NaN;
        end
    end    
end


e=0.0001; %閾値
for i=1:p+1
    for j=1:p+1
        if abs(C12(i,j))<e
            D12(i,j)=1;
        end
        if abs(C13(i,j))<e
            D13(i,j)=1;
        end
        if abs(C23(i,j))<e
            D23(i,j)=1;
        end
    end
end

for i=1:p+1
    %a(i,1) %x1=0
    if C12(i,1)<0 %s12<s21
        D12(i,1)=1;
    end
    if C13(i,1)<0 %s13<s31
        D13(i,1)=1;
    end
    %b(1,i) %x2=0
    if C12(1,i)>0 %s12>s21
        D12(1,i)=1;
    end
    if C23(1,i)<0 %s23<s32
        D23(1,i)=1;
    end
    %vN-a(i,p-i+2)-b(i,p-i+2) %x3=0
    if C13(i,p-i+2)>0 %s13>s31
        D13(i,p-i+2)=1;
    end
    if C23(i,p-i+2)>0 %s23>s32
        D23(i,p-i+2)=1;
    end
end

%{
C12
D12
C13
D13
C23
D23
%}

S12x=[];
S12y=[];
S13x=[];
S13y=[];
S23x=[];
S23y=[];

%plotのための準備
for i=1:p+1
    for j=1:p+1
        if D12(i,j)==1;
            S12x=[S12x,(j-1)*vN/p];
            S12y=[S12y,(i-1)*vN/p];
        end
        if D13(i,j)==1;
            S13x=[S13x,(j-1)*vN/p];
            S13y=[S13y,(i-1)*vN/p];
        end
        if D23(i,j)==1;
            S23x=[S23x,(j-1)*vN/p];
            S23y=[S23y,(i-1)*vN/p];
        end
    end
end

x=-0.05*vN:vN/100:1.05*vN;
figure
hold on
scatter(S13x,S13y,25,'g','filled')
scatter(S23x,S23y,10,'b','filled')
scatter(S12x,S12y,6,'r','filled')
plot(x,vN-x)
hold off
title(['Kernel: N=3, v(12) = ',num2str(v(3)),...
    ', v(13) = ',num2str(v(4)),', v(23) = ',num2str(v(5)),', v(N) = ',num2str(v(6))])
xlabel('x_1')
ylabel('x_2')
xlim([-0.05*vN,1.05*vN])
ylim([-0.05*vN,1.05*vN])
legend('x_1\sim_x x_3','x_2\sim_x x_3','x_1\sim_x x_2')


Kx=[];
Ky=[];
for i=1:p+1
    for j=1:p+1
        if D12(i,j)==1&D13(i,j)==1&D23(i,j)==1
            Kx=[Kx,(j-1)*vN/p];
            Ky=[Ky,(i-1)*vN/p];
        end
    end
end
if ~isempty(Kx)
    disp('カーネルは')
    for i=1:length(Kx)
        disp(['(',num2str(Kx(i)),',',num2str(Ky(i)),',',num2str(vN-Kx(i)-Ky(i)),')'])
    end
    disp(['の',num2str(length(Kx)),'点です。'])
else
    disp('解が見つからなかったので近似解を探します')
end
    
%解が求まらなかった場合の近似解の導出
if isempty(Kx)
    for i=2:p
        for j=2:p
            if D12(i-1,j-1)==1||D12(i-1,j)==1||D12(i-1,j+1)==1||D12(i,j-1)==1||D12(i,j)==1||D12(i,j+1)==1||D12(i+1,j-1)==1||D12(i+1,j)==1||D12(i+1,j+1)==1
                if D13(i-1,j-1)==1||D13(i-1,j)==1||D13(i-1,j+1)==1||D13(i,j-1)==1||D13(i,j)==1||D13(i,j+1)==1||D13(i+1,j-1)==1||D13(i+1,j)==1||D13(i+1,j+1)==1
                    if D23(i-1,j-1)==1||D23(i-1,j)==1||D23(i-1,j+1)==1||D23(i,j-1)==1||D23(i,j)==1||D23(i,j+1)==1||D23(i+1,j-1)==1||D23(i+1,j)==1||D23(i+1,j+1)==1
                        Kx=[Kx,(j-1)*vN/p];
                        Ky=[Ky,(i-1)*vN/p];
                    end
                end
            end
        end
    end
    if ~isempty(Kx)
        disp('カーネルの近似解は')
        disp(['(',num2str(sum(Kx)/length(Kx)),',',num2str(sum(Ky)/length(Ky)),',',num2str(vN-sum(Kx)/length(Kx)-sum(Ky)/length(Ky)),')'])
        disp('です。')
    else
        disp('解が求まりませんでした。特性関数やp,eの値を確認してください')
    end
end

end
