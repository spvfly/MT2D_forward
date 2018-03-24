clear;
clc;
%geo
X_length = 1000;
Y_length = -1000;
NXTM = 55; % X����  ������ �����
NYTM = 50; % Y����  ������ �����
NETM = NXTM * NYTM; % �ܵ�Ԫ����
NPSUMTM = (NXTM + 1) * (NYTM + 1); %�ܽڵ����
a = 100; % ��Ԫ��  ���򳤶�
h = 100; % % ��Ԫ��  ���򳤶�
PaTM = [100, 1500]; % ��Ҫ�趨�ĵ�����ֵ
f=logspace(3, -3, 20); % Ƶ��
xTM = zeros(NPSUMTM, length(f));
X = linspace(0, 1000, NXTM+1);
Y = linspace(0, -1000, NYTM+1);
XY = zeros(2, NPSUMTM); %�ڵ�������
for ix = 1:NXTM+1
    for iy = 1:NYTM+1
        N = (ix-1)*(NYTM+1) + iy;
        XY(1, N) = X(ix);
        XY(2, N) = Y(iy);
    end
end

BHXY = zeros(4, NETM); %��Ԫ��  �ڵ�������
for IX=1:NXTM
    for IY=1:NYTM
        N=(IX-1)*NYTM+IY;%��Ԫ���
        N1=N+(IX-1);%ÿ����Ԫ�ĵ�һ���ڵ���
        BHXY(1,N)=N1;%ÿ����Ԫ�ĵ�һ���ڵ���
        BHXY(2,N)=N1+1;%ÿ����Ԫ�ĵڶ����ڵ���
        BHXY(3,N)=N1+NYTM+2;%ÿ����Ԫ�ĵ������ڵ���
        BHXY(4,N)=N1+NYTM+1;%ÿ����Ԫ�ĵ��ĸ��ڵ���
    end
end
%%%%%%%�����ص�ģ��%%%%%%
rhoTM = zeros(NETM, 1); % ��Ԫ�� ��������Ϣ
% ����Ϊ ����ص�ģ�� �ʲ�״�ṹ
%  for i=1:30
%      for j=1:NYTM
%          m=(i-1)*NYTM+j;%��Ԫ���
%          rhoTM(m)=PaTM(1);
%      end
%  end
%  for i=31:60
%      for j=1:NYTM
%          m=(i-1)*NYTM+j;%��Ԫ���
%          rhoTM(m)=PaTM(2);
%      end
%  end
%  for i=1:NXTM
%      for j=1:NYTM
%          m=(i-1)*NYTM+j;
%          rhoTM(m)=100;
%      end
%  end
% for i=16:24%�쳣��
%      for j=2:18
%          m=(i-1)*NYTM+j;
%          rhoTM(m)=400;
%      end
% end
% for i=26:30%�쳣��
%      for j=1:10
%          m=(i-1)*NYTM+j;
%          rhoTM(m)=800;
%      end
%  end
for i = 1:NXTM
    for j = 1:5
        m=(i-1)*NYTM+j;
        rhoTM(m) = 10;
    end
end
for i = 1:NXTM
    for j = 6:15
        m=(i-1)*NYTM+j;
        rhoTM(m) = 100;
    end
end
for i = 1:NXTM
    for j = 16:40
        m=(i-1)*NYTM+j;
        rhoTM(m) = 10;
    end
end
for i = 1:NXTM
    for j = 41:NYTM
        m=(i-1)*NYTM+j;
        rhoTM(m) = 100;
    end
end

        
 %%%%%%%%%%%%%%
 for ff=1:size(f,2)
     K1TM=zeros(NPSUMTM); %
     K2TM=zeros(NPSUMTM); %
     K3TM=zeros(NPSUMTM); %
     PTM=zeros(NPSUMTM,1); %
     u=(4e-7)*pi;%�硢���Բ���
     w=2*pi*f(ff);
     TAOTM=rhoTM;
     LMTM=sqrt(-1)*w*u;
     %%%% K1
     for m=1:NETM
         A=(h/a)/6;
         B=(a/h)/6;
         K1e=A*[2 1 -1 -2;1 2 -2 -1;-1 -2 2 1;-2 -1 1 2]+B*[2 -2 -1 1;-2 2 1 -1;-1 1 2 -2;1 -1 -2 2];   
         for j=1:4
             NJ=BHXY(j,m);
             for k=1:4
                 NK=BHXY(k,m);
                 K1TM(NJ,NK)=K1TM(NJ,NK)+K1e(j,k)*TAOTM(m);
             end
         end
     end
      %%%% K2
      for m=1:NETM
           K2e=[4 2 1 2;2 4 2 1;1 2 4 2;2 1 2 4];
           for j=1:4
               NJ=BHXY(j,m);
               for k=1:4
                   NK=BHXY(k,m);
                   K2TM(NJ,NK)=K2TM(NJ,NK)+K2e(j,k)*LMTM*a*h/36;
               end
           end
      end
       %%%%%  K3
       for n1=NYTM:NYTM:NETM
           i=BHXY(1,n1); j=BHXY(2,n1);
           kn=sqrt(-sqrt(-1)*w*u/rhoTM(n1));
           mk=TAOTM(n1)*kn*a/6;%???????????????????????????
           Kjj=2*mk;
           K3TM(j,j)=K3TM(j,j)+Kjj;
           K3TM(i,i)=K3TM(i,i)+Kjj;
           Kji=1*mk;
           K3TM(j,i)=K3TM(j,i)+Kji;
           K3TM(i,j)=K3TM(i,j)+Kji;
       end
       %%%%%  ��װ����նȾ���
       vTM=sparse(K1TM-K2TM+K3TM);
       
       PTM=zeros(NPSUMTM,1);
       for i=1:NXTM+1
           j=1+(i-1)*(NYTM+1);
           vTM(j,j)=vTM(j,j)*10^10;
           PTM(j)=vTM(j,j)*1;
       end
       vvTM=vTM;
       %xTM=zeros(size(PTM),2);
       tol=1e-15;
       maxsteps=100;
       [LTM,UTM]=luinc(vvTM,'0');
       xTM(:,ff)=bicgstab(vvTM,PTM,tol,maxsteps,LTM,UTM);
 end
 %%%%%%
vect=1:NYTM+1:NPSUMTM-NYTM;%���ǵر�ڵ�ı��
% du/dy = (-11*u1 + 18*u2 - 9u3 + 2u4)
u1=zeros(size(vect,2),size(f,2));
u2=zeros(size(vect,2),size(f,2));
u3=zeros(size(vect,2),size(f,2));
u4=zeros(size(vect,2),size(f,2));
z=zeros(size(vect,2),size(f,2));     %%�迹
pc=zeros(size(vect,2),size(f,2));    %�ӵ�����
phase=zeros(size(vect,2),size(f,2)); %�迹��λ
ve=[1:NYTM:(NYTM*NXTM-NYTM+1),NYTM*NXTM-NYTM+1];
%���һ����Ԫ��Ϊ�˲�������ve�ĳ�����vect��ͬ������������
% ��ȡ�ر�� �迹 �ӵ����� �迹��λ 
for ff=1:size(f,2)
       for i=1:size(vect,2)
           u=(4e-7)*pi;
           w=2*pi*f(ff);   
           u1(i,ff)=xTM(vect(i),ff);
           u2(i,ff)=xTM(vect(i)+1,ff);
           u3(i,ff)=xTM(vect(i)+2,ff);
           u4(i,ff)=xTM(vect(i)+3,ff);
         %�ĸ��ڵ�ı߳�
           l=3*h;
           ux(i,ff)=(-11*u1(i,ff)+18*u2(i,ff)-9*u3(i,ff)+2*u4(i,ff))/(2*l); %��������
           z(i,ff)=(rhoTM(ve(i))*ux(i,ff))/xTM(vect(i),ff);%�迹
           pc(i,ff)=(abs((rhoTM(ve(i))*ux(i,ff))/xTM(vect(i),ff))^2)/(w*u);    %�ӵ�����
           phase(i,ff)=-atan(imag(z(i,ff))/real(z(i,ff)))*180/pi;                %�迹��λ
           HEY(i,ff)=ux(i,ff)*rhoTM(ve(i));%EY
       end  
end
% semilogx��x1��y1��ѡ��1��x2��y2��ѡ��2������
% plot(f,pc(20,:));
semilogx(1./f, pc(20, :));

       
