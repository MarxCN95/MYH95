%%������϶����ģ�͡������������ؿ��巨��������϶����ģ�ͼ���͸���������о���������������ҵ��ѧ
% �ɸ���ɺ��������û�ֱ�����������ж�������ֵ��֮����ú���

function[LXDD] = frac_network_build(N,chang,kuan,midu,zouxjz,zouxbzc,jicjz,jicbzc,xikanjz,xikanbzc)
%N=input('��϶����=');
%chang=input('������=');
%kuan=input('�������=');
mianji=chang*kuan; %���
LXDD=[];
for k=1:N
%midu=input('�ܶ�=');  %�ܶȵ�ȡֵ��Χ������ν��б궨�ģ�
%zouxjz=input('�����ֵ=');
%zouxbzc=input('�����׼��=');
%jicjz=input('������ֵ=');
%jicbzc=input('������׼��=');
%xikanjz=input('϶���ֵ=');
%xikanbzc=input('϶���׼��=');
n=mianji*midu; %����
zhongxidian=unifrnd(0,chang,n,2); %���ĵ� (-chang/2,chang/2,n,2)��Ϊ(0,chang,n,2)
zouxiang=normrnd(zouxjz,zouxbzc,n,1); %����  ~�ɸ�Ϊ�����еĲ�������
m=jicjz;
v=jicbzc^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
jichang=lognrnd(mu,sigma,n,1); %����
xikuang=normrnd(xikanjz,xikanbzc,n,1); %϶��
x1=zhongxidian(:,1)+(jichang./2).*cos(zouxiang*pi/180);
y1=zhongxidian(:,2)+(jichang./2).*sin(zouxiang*pi/180);
x2=zhongxidian(:,1)+(jichang./2).*cos((zouxiang+180)*pi/180);
y2=zhongxidian(:,2)+(jichang./2).*sin((zouxiang+180)*pi/180);
for i=1:n
if x1(i)>x2(i)

w=x2(i);
x2(i)=x1(i);
q=y2(i);
y2(i)=y1(i);
x1(i)=w;
y1(i)=q;
end
end
%LXDD1=[x1,y1,x2,y2,zouxiang,xikuang]; ���Խ�zouxiang��xikuangȥ�����ٽ��о������װ
LXDD1=[x1,y1,x2,y2];
if isempty(LXDD)
LXDD=LXDD1;
else
LXDD=[LXDD;LXDD1];
end
for i=1:n
plot([x1(i),x2(i)],[y1(i),y2(i)]);
xlim([-chang/2,chang/2]);
ylim([-kuan/2,kuan/2]);
hold on
end
end



%model.transfer_model_object = LinearSaturationTransferFunction(1e-9,0.8,0.3,act);
% ��������1 frac_network_build(1,400,400,0.001,45,20,25,10,15,5)
% ��������2 frac_network_build(1,1200,600,0.001,45,20,25,10,15,5)