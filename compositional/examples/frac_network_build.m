%%建立裂隙网络模型——《基于蒙特卡洛法的岩体裂隙网络模型及渗透率张量的研究》——哈尔滨工业大学
% 可改造成函数，由用户直接在主程序中定义输入值，之后调用函数

function[LXDD] = frac_network_build(N,chang,kuan,midu,zouxjz,zouxbzc,jicjz,jicbzc,xikanjz,xikanbzc)
%N=input('裂隙组数=');
%chang=input('生成域长=');
%kuan=input('生成域宽=');
mianji=chang*kuan; %面积
LXDD=[];
for k=1:N
%midu=input('密度=');  %密度的取值范围，是如何进行标定的？
%zouxjz=input('走向均值=');
%zouxbzc=input('走向标准差=');
%jicjz=input('迹长均值=');
%jicbzc=input('迹长标准差=');
%xikanjz=input('隙宽均值=');
%xikanbzc=input('隙宽标准差=');
n=mianji*midu; %条数
zhongxidian=unifrnd(0,chang,n,2); %中心点 (-chang/2,chang/2,n,2)改为(0,chang,n,2)
zouxiang=normrnd(zouxjz,zouxbzc,n,1); %走向  ~可改为文献中的参数类型
m=jicjz;
v=jicbzc^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
jichang=lognrnd(mu,sigma,n,1); %迹长
xikuang=normrnd(xikanjz,xikanbzc,n,1); %隙宽
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
%LXDD1=[x1,y1,x2,y2,zouxiang,xikuang]; 可以将zouxiang和xikuang去掉后，再进行矩阵的组装
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
% 试验数据1 frac_network_build(1,400,400,0.001,45,20,25,10,15,5)
% 试验数据2 frac_network_build(1,1200,600,0.001,45,20,25,10,15,5)