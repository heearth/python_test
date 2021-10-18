%% path������룺�͹۲���µķ���·������ê������STê����K����:�ޣ�

clear;clc;close all;
mu=0.05; r=0.03; s0=100; K=100; rho=[-0.1;0.1];
bs_sigma=0.3;
% kappa=2; theta=0.1; v0=0.1; zeta=-0.1; sigma=0.2; %����feller����
T=3; delta_t=1/250; tt=0:delta_t:T; steps=length(tt)-1; M=length(rho); aaa=inf; 
kappa=1; theta=0.01; v0=0.1; sigma=0.25; zeta=-0.1; %������feller������
% kappa=1; theta=0.5; v0=0.1; sigma=0.25; zeta=-0.1; %����feller������
z_s=zeros(length(rho),T/delta_t); z_v=zeros(length(rho),T/delta_t);
% mu_feedback=zeros(length(rho),T/delta_t);
% sigma_feedback=zeros(length(rho),T/delta_t);
%% ��ʼֵstart %%%%%%%%%%%%%%%%%%%%%%%%
s(:,1)=s0.*ones(length(rho),1); x(:,1)=log(s0).*ones(length(rho),1);
v(:,1)=v0.*ones(length(rho),1); y(:,1)=log(v0).*ones(length(rho),1);
tau(1)=T-0;

heston_s(:,1)=s0.*ones(length(rho),1); heston_x(:,1)=log(s0).*ones(length(rho),1); 
heston_v(:,1)=v0.*ones(length(rho),1);  heston_y(:,1)=log(v0).*ones(length(rho),1);

bsf_s(:,1)=s0.*ones(length(rho),1); bsf_x(:,1)=log(s0).*ones(length(rho),1);
%% t0ʱ��ϣ����ĸ��ʼֵ��i=0
i=0;
%% BS greeks
d1(:,i+1)=(log(s(:,i+1)./K)+(r+0.5*bs_sigma^2)*tau(i+1))/(bs_sigma*sqrt(tau(i+1)));
d2(:,i+1)=(log(s(:,i+1)./K)+(r-0.5*bs_sigma^2)*tau(i+1))/(bs_sigma*sqrt(tau(i+1)));
delta_c0(:,i+1)=normcdf(d1(:,i+1));
gamma_c0(:,i+1)=exp(-d1(:,i+1).^2/2)./(s(:,i+1)*bs_sigma*sqrt(2*pi*tau(i+1)));
charm_c0(:,i+1)=-exp(-d1(:,i+1).^2/2)/sqrt(2*pi).*(2*r*tau(i+1)-d2(:,i+1).*bs_sigma*sqrt(tau(i+1)))/(2*tau(i+1)*sigma*sqrt(tau(i+1)));
speed_c0(:,i+1)=-exp(-d1(:,i+1).^2/2)/sqrt(2*pi)./(s(:,i+1).^2*bs_sigma*sqrt(tau(i+1))).*...
    (d1(:,i+1)./(bs_sigma*sqrt(tau(i+1)))+1);
%% ����
delta(:,i+1)=delta_c0(:,i+1);
gamma(:,i+1)=gamma_c0(:,i+1);
charm(:,i+1)=-charm_c0(:,i+1);%delta��t��=-delta��tau��
speed(:,i+1)=speed_c0(:,i+1);

sigma_feedback(:,i+1)=sqrt(v(:,i+1))./(1-rho.*s(:,i+1).*gamma(:,i+1));
mu_feedback(:,i+1)=rho./(1-rho.*s(:,i+1).*gamma(:,i+1)).*...
    (charm(:,i+1)+v(:,i+1).*s(:,i+1).^2.*speed(:,i+1)./(2.*(1-rho.*s(:,i+1).*gamma(:,i+1)).^2))+...
    mu./(1-rho.*s(:,i+1).*gamma(:,i+1));
%% bsfģ��
bsf_sigma_feedback(:,i+1)=sigma/(1-rho.*s(:,i+1).*gamma(:,i+1));
bsf_mu_feedback(:,i+1)=rho./(1-rho.*s(:,i+1).*gamma(:,i+1)).*...
    (charm(:,i+1)+sigma^2*s(:,i+1).^2.*speed(:,i+1)./(2*(1-rho.*s(:,i+1).*gamma(:,i+1)).^2))+...
    mu./(1-rho.*s(:,i+1).*gamma(:,i+1));
%% ��ʼֵend %%%%%%%%%%%%%%%%%%%%%%%%

%% ����·�� %%%%%%%%%%%%%%%%%%%%%%%%
for i=1:steps
    %% feedback��c1:ŷ����ɢ
%     load z_s2 %rho=0.6
%    load z_s_neg %rho=-0.1
%    z_s=z_s.*ones(length(rho),1);
     z_s(:,i)=random('norm',0,1)*ones(length(rho),1);
    %����������steps���ļ�:save(strcat('C:\Users\Administrator.DESKTOP-5PPQTMV\Desktop\�о����׶�\Ͷ��\ϵͳ����ѧ��\�������\',num2str(:,i)),'z_s');
%     load ��z_v22222 %%rho=0.6
%    load z_v_neg %rho=-0.1
%    z_v=z_v.*ones(length(rho),1);
     z_v(:,i)=random('norm',0,1)*ones(length(rho),1);
    tau(i+1)=T-i*delta_t;
    %% hestonģ��·��
    %z1=random('norm',0,1,M,1);%����z�������ֲ�����С
    heston_x(:,i+1)=heston_x(:,i)+(mu-0.5*heston_v(:,i))*delta_t+sqrt(heston_v(:,i)).*z_s(:,i)*sqrt(delta_t);
    heston_xx(:,i+1)=heston_x(:,i+1)-heston_x(:,i); 
    heston_s(:,i+1)=exp(heston_x(:,i+1));
    heston_y(:,i+1)=heston_y(:,i)+(kappa*(theta-heston_v(:,i))./heston_v(:,i)-0.5*sigma^2*heston_v(:,i).^(-1))*delta_t+...
        sigma*heston_v(:,i).^(-1/2).*(zeta*z_s(:,i)+sqrt(1-zeta^2)*z_v(:,i))*sqrt(delta_t);
    heston_v(:,i+1)=exp(heston_y(:,i+1));
    %% =======================================================
    %v��ʼֵ��bs��gamma
    %�������Բ����
    %     x(:,i+1)=x(:,i)+(r-0.5*v(:,i)^2)*delta_t+v(:,i)*z*sqrt(delta_t);
    %��͹۲���¼۸����ߣ�b=ԭ��b/s
    x(:,i+1)=x(:,i)+(mu_feedback(:,i)-0.5.*sigma_feedback(:,i).^2).*delta_t+sigma_feedback(:,i).*z_s(:,i).*sqrt(delta_t);
    xx(:,i+1)=x(:,i+1)-x(:,i); s(:,i+1)=exp(x(:,i+1));
    y(:,i+1)=y(:,i)+(kappa*(theta-v(:,i))./v(:,i)-0.5*sigma^2*v(:,i).^(-1))*delta_t+...
        sigma*v(:,i).^(-1/2).*(zeta*z_s(:,i)+sqrt(1-zeta^2).*z_v(:,i))*sqrt(delta_t);
    v(:,i+1)=exp(y(:,i+1));
    %% ����ģ�����ݴ���bsϣ����ĸ��ʽ��c0
    %% ������Ȩ
    d1(:,i+1)=(log(s(:,i+1)./K)+(r+0.5*bs_sigma^2)*tau(i+1))/(bs_sigma*sqrt(tau(i+1)));
    d2(:,i+1)=(log(s(:,i+1)./K)+(r-0.5*bs_sigma^2)*tau(i+1))/(bs_sigma*sqrt(tau(i+1)));
    delta_c0(:,i+1)=normcdf(d1(:,i+1));
    gamma_c0(:,i+1)=exp(-d1(:,i+1).^2/2)./(s(:,i+1)*bs_sigma*sqrt(2*pi*tau(i+1)));
    charm_c0(:,i+1)=-exp(-d1(:,i+1).^2/2)/sqrt(2*pi).*(2*r*tau(i+1)-d2(:,i+1).*bs_sigma*sqrt(tau(i+1)))/(2*tau(i+1)*sigma*sqrt(tau(i+1)));
    speed_c0(:,i+1)=-exp(-d1(:,i+1).^2/2)/sqrt(2*pi)./(s(:,i+1).^2*bs_sigma*sqrt(tau(i+1))).*...
        (d1(:,i+1)./(bs_sigma*sqrt(tau(i+1)))+1);
    %% ����ϣ����ĸ
    delta(:,i+1)=delta_c0(:,i+1);
    gamma(:,i+1)=gamma_c0(:,i+1);
    charm(:,i+1)=-charm_c0(:,i+1);%delta��t��=-delta��tau��
    speed(:,i+1)=speed_c0(:,i+1);
    %% �����ʡ�Ư����
    %% svfģ��
    sigma_feedback(:,i+1)=sqrt(v(:,i+1))./(1-rho.*s(:,i+1).*gamma(:,i+1));
    mu_feedback(:,i+1)=rho./(1-rho.*s(:,i+1).*gamma(:,i+1)).*...
        (charm(:,i+1)+v(:,i+1).*s(:,i+1).^2.*speed(:,i+1)./(2*(1-rho.*s(:,i+1).*gamma(:,i+1)).^2))+...
        mu./(1-rho.*s(:,i+1).*gamma(:,i+1));
    %% =======================================================
    %% bsfģ��
    bsf_x(:,i+1)=bsf_x(:,i)+(bsf_mu_feedback(:,i)-0.5.*bsf_sigma_feedback(:,i).^2).*delta_t+bsf_sigma_feedback(:,i).*z_s(:,i).*sqrt(delta_t);
    bsf_xx(:,i+1)=bsf_x(:,i+1)-bsf_x(:,i);
    bsf_s(:,i+1)=exp(bsf_x(:,i+1));
    %% ����ģ�����ݴ���bsϣ����ĸ��ʽ��c0
    bsf_d1(:,i+1)=(log(bsf_s(:,i+1)./K)+(r+0.5*sigma^2)*tau(i+1))./(sigma*sqrt(tau(i+1)));
    bsf_d2(:,i+1)=(log(bsf_s(:,i+1)./K)+(r-0.5*sigma^2)*tau(i+1))./(sigma*sqrt(tau(i+1)));
    %% ������Ȩ
    %call
    bsf_delta_c0(:,i+1)=normcdf(bsf_d1(:,i+1));
    %put
    bsf_delta_c0(:,i+1)=-normcdf(-bsf_d1(:,i+1));
    %ͬ
    bsf_gamma_c0(:,i+1)=exp(-bsf_d1(:,i+1).^2./2)./(bsf_s(:,i+1)*sigma*sqrt(2*pi*tau(i+1)));
    bsf_charm_c0(:,i+1)=-exp(-bsf_d1(:,i+1).^2./2)./sqrt(2*pi).*...
        (2*r*tau(i+1)-bsf_d2(:,i+1).*sigma.*sqrt(tau(i+1)))./(2*tau(i+1)*sigma*sqrt(tau(i+1)));
    bsf_speed_c0(:,i+1)=-exp(-bsf_d1(:,i+1).^2./2)./sqrt(2*pi)./(bsf_s(:,i+1).^2*sigma*sqrt(tau(i+1))).*...
        (bsf_d1(:,i+1)./(sigma*sqrt(tau(i+1)))+1);
    %% ������Ȩ
    %call
    bsf_delta_c0(:,i+1)=-normcdf(bsf_d1(:,i+1));
    %put
    bsf_delta_c0(:,i+1)=normcdf(-bsf_d1(:,i+1));
    %ͬ
    bsf_gamma_c0(:,i+1)=-(exp(-bsf_d1(:,i+1).^2/2)./(bsf_s(:,i+1).*sigma.*sqrt(2*pi*tau(i+1))));
    bsf_charm_c0(:,i+1)=-(-exp(-bsf_d1(:,i+1).^2/2)./sqrt(2*pi).*...
        (2*r*tau(i+1)-bsf_d2(:,i+1).*sigma.*sqrt(tau(i+1)))./(2*tau(i+1)*sigma*sqrt(tau(i+1))));
    bsf_speed_c0(:,i+1)=-(-exp(-bsf_d1(:,i+1).^2/2)./sqrt(2*pi)./(bsf_s(:,i+1).^2*sigma*sqrt(tau(i+1))).*...
        (bsf_d1(:,i+1)./(sigma*sqrt(tau(i+1)))+1));
    %% ����ϣ����ĸ
    bsf_delta(:,i+1)=bsf_delta_c0(:,i+1);
    bsf_gamma(:,i+1)=bsf_gamma_c0(:,i+1);
    bsf_charm(:,i+1)=-bsf_charm_c0(:,i+1);%delta��t��=-delta��tau��
    bsf_speed(:,i+1)=bsf_speed_c0(:,i+1);
    %% �����ʡ�Ư����
    bsf_sigma_feedback(:,i+1)=sigma./(1-rho.*bsf_s(:,i+1).*bsf_gamma(:,i+1));
    bsf_mu_feedback(:,i+1)=rho./(1-rho.*bsf_s(:,i+1).*bsf_gamma(:,i+1)).*...
        (bsf_charm(:,i+1)+sigma^2.*bsf_s(:,i+1).^2.*bsf_speed(:,i+1)./(2.*(1-rho.*bsf_s(:,i+1).*bsf_gamma(:,i+1)).^2))+...
        mu./(1-rho.*bsf_s(:,i+1).*bsf_gamma(:,i+1));
end

%% ������������ʵ�ƫ�ȷ�� %%%%%%%%%%%%%%%%%%%%%%%%
feedback_kurtois = kurtosis(xx,1,2);                    %���Ͷ�=���
feedback_sknew = skewness(xx,1,2);                      %ƫб��=ƫ��
heston_kurtois = kurtosis(heston_xx,1,2);
heston_sknew = skewness(heston_xx,1,2);
bsf_kurtois = kurtosis(bsf_xx,1,2);                    %���Ͷ�=���
bsf_sknew = skewness(bsf_xx,1,2);
a=["ģ��" "rho=-0.1�ķ��" "rho=-0.1��ƫ��" "rho=0.1�ķ��" "rho=0.1��ƫ��"
    "SVF" feedback_kurtois(1) feedback_sknew(1) feedback_kurtois(2) feedback_sknew(2)
    "HESTON" heston_kurtois(1)  heston_sknew(1) heston_kurtois(2) heston_sknew(2)
    "GBMF" bsf_kurtois(1) bsf_sknew(1) bsf_kurtois(2) bsf_sknew(2)];
%% ��ͼ %%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot ( tt , s,'b--');hold on;
% plot ( tt , heston_s,'r--');hold on;
% plot ([0,T],[K,K]);%������
% % hold on;plot ([aaa*delta_t,aaa*delta_t],[min(s),max(s)]);%������
% xlabel ( 't' ) ;ylabel ( 'S' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;axis tight;
% legend('������ģ��','hestonģ��');

% figure
% plot ( tt , x,'b--');hold on;
% plot ( tt , bs_x,'r--');hold on;
% plot ([0,T],[log(K),log(K)]);
% xlabel ( 't' ) ;ylabel ( 'X' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;axis tight;
% legend('������ģ��','BSģ��');

% figure
% plot ( tt+delta_t , xx,'b--');hold on;
% % plot ( tt+delta_t , bs_xx,'r--');hold on;
% plot ([0,T],[0,0]);
% xlabel ( 'ʱ��' ) ;ylabel ( '����������' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;axis tight;
% % legend('������ģ��','BSģ��');


% figure
% % subplot(2,2,1);plot (tt,b);xlabel ( 't' ) ;ylabel ( 'b') ;axis tight;
% subplot(2,2,1);
% plot (tt,delta);
% % hold on;plot (tt,rho*delta_c11,'r-');
% % hold on;plot (tt,delta_c0,'b-');
% xlabel ( 't' ) ;ylabel ( 'delta') ;axis tight;
% subplot(2,2,2);
% plot (tt,gamma);
% % hold on;plot (tt,rho*gamma_c1,'r-');
% % hold on;plot (tt,gamma_c0,'b-');
% xlabel ( 't' ) ;ylabel ( 'gamma') ;axis tight;
% subplot(2,2,3);
% plot (tt,charm);
% % hold on;plot (tt,rho*charm_c1,'r-');
% % hold on;plot (tt,-charm_c0,'b-');
% xlabel ( 't' ) ;ylabel ( 'charm') ;axis tight;
% subplot(2,2,4);
% plot (tt,speed);
% % hold on;plot (tt,rho*speed_c1,'r-');
% % hold on;plot (tt,speed_c0,'b-');
% xlabel ( 't' ) ;ylabel ( 'speed') ;axis tight;
