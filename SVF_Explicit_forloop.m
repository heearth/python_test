%% implicit finite difference scheme, for call option
clear;close all;clc;
%params set   1/12=0.0833
% sigma=0.4;  K=100;  r=0.03;  q=0;  Smax=300;  T=1/12;      rho=0.005; ds=1.875;  dtau=5.21e-04;  sigma2=sigma.*sigma;
%% ��ϵ��
% r=0.03; s0=10; K=10; q=0; kappa=0.5; theta=0.16; v0=0.0625; sigma=0.9; zeta=0.1;
% r=0.02; s0=100; K=100; q=0; kappa=1.5; theta=0.04; v0=0.05; sigma=0.3; zeta=-0.9;
r=0.03; s0=100; K=85; q=0; 
kappa=1; theta=0.01; v0=0.1; sigma=0.25; zeta=-0.1;%������������������ȷbut������feller����
% kappa=2; theta=0.05; v0=0.1; zeta=-0.1; sigma=0.1; %����feller����
%hestonģ��feller������2*kappa*theta>sigma^2
rho=-0.05; %����ϵ��
%Ϊ+0.01*��λ��ʱ����Ȩ�۸�Ϊ��������hestonģ����ͬ=38��
%Ϊ-0.01*��λ��ʱ����Ȩ�۸�Ϊ������=38��
% %��¶��ʦ���������
% r=0.1; s0=10; K=10; q=0; kappa=0.5; theta=0.16; v0=0.0625; sigma=0.9; zeta=0.1;
% rho=0.2; %�г�������ϵ��
%% ��������
Smin=0; Smax=200; ds=10;%0.1��10
Vmin=0; Vmax=1; dv=0.055;%85��0.035������0.055
%��¶��ʦ���������
% Smin=6; Smax=570; ds=10;
% Vmin=0.0263; Vmax=1; dv=0.1;
Tmin=0; T=1; dt=1/250;
% Tmin=0; T=0.1; dt=1/250;
N=round(T/dt);  %ʱ�� N
M=round(Smax/ds); %�ռ� M��S
V=round(Vmax/dv); %�ռ� V��V
ds=Smax/M; dt=T/N;
vetS=linspace(Smin,Smax,M+1);%vetS=ds.*veti;
vetV=linspace(Vmin,Vmax,V+1);%vetV=dV.*vetj;
vetn=0:N;  %ʱ�� n
veti=0:M;  %�ռ� i:S
vetj=0:V;  %�ռ� j:V
%% ��ʼ��(S,V,t)
c=zeros(M+1,V+1,N+1);  %��ʼ����ά����S��V��Tҳ
%% �߽�����
% ��V=0ʱ��һ��ʽ��=0����������Ȩ�۸���ڵ����֮�󣬼���˱߽�����,���ڵڶ������֮ǰ
c(:,end,:)=repmat(vetS',1,N+1);              % ��VΪVmaxʱ����Ȩ�۸�=S
c(1,:,:)=0;                                 % ��S=0ʱ����Ȩ�۸�
c(end,:,:)=repmat(Smax-K.*exp(-r.*(T-dt.*vetn)),V+1,1);%��S=Smaxʱ����Ȩ�۸� = ��SΪ�����ʱ����Ȩ�Լ۸��һ�׵�=1
c(:,:,end)=repmat(max(vetS'-K,0),1,V+1);    % ��ʱ�䵽��ʱ����Ȩ�۸�t=T����ӦN+1λ��
%% ==============================
%% ����Ŀ�ģ����tau=T��t=0����Ȩ�۸��棬ʱ���λ����1����c(:,:,1)
%% ==============================

%% ------------------------------ʹ��SVF��ϣ����ĸ����������г�
%% �ڵ��ϵ���������
for n=N+1:-1:2 %�ƽ���һ��ʱ���
    for i=2:M %��T-1��V=0�ı߽�
        for j=2:V %�����ڵ�=��Ȩ�۸�
            %��n+1��DerSS�ɵ�n+1��ĵ����
            DerSS=(c(i+1,j,n) - 2*c(i,j,n) + c(i-1,j,n))  ./ ds.^2;
            aaa(i,j,n)=DerSS;
            Volatility_VARYING=sqrt(vetV(j)) / ( 1-( rho.* vetS(i)*DerSS ) );
            aaaaaa(i,j,n)=Volatility_VARYING;
            Volatility_1=sqrt(vetV(j)) * Volatility_VARYING;
            Volatility_2=Volatility_VARYING * Volatility_VARYING;
            %ϵ������ֵ���Ǿ������¸�ֵ
            A = zeta*sigma*Volatility_1.*vetS(i) / ( 4*ds*dv );
            B = Volatility_2*( vetS(i).*vetS(i) ) / ( 2*ds^2 );
            C = -r*vetS(i) / ( 2*ds );
            AA = sigma^2*vetV(j) / ( 2*dv^2 );
            CC = -kappa*(theta-vetV(j)) / ( 2*dv );
            BB = - ( B*2 +  AA*2 + r - 1/dt );
            BC1 = B+C; BC2 = B-C;
            AACC1 = AA+CC; AACC2 = AA-CC;
            c(i,j,n-1)=dt*(   A*c(i-1,j-1,n) + BC1*c(i-1,j,n) + (-A)*c(i-1,j+1,n)+...
                AACC1*c(i,j-1,n) + BB*c(i,j,n) + AACC2*c(i,j+1,n)+...
                (-A)*c(i+1,j-1,n) + BC2*c(i+1,j,n) + A*c(i+1,j+1,n)   );
        end
        % ��V=0ʱ��һ��ʽ��=0����������Ȩ�۸���ڵ����֮�󣬼���˱߽�����,���ڵڶ������֮ǰ
        c(i,1,n-1)=( r*vetS(i)/ds + kappa*theta/dv + r + 1/dt )^-1 * ...
            ( r*vetS(i)/ds*c(i+1,1,n-1) + kappa*theta/dv*c(i,2,n-1) + 1/dt*c(i,1,n) );
    end
end
%% ------------------------------ʹ��SVF��ϣ����ĸ����������г�

%% ------------------------------ʹ��BS��ϣ����ĸ���ǳ�������г�
% %% �ڵ��ϵ���������
% for n=N+1:-1:2 %�ƽ���һ��ʱ���
%     for i=2:M %��T-1��V=0�ı߽�
%         for j=2:V %�����ڵ�=��Ȩ�۸�
%             bs_sigma=0.3; tau=T-(n-2)*dt;%n=N+1����1��Ҫ��2
%             d1=(log(vetS(i)/K)+(r+0.5*bs_sigma^2)*tau)/(bs_sigma*sqrt(tau));
%             gamma=exp(-d1^2/2)/(vetS(i)*bs_sigma*sqrt(2*pi*tau));
%             DerSS=gamma;
%             bbb(i,j,n)=DerSS;
%             Volatility_VARYING=sqrt(vetV(j)) / ( 1-( rho.* vetS(i)*DerSS ) );
%             aaaaaa(i,j,n)=Volatility_VARYING;
%             Volatility_1=sqrt(vetV(j)) * Volatility_VARYING;
%             Volatility_2=Volatility_VARYING * Volatility_VARYING;
%             %ϵ������ֵ���Ǿ������¸�ֵ
%             A = zeta*sigma*Volatility_1.*vetS(i) / ( 4*ds*dv );
%             B = Volatility_2*( vetS(i).*vetS(i) ) / ( 2*ds^2 );
%             C = -r*vetS(i) / ( 2*ds );
%             AA = sigma^2*vetV(j) / ( 2*dv^2 );
%             CC = -kappa*(theta-vetV(j)) / ( 2*dv );
%             BB = - ( B*2 +  AA*2 + r - 1/dt );
%             BC1 = B+C; BC2 = B-C;
%             AACC1 = AA+CC; AACC2 = AA-CC;
%             c(i,j,n-1)=dt*(   A*c(i-1,j-1,n) + BC1*c(i-1,j,n) + (-A)*c(i-1,j+1,n)+...
%                 AACC1*c(i,j-1,n) + BB*c(i,j,n) + AACC2*c(i,j+1,n)+...
%                 (-A)*c(i+1,j-1,n) + BC2*c(i+1,j,n) + A*c(i+1,j+1,n)   );
%         end
%         % ��V=0ʱ��һ��ʽ��=0����������Ȩ�۸���ڵ����֮�󣬼���˱߽�����,���ڵڶ������֮ǰ
%         c(i,1,n-1)=( r*vetS(i)/ds + kappa*theta/dv + r + 1/dt )^-1 * ...
%             ( r*vetS(i)/ds*c(i+1,1,n-1) + kappa*theta/dv*c(i,2,n-1) + 1/dt*c(i,1,n) );
%     end
% end
%% ------------------------------ʹ��BS��ϣ����ĸ���ǳ�������г�

%% Calculate the option price------ʹ��SVF��ϣ����ĸ����������г�
option_price=interp2(vetS,vetV,c(:,:,1)',s0,v0)
DSS=interp2(vetS(1:end-1),vetV(1:end-1),aaa(:,:,251)',s0,v0);
% DSS=interp2(vetS(1:end-1),vetV(1:end-1),aaa(:,:,26)',s0,v0);
% sum(sum(c(:,:,1)<0))
% s0
% K
IV=blsimpv(s0, K, r, T, option_price)
%ʦ�㣺price =18.2261
%% Calculate the option price------ʹ��BS��ϣ����ĸ���ǳ�������г�
% option_price=interp2(vetS,vetV,c(:,:,1)',s0,v0)
% DSS=interp2(vetS(1:end-1),vetV(1:end-1),bbb(:,:,251)',s0,v0);
% % DSS=interp2(vetS(1:end-1),vetV(1:end-1),bbb(:,:,26)',s0,v0);
% % sum(sum(c(:,:,1)<0))
% % s0
% % K
% IV=blsimpv(s0, K, r, 1, option_price)
% % blsimpv(s0, K, r, 0.1, option_price)
%% plot
% figure
% % ('Renderer','zbuffer','Color',[1 1 1]);
% [XX,YY]=meshgrid(vetS,vetV); surf(XX,YY,c(:,:,1)');
% % shading interp;light;lighting gouraud;
% % colorbar
% % axis equal;
% xlabel('S');ylabel('V');zlabel('option value(S,V,1)');
%% ��������������blsimpv
% Volatility = blsimpv(Price, Strike, Rate, Time, Value, Limit, Yield, Tolerance, Class)
% K=[80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120
%     ];
% option_price_heston=[23.6677    22.9016	22.1354	21.3693	20.6031	19.837	19.0708	18.3046	17.5385	16.7723	16.0062	15.24	14.4739	13.7077	12.9416	12.1754	11.4092	10.6431	9.8769	9.1108	8.3446	8.1036	7.8627	7.6217	7.3807	7.1397	6.8988	6.6578	6.4168	6.1759	5.9349	5.6939	5.4529	5.212	4.971	4.73	4.489	4.2481	4.0071	3.7661	3.5251
%     ];
% option_price_HFB=[24.9755	24.2782	23.5941	22.9232	22.2658	21.6219	20.9915	20.3745	19.7708	19.1803	18.6026	18.0376	17.4847	16.9437	16.4141	15.8955	15.3874	14.8892	14.4006	13.9208	13.4496	13.0171	12.593	12.1777	11.7717	11.3753	10.989	10.6131	10.2477	9.8932	9.5495	9.2167	8.8947	8.5833	8.2824	7.9915	7.7103	7.4384	7.1753	6.9206	6.674
%     ];
% IV_heston=blsimpv(s0, K, r, T, option_price_heston)
% IV_HFB=blsimpv(s0, K, r, T, option_price_HFB)

% %% tau=1��
% K=[80	85	90	95	100	105	110	115	120  ];
% IV_heston=[0.3433	0.337	0.3229	0.323	0.315	0.3199	0.3177	0.3257	0.3282];
% IV_HFB_posrho=[0.3616	0.3551	0.3473	0.3456	0.3426	0.3442	0.3449	0.3489	0.3526];
% %ds=10
% K=[80	90	100	110	120  ];
% IV_heston=[0.3433	0.3229	0.315	0.3177	0.3282];
% IV_HFB_posrho=[0.3616	0.3473	0.3426	0.3449	0.3526];
% %ds=20
% IV_heston=[0.3406	0.3343	0.298	0.323	0.3218];
% IV_HFB_posrho=[0.3572	0.3518	0.3296	0.3443	0.3454];
% IV_HFB_negrho=[0.3165	0.307	0.2535	0.2935	0.293];

% %% tau=1��
% % K=[81  85	95  100	105	115 119];
% % IV_heston=[0.3419	0.337	0.323   0.298	0.3199	0.3257	0.323];
% % IV_HFB_posrho=[0.3577	0.3571	0.3424  0.3296	0.339	0.3462	0.3457];
% % IV_HFB_negrho=[0.3254    0.3188	0.2848  0.2649	0.2792	0.298   0.3028];
% % ��ȡ
% K=[85	95  100	105	115];
% IV_heston=[0.337	0.323   0.298	0.3199	0.3257];
% IV_HFB_posrho=[0.3571	0.3424  0.3296	0.339	0.3462];
% IV_HFB_negrho=[0.3188	0.2848  0.2649	0.2792	0.298];
% %% ��ͼ��
% figure
% K1=linspace(min(K),max(K));
% IV_HFB_posrho1=interp1(K,IV_HFB_posrho,K1,'cubic');
% IV_heston1=interp1(K,IV_heston,K1,'cubic');
% IV_HFB_negrho1=interp1(K,IV_HFB_negrho,K1,'cubic');
% plot ( K1 , IV_HFB_posrho1); hold on; plot ( K1 , IV_heston1,'-.'); hold on; plot ( K1 , IV_HFB_negrho1,'k--');
% % hold on; line([100 100],[0.2649 0.3296],'color','k');%vlines (100, ymin, ymax)
% legend('rho=0.05','rho=0','rho=-0.05');
% axis tight;xlabel ( 'ִ�м۸�' ) ;ylabel ( '����������' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;
% 
% %% tau=0.1��
% K=[81  85	95  100	105	115 119];
% IV_heston=[1.2387	1.2218	1.1197  1.041	1.0819	1.1021	1.0899];
% IV_HFB_posrho=[1.284	1.2657	1.1908  1.1396	1.1597	1.1686	1.1628];
% IV_HFB_negrho=[1.2014	1.1847	1.049   0.9383	1.0046	1.0428	1.0253];
% % ��ȡ
% K=[85	95  100	105	115];
% IV_heston=[1.2218	1.1197  1.041	1.0819	1.1021];
% IV_HFB_posrho=[1.2657	1.1908  1.1396	1.1597	1.1686];
% IV_HFB_negrho=[1.1847	1.049   0.9383	1.0046	1.0428];
% %% ��ͼ��
% figure
% K1=linspace(min(K),max(K));
% IV_HFB_posrho1=interp1(K,IV_HFB_posrho,K1,'cubic');
% IV_heston1=interp1(K,IV_heston,K1,'cubic');
% IV_HFB_negrho1=interp1(K,IV_HFB_negrho,K1,'cubic');
% plot ( K1 , IV_HFB_posrho1); hold on; plot ( K1 , IV_heston1,'-.'); hold on; plot ( K1 , IV_HFB_negrho1,'k--');
% legend('rho=0.05','rho=0','rho=-0.05');
% axis tight;xlabel ( 'ִ�м۸�' ) ;ylabel ( '����������' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;

% figure
% plot ( K2 , IV_HFB_posrho2); hold on; plot ( K2 , IV_heston2,'-.'); hold on; plot ( K2 , IV_HFB_negrho2,'k--');
% legend('rho>0','rho=0','rho<0');
% axis tight;xlabel ( 'ִ�м۸�' ) ;ylabel ( '����������' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;
%% ƽ����
% K1=linspace(min(K),max(K));
% IV_heston1=interp1(K,IV_heston,K1,'cubic');
% plot ( K1 , IV_heston1,'-'); 
% IV_HFB1=interp1(K,IV_HFB,K1,'cubic');
% hold on; plot ( K1 , IV_HFB1,'+-');

% % plot ( K(12:27) , IV_heston(12:27),'-'); hold on; plot ( K(12:27) , IV_HFB(12:27),'+-');
% % values_h = spcrv([[K(1) K K(end)];[IV_heston(1) IV_heston IV_heston(end)]],3);
% % values_HFB = spcrv([[K(1) K K(end)];[IV_HFB(1) IV_HFB IV_HFB(end)]],3);
% % hold on; plot(values_h(1,:),values_h(2,:), '-'); 
% % hold on; plot(values_HFB(1,:),values_HFB(2,:), '+-');
% legend('Heston','HFB');
% axis tight;xlabel ( 'K' ) ;ylabel ( '����������' ) ;set ( 0 , 'defaultfigurecolor' , 'w' ) ;



