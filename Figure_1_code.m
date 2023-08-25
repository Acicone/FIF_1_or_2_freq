% Test artifical example

% We go back to the standard stopping criterion, but we shift the
% eigenvalues of fft(a) to ensure that we have a 0 correpsonding to the 1
% freq (position 200 in the fft)

clear all
close all

for nn=10%:10
    
    T=20*nn; % Period
    N=80000; % Number of sample points
    dt=T/N;
    t=dt:dt:T;
    
    if N/T<4 % checking for the Nyquest rate
        disp('Not enough sample points per period')
        break
    end
    df=0.01;
    da=0.04;
    dphi=1;
    x1=cos(2*pi*t);
    
    Samples_per_Period=N/T;
    
    %% FIF_v2 using the prefixed_double_filter
    
    xi=2*207/200;
    
    c_1=zeros(length(df:df:0.999),length(-2:da:2),length(0:dphi:2*pi));
    err=zeros(length(df:df:0.999),length(-2:da:2),length(0:dphi:2*pi));
    ii=0;
    for f = df:df:1
        ii=ii+1
        jj=0;
        for exp_a = -2:da:2
            jj=jj+1;
            kk=0;
            for phi = 0:dphi:2*pi
                kk=kk+1;
                
                x2=10^exp_a*cos(2*pi*f*t+phi);
                
                x=x1+x2;
                
                opts=Settings_IF_v1('IF.Xi',xi,'IF.alpha',1,'IF.NIMFs',1,'IF.MaxInner',10000000,'verbose',0,'IF.delta',10^-35,'plots',0);
                
                [IMF,logM] = FIF_v2_5(x,opts,Samples_per_Period*xi);
                
                c_1(ii,jj,kk)=norm(IMF(1,:)-x1,2)/norm(x2,2);
                
                err(ii,jj,kk)=norm(IMF(1,:)-x1,2);
                
            end
        end
    end
    
    save(['Artifical_example_StopCrit_5_ver2'])
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df:df:1);
    surf(X,Y,mean(c_1(:,:,:),3),'edgecolor','none');
   xlabel('$\log_{10}(a)$','Interpreter','latex');
   ylabel('$f$','Interpreter','latex');
   
   shading interp
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gca,'fontsize', 30);
%     saveas(fig,'Artifical_example_StopCrit_5_ver2_zoomed_in', 'png')
%     saveas(fig,'Artifical_example_StopCrit_5_ver2_zoomed_in', 'epsc')
% 
%
%     
     axis([-2,2,df,1,0,1])
    
    saveas(fig,'Artifical_example_StopCrit_5_ver2', 'png')
    saveas(fig,'Artifical_example_StopCrit_5_ver2', 'fig')
    saveas(fig,'Artifical_example_StopCrit_5_ver2', 'epsc')
    %%
    set(gca,'zscale','log')
    saveas(fig,'Artifical_example_StopCrit_5_ver2_log_scale', 'epsc')
    
    
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df:df:1);
    surf(X,Y,mean(err(:,:,:),3),'edgecolor','none');
   xlabel('$\log_{10}(a)$','Interpreter','latex');
   ylabel('$f$','Interpreter','latex');
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gca,'fontsize', 30);
    saveas(fig,'Artifical_example_StopCrit_5_ver2_New_Err', 'png')
    saveas(fig,'Artifical_example_StopCrit_5_ver2_New_Err', 'epsc')

%%
    
    
end