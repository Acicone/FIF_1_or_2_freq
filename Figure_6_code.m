% Test artifical example

% We go back to the standard stopping criterion, but we shift the
% eigenvalues of fft(a) to ensure that we have a 0 correpsonding to the 1
% freq (position automatically identified by the algorithm in the fft)

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
    df=sqrt(2)/140;
    da=0.05;
    dphi=1;
    x1=cos(2*pi*t);
    
    Samples_per_Period=N/T;
    
    %% FIF_v2 using the prefixed_double_filter
    
    xi=2*207/200;
    
    Nf=99;
    
    c_1=zeros(Nf,length(-2:da:2),length(0:dphi:2*pi));
    err=zeros(Nf,length(-2:da:2),length(0:dphi:2*pi));
    ii=0;
    for nf = 1:Nf% 50
        f=df*nf;
        ii=ii+1
        jj=0;
        for exp_a = -2:da:2%-0.5
            jj=jj+1;
            kk=0;
            for phi = 0:dphi:2*pi
                kk=kk+1;
                
                x2=10^exp_a*cos(2*pi*f*t+phi);
                
                x=x1+x2;
                
                opts=Settings_IF_v1('IF.Xi',xi,'IF.alpha','ave','IF.NIMFs',1,'IF.MaxInner',10000000,'verbose',0,'IF.delta',10^-20,'plots',0);
                
                [IMF,posF(ii,jj,kk)] = FIF_v2_6_ver_1_or_2_freq(x,opts);
                
                c_1(ii,jj,kk)=norm(IMF(1,:)-x1,2)/norm(x2,2);
                
                err(ii,jj,kk)=norm(IMF(1,:)-x1,2);
                
            end
        end
    end
    %%
    
    
    save(['Artifical_example_StopCrit_18_ver2'])
    %%
    for iii=4
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df*(1:Nf));
    surf(X,Y,c_1(1:end,:,iii),'edgecolor','none');
    %title(['Case T = ' num2str(nn*60) ' delta = 10^-' num2str(abs(log10(opts.IF.delta)))])
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %hold on
    %view(0,90)
%     plot3(-2:da:2,1./10.^(-2:da:2),0*ones(1,101),'r','linewidth',3)
%     axis([-2,2,df,Df,-Inf,Inf])
%     plot3(-2:da:2,1./sqrt(10.^(-2:da:2)),1*ones(1,101),'k','linewidth',3)
     xlabel('$\log_{10}(a)$','Interpreter','latex');
    ylabel('$f$','Interpreter','latex');
    set(gca,'fontsize', 30);
    %view(0,90)
    %pause(1)
    shading interp
    end
    saveas(fig,'Artifical_example_StopCrit_18_ver2', 'epsc')
    saveas(fig,'Artifical_example_StopCrit_18_ver2', 'png')
    saveas(fig,'Artifical_example_StopCrit_18_ver2', 'fig')
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df*(1:Nf));
    surf(X,Y,c_1(1:end,:,4),'edgecolor','none');
    %title(['Case T = ' num2str(nn*60) ' delta = 10^-' num2str(abs(log10(opts.IF.delta)))])
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    xlabel('$\log_{10}(a)$','Interpreter','latex');
   ylabel('$f$','Interpreter','latex');
    set(gca,'fontsize', 30);
    %view(0,90)
    %colormap(fig,parula(200000));
    hold on
    plot3(-2:da:2,1./10.^(-2:da:2),2*ones(1,81),'r','linewidth',3)
    axis([-2,2,df,df*Nf,-Inf,Inf])
    plot3(-2:da:2,1./sqrt(10.^(-2:da:2)),2*ones(1,81),'k','linewidth',3)
    %axis([-2,2,df,Nf*df,0,1])
    view(0,90)
    shading interp
    saveas(fig,'Artifical_example_StopCrit_18_ver2_top', 'epsc')
    saveas(fig,'Artifical_example_StopCrit_18_ver2_top', 'png')
    saveas(fig,'Artifical_example_StopCrit_18_ver2_top', 'fig')
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df*(1:Nf));
    surf(X,Y,mean(err(1:end,:,:),3),'edgecolor','none');
    %title(['Case T = ' num2str(nn*60) ' delta = 10^-' num2str(abs(log10(opts.IF.delta)))])
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    set(gca,'fontsize', 30);
    saveas(fig,'Artifical_example_StopCrit_18_ver2_New_Err', 'epsc')
    saveas(fig,'Artifical_example_StopCrit_18_ver2_New_Err', 'png')
    saveas(fig,'Artifical_example_StopCrit_18_ver2_New_Err', 'fig')
    
end