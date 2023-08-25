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
    df=sqrt(2)/140;
    da=0.04;
    dphi=1;
    x1=cos(2*pi*t);
    
    Samples_per_Period=N/T;
    
    %% FIF_v2 using the prefixed_double_filter
    
    xi=2*207/200;
    Nf=99;
    
    c_1=zeros(Nf,length(-2:da:2),length(0:dphi:2*pi));
    err=zeros(Nf,length(-2:da:2),length(0:dphi:2*pi));
    ii=0;
    for nf = 1:Nf
        f=df*nf;
        ii=ii+1
        jj=0;
        for exp_a = -2:da:2
            jj=jj+1;
            kk=0;
            for phi = 0:dphi:2*pi
                kk=kk+1;
                
                x2=10^exp_a*cos(2*pi*f*t+phi);
                
                x=x1+x2;
                
                opts=Settings_IF_v1('IF.Xi',xi,'IF.alpha',1,'IF.NIMFs',1,'IF.MaxInner',10000000,'verbose',0,'IF.delta',10^-22,'plots',0);
                
                [IMF,logM] = FIF_v2_5(x,opts,Samples_per_Period*xi);
                
                c_1(ii,jj,kk)=norm(IMF(1,:)-x1,2)/norm(x2,2);
                
                err(ii,jj,kk)=norm(IMF(1,:)-x1,2);
                
            end
        end
    end
    
    save(['Artifical_example_StopCrit_15'])
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df*(1:Nf-1));
    surf(X,Y,mean(c_1(1:end-1,:,:),3),'edgecolor','none');
    set(gca,'fontsize', 30);
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    xlabel('$\log_{10}(a)$','Interpreter','latex');
   ylabel('$f$','Interpreter','latex');
    saveas(fig,'Artifical_example_StopCrit_15_zoomed_in', 'png')
    saveas(fig,'Artifical_example_StopCrit_15_zoomed_in', 'epsc')
    
    %%
    
    axis([-2,2,df,df*(Nf-1),0,1])
    saveas(fig,'Artifical_example_StopCrit_15', 'png')
    saveas(fig,'Artifical_example_StopCrit_15', 'fig')
    saveas(fig,'Artifical_example_StopCrit_15', 'epsc')
    
    %%
    fig=figure;
    [X,Y] = meshgrid(-2:da:2,df*(1:Nf-1));
    surf(X,Y,mean(err(1:end-1,:,:),3),'edgecolor','none');
    xlabel('$\log_{10}(a)$','Interpreter','latex');
   ylabel('$f$','Interpreter','latex');
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gca,'fontsize', 30);
    saveas(fig,'Artifical_example_StopCrit_15_New_Err', 'png')
    saveas(fig,'Artifical_example_StopCrit_15_New_Err', 'epsc')
    
    %%
    
    
end