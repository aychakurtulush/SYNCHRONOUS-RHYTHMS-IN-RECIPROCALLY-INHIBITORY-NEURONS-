function trial
%%parameters%%
g_pir=0.3;
k_syn=2;
g_syn=0.3;
V_pir=120;
V_l=-60;
V_syn=-80;
C=1;
g_l=0.1;
phi=3;
v10=-80;
h10=.1;
 
v20=-10;
h20=0.2;
 
y0=[v10;h10;v20;h20];
 
[t,y] = ode23(@InhibitoryNeurons, [0 250], y0);
 
v1=y(:,1);
h1=y(:,2);
v2=y(:,3);
h2=y(:,4);
%% figure (1)
figure(1);
clf;
plot(t,v1,'r-');
hold on;
plot(t,v2,'g-');
%% figure (2)
figure(2);
clf;
hold on;
plot(v1(100:148),h1(100:148));
V1=-120:20;
plot(V1,hinf(V1));
plot(V1,(-g_l*(V1-V_l)-g_syn*(V1-V_syn))./(g_pir*minf(V1).^3.*(V1-V_pir)));
plot(V1,-g_l*(V1-V_l)./(g_pir*minf(V1).^3.*(V1-V_pir)));
ylim([0 1]);
V01=fzero(@(V0) (-g_l*(V0-V_l)-g_syn*(V0-V_syn))./(g_pir*minf(V0).^3.*(V0-V_pir))-hinf(V0), -70);
scatter(V01,hinf(V01),'filled');
text(V01-6,hinf(V01),'V_I');
V02=fzero(@(V0) (-g_l*(V0-V_l))./(g_pir*minf(V0).^3.*(V0-V_pir))-hinf(V0), -50);
scatter(V02,hinf(V02),'filled');
text(V02-7,hinf(V02),'V_F');
xlabel('V(mV)');
ylabel('h');
%% function
    function dVdt = InhibitoryNeurons(t,y)




    v1=y(1);
    h1=y(2);
    v2=y(3);
    h2=y(4);

    dVdt = [-g_pir*minf(v1)^3*h1*(v1-V_pir)-g_l*(v1-V_l)-g_syn*sinf(v2)*(v1-V_syn); phi*(hinf(v1)-h1)/tauh(v1); -g_pir*minf(v2)^3*h2*(v2-V_pir)-g_l*(v2-V_l)-g_syn*sinf(v1)*(v2-V_syn); phi*(hinf(v2)-h2)/tauh(v2)];

    end

    function minf = minf(v)

    minf = 1./(1+ exp(-(v+65)/7.8));

    end

    function hinf = hinf(v)

    hinf = 1./(1+exp((v+81)/11));

    end

    function tauh = tauh(v)

    tauh = hinf(v)*exp((v+162.3)/17.8);

    end

    function sinf = sinf(v)
    theta_syn=-44;
    k_syn=2;

    sinf = 1/(1+exp(-(v-theta_syn)/k_syn));

    end
end