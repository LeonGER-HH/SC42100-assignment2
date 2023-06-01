% my student ID: 4834186, thus:
a = 4; b = 3; c = 6;

A = [a - b, 0.5 - c; 0, 1];
B = [0; 1];

%% controller from prev. assignment
p1 = [-1 + 2j, -1 - 2j];
p2 = [-1, -2];
K1 = place(A, B, p1);
K2 = place(A, B, p2);

%% packet loss
syms h s
F_h = expm(A * h);
G_h = int(expm(A*s)*B, s, [0 h]);

%x^e = x1;x2;u_k-1
%to hold approach
Fh_0 = [F_h-G_h*K1, [0;0]; -1*K1, 0];
Fh_1 = [F_h, G_h; [0,0], 1];
%to zero approach
Fz_0 = F_h-G_h*K1;
Fz_1 = F_h;

%% Q1 stability
hs = 0.1:0.005:0.6;
Fh_seq = Fh_1 * Fh_0 * Fh_0;
Fz_seq = Fz_1 * Fz_0 * Fz_0;

rho_to_hold = zeros(1,length(hs));
rho_to_zero = zeros(1,length(hs));
for i=1:length(hs)
    Fh_seq_n = double(subs(Fh_seq, hs(i)));
    Fz_seq_n = double(subs(Fz_seq, hs(i)));
    rho_to_hold(1,i) = max(abs(eig(Fh_seq_n)));
    rho_to_zero(1,i) = max(abs(eig(Fz_seq_n)));
end
figure(1)
plot(hs, rho_to_hold, hs, rho_to_zero, hs, ones(1,length(hs)), "--")
grid("on")
ylim([0,2])
xlabel("sampling time h[s]")
ylabel("\rho[-]")
legend("to-hold", "to-zero")
saveas(gcf,'figures4report/q1-stability.png')
%% Q2 non-deterministic dropout
hs = 0.1:0.005:0.6;
%to hold
Fh_seq_0 = Fh_0*Fh_0*Fh_0;
Fh_seq_1 = Fh_1*Fh_0*Fh_0;
%to zero
Fz_seq_0 = Fz_0*Fz_0*Fz_0;
Fz_seq_1 = Fz_1*Fz_0*Fz_0;

Q = 0.001*eye(3);
rhoh = ones(1, length(hs));
for i=1:length(hs) % to hold solving
    P = sdpvar(3, 3);
    F_seq_0_i = double(subs(Fh_seq_0, hs(i)));
    F_seq_1_i = double(subs(Fh_seq_1, hs(i)));
    cons = [P >= 0, F_seq_1_i'*P*F_seq_1_i-P <= -Q, F_seq_0_i'*P*F_seq_0_i-P <= -Q];
    obj = 0;
    options = sdpsettings('verbose',0,'solver', "sedumi");
    result = optimize(cons, obj, options);
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        rhoh(i) = 0;
    end
end

Q = 0.001*eye(2);
rhoz = ones(1, length(hs));
for i=1:length(hs) % to zero solving
    P = sdpvar(2, 2);
    F_seq_0_i = double(subs(Fz_seq_0, hs(i)));
    F_seq_1_i = double(subs(Fz_seq_1, hs(i)));
    cons = [P >= 0, F_seq_1_i'*P*F_seq_1_i-P <= -Q, F_seq_0_i'*P*F_seq_0_i-P <= -Q];
    obj = 0;
    options = sdpsettings('verbose',0,'solver', "sedumi");
    result = optimize(cons, obj, options);
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        rhoz(i) = 0;
    end
end

%% Question 2 - as a function of h
figure(2)
plot(hs,rho_to_hold, hs, rhoh)
ylim([0,2])
xlabel("sampling time h[s]")
ylabel("\rho[-] / binary")
legend("deterministic", "non-deterministic")
saveas(gcf,'figures4report/q2-LMI-stability-to_hold.png')

figure(3)
plot(hs,rho_to_zero, hs, rhoz)
ylim([0,2])
xlabel("sampling time h[s]")
ylabel("\rho[-] / binary")
legend("deterministic", "non-deterministic")
saveas(gcf,'figures4report/q2-LMI-stability-to_zero.png')

%% Question 3 - probabilistic dropout
hs = 0.1:0.05:0.6;

p = [0, 0.99, 0, 0, 0.01, 0;
    0, 0, 0.49, 0, 0, 0.51;
    0.99, 0, 0, 0.01, 0, 0;
    0, 0.99, 0, 0, 0.01, 0;
    0, 0, 0.49, 0, 0, 0.51;
    0.99, 0, 0, 0.01, 0, 0];
N = length(P);
%to zero approach
rho_z = ones(1,length(hs));
Q = 0.001*eye(2);
rho_cnt = 1;
for h=hs
    Ps = {sdpvar(2, 2),sdpvar(2, 2),sdpvar(2, 2),sdpvar(2, 2),sdpvar(2, 2),sdpvar(2, 2)};
    Fz_0_i = double(subs(Fz_0, h));
    Fz_1_i = double(subs(Fz_1, h));
    Fs = {Fz_0_i, Fz_0_i, Fz_0_i, Fz_1_i, Fz_1_i, Fz_1_i};
    cons = [];
    for i=1:N
        cons = [cons, Ps{i} >= Q];
        summation = 0;
        for j=1:N
            summation = summation + p(i,j)*Ps{j};
        end
        cons = [cons, -1*Ps{i} + Fs{i}' * (summation) * Fs{i} <= -1*Q];
    end
    obj = 0;
    options = sdpsettings('verbose',0,'solver', "sedumi");
    result = optimize(cons, obj, options);
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        rho_z(rho_cnt) = 0;
    end
    rho_cnt = rho_cnt + 1;
end

%to hold approach
rho_h = ones(1,length(hs));
Q = 0.001*eye(3);
rho_cnt = 1;
for h=hs
    Ps = {sdpvar(3, 3),sdpvar(3, 3),sdpvar(3, 3),sdpvar(3, 3),sdpvar(3, 3),sdpvar(3, 3)};
    Fh_0_i = double(subs(Fh_0, h));
    Fh_1_i = double(subs(Fh_1, h));
    Fs = {Fh_0_i, Fh_0_i, Fh_0_i, Fh_1_i, Fh_1_i, Fh_1_i};
    cons = [];
    for i=1:N
        cons = [cons, Ps{i} >= Q];
        summation = 0;
        for j=1:N
            summation = summation + p(i,j)*Ps{j};
        end
        cons = [cons, -1*Ps{i} + Fs{i}' * (summation) * Fs{i} <= -1*Q];
    end
    obj = 0;
    options = sdpsettings('verbose',0,'solver', "sedumi");
    result = optimize(cons, obj, options);
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        rho_h(rho_cnt) = 0;
    end
    rho_cnt = rho_cnt + 1;
end

figure(4)
plot(hs, rho_h, hs, rho_z)
legend("to-hold", "to-zero")
xlabel("sampling time h[s]")
ylabel("instability")
saveas(gcf,'figures4report/q3-LMI-stability.png')

%% q4
% matrices depending on tau
syms s tau h
Q = [1, 0;0,-5.5];
S = inv(Q)*exp(s) * [1, s; 0, 1]*Q;
intgr1 = int(S, s, [h-tau, h])*B;
intgr2 = int(S, s, [0, h-tau])*B;

alpha1 = exp(h-tau)*(tau - h + 1);
alpha2 = exp(h-tau);

taus = 0:0.05:0.7;
hs = 0.01:0.05:0.7;
figure(5)
for i=1:length(hs)
    alpha1_i = subs(alpha1, h, hs(i));
    alpha2_i = subs(alpha2, h, hs(i));
    a1s = zeros(1,length(taus));
    a2s = zeros(1,length(taus));
    taus = hs(1:i-1);
    for j=1:length(taus)
        if taus(j) < hs(i)
            a1s(1,j) = double(subs(alpha1_i, tau, taus(j)));
            a2s(1,j) = double(subs(alpha2_i, tau, taus(j)));
        end
    end
    plot(a1s, a2s, "blue",'LineWidth',2)
    xlim([0.65,1.025])
    ylim([1,2])
    xlabel("\alpha_1(\tau)")
    ylabel("\alpha_2(\tau)")
    hold on
%     pgon = polyshape([max(a1s), max(a1s), min(a1s), min(a1s)], [min(a2s), max(a2s), max(a2s), min(a2s)]);
    pgon = polyshape([max(a1s), max(a1s), min(a1s)], [min(a2s), max(a2s), max(a2s)]);
    plot(pgon,'FaceColor','green','FaceAlpha',0.01)
    pause(0.05)
end
saveas(gcf,'figures4report/q4-polytope-overaproximation-triangle.png')
%% q4
% matrices depending on tau
syms s tau h
EV = [1, 0;0,-5.5];
S = inv(EV)*exp(s) * [1, s; 0, 1]*EV;
intgr1 = int(S, s, [h-tau, h])*B;
intgr2 = int(S, s, [0, h-tau])*B;

%%

alpha1 = exp(h-tau)*(tau - h + 1);
alpha2 = exp(h-tau);

F0_poly = [exp(h),-11/2*h*exp(h),-11/2*exp(h)*(h-1);0,exp(h),exp(h);0,0,0];
F1_poly = [0, 0, -11/2; zeros(2,3)];
F2_poly = [zeros(1,3);0,0,-1;zeros(1,3)];
G0_poly = [-11/2;-1;1];
G1_poly = [11/2;0;0];
G2_poly = [0;1;0];

Q = 1E-3*eye(3);

taus = 0:0.05:0.7;
hs = 0.01:0.01:0.2;

rho4 = ones(1,length(hs));
for i=1:length(hs)
    alpha1_i = subs(alpha1, h, hs(i));
    alpha2_i = subs(alpha2, h, hs(i));
    a1s = zeros(1,length(taus));
    a2s = zeros(1,length(taus));
    for j=1:length(taus)
        if taus(j) < hs(i)
            a1s(1,j) = double(subs(alpha1_i, tau, taus(j)));
            a2s(1,j) = double(subs(alpha2_i, tau, taus(j)));
        end
    end
    %assess stability
    P = sdpvar(3,3,'symmetric');
    gamma = 0.0001;
    a1_min = double(subs(alpha1_i, tau, 0));
    a1_max = double(subs(alpha1_i, tau, hs(i)));
    a2_min = double(subs(alpha2_i, tau, hs(i)));
    a2_max = double(subs(alpha2_i, tau, 0));


    F0_poly_i = double(subs(F0_poly, h,hs(i)));
    H_F1 = F0_poly_i + a1_min * F1_poly + a2_max * F2_poly;
    H_F2 = F0_poly_i + a1_max * F1_poly + a2_max * F2_poly;
    H_F3 = F0_poly_i + a1_min * F1_poly + a2_min * F2_poly;
    H_F4 = F0_poly_i + a1_max * F1_poly + a2_min * F2_poly;
    H_G1 = G0_poly + a1_min * G1_poly + a2_max * G2_poly;
    H_G2 = G0_poly + a1_max * G1_poly + a2_max * G2_poly;
    H_G3 = G0_poly + a1_min * G1_poly + a2_min * G2_poly;
    H_G4 = G0_poly + a1_max * G1_poly + a2_min * G2_poly;
    cons0 = P >= Q;
    K = [K1,0];
    cons1 = (H_F1-H_G1*K)' * P * (H_F1-H_G1*K) - P <= -1*gamma*P;
    cons2 = (H_F2-H_G2*K)' * P * (H_F2-H_G2*K) - P <= -1*gamma*P;
    cons3 = (H_F4-H_G4*K)' * P * (H_F4-H_G4*K) - P <= -1*gamma*P;
    cons4 = (H_F3-H_G3*K)' * P * (H_F3-H_G3*K) - P <= -1*gamma*P;

    cons = [cons0, cons1,cons2,cons3];
    options = sdpsettings('verbose',0,'solver', "sedumi");
    obj = 0;
    result = optimize(cons, obj, options);
    result.info
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        disp("STABLE")
        rho4(1,i) = 0;
    else
        result.info
        rho4(1,i) = 1;
    end
end

scatter(hs,rho4, "ro")

%% Q4 better approximation
alpha1 = exp(h-tau)*(tau - h + 1);
alpha2 = exp(h-tau);

taus = 0:0.05:0.7;
hs = 0.1:0.1:0.7;
for i=1:length(hs)
    alpha1_i = subs(alpha1, h, hs(i));
    alpha2_i = subs(alpha2, h, hs(i));
    a1s = zeros(1,length(taus));
    a2s = zeros(1,length(taus));
    taus = hs(1:i-1);
    for j=1:length(taus)
        if taus(j) < hs(i)
            a1s(1,j) = double(subs(alpha1_i, tau, taus(j)));
            a2s(1,j) = double(subs(alpha2_i, tau, taus(j)));
        end
    end
    plot(a1s, a2s, "blue",'LineWidth',2)
    xlim([0.65,1.025])
    ylim([1,2])
    xlabel("\alpha_1(\tau)")
    ylabel("\alpha_2(\tau)")
    hold on
    a1_min = min(a1s); %double(subs(alpha1_i, tau, 0));
    a1_max = max(a1s); %double(subs(alpha1_i, tau, hs(i)));
    a2_min = min(a2s); %double(subs(alpha2_i, tau, hs(i)));
    a2_max = max(a2s); %double(subs(alpha2_i, tau, 0));
%     pgon = polyshape([max(a1s), max(a1s), min(a1s), min(a1s)], [min(a2s), max(a2s), max(a2s), min(a2s)]);
    pgon = polyshape([a1_min, a1_min+(a1_max-a1_min)/5, a1_max, a1_max],[a2_max, a2_max, a2_min+(a2_max-a2_min)/5, a2_min]);
    plot(pgon,'FaceColor','green','FaceAlpha',0.01)
    pause(0.05)
end
%% stability analysis better approximation

alpha1 = exp(h-tau)*(tau - h + 1);
alpha2 = exp(h-tau);

F0_poly = [exp(h),-11/2*h*exp(h),-11/2*exp(h)*(h-1);0,exp(h),exp(h);0,0,0];
F1_poly = [0, 0, -11/2; zeros(2,3)];
F2_poly = [zeros(1,3);0,0,-1;zeros(1,3)];
G0_poly = [-11/2;-1;1];
G1_poly = [11/2;0;0];
G2_poly = [0;1;0];

Q = 1E-3*eye(3);

taus = 0:0.05:0.7;
hs = 0.01:0.01:0.2;

rho4 = ones(1,length(hs));
for i=1:length(hs)
    alpha1_i = subs(alpha1, h, hs(i));
    alpha2_i = subs(alpha2, h, hs(i));
    a1s = zeros(1,length(taus));
    a2s = zeros(1,length(taus));
    for j=1:length(taus)
        if taus(j) < hs(i)
            a1s(1,j) = double(subs(alpha1_i, tau, taus(j)));
            a2s(1,j) = double(subs(alpha2_i, tau, taus(j)));
        end
    end
    %assess stability
    P = sdpvar(3,3,'symmetric');
    gamma = 0.0001;
    a1_min = double(subs(alpha1_i, tau, 0));
    a1_max = double(subs(alpha1_i, tau, hs(i)));
    a2_min = double(subs(alpha2_i, tau, hs(i)));
    a2_max = double(subs(alpha2_i, tau, 0));


    F0_poly_i = double(subs(F0_poly, h,hs(i)));
%     polyshape([a1_min, a1_min+(a1_max-a1_min)/5, a1_max, a1_max],[a2_max, a2_max, a2_min+(a2_max-a2_min)/5, a2_min]);
    H_F1 = F0_poly_i + a1_min * F1_poly + a2_max * F2_poly;
    H_F2 = F0_poly_i + (a1_min+(a1_max-a1_min)/5) * F1_poly + a2_max * F2_poly;
    H_F3 = F0_poly_i + a1_max * F1_poly + (a2_min+(a2_max-a2_min)/5) * F2_poly;
    H_F4 = F0_poly_i + a1_max * F1_poly + a2_min * F2_poly;
    H_G1 = G0_poly + a1_min * G1_poly + a2_max * G2_poly;
    H_G2 = G0_poly + (a1_min+(a1_max-a1_min)/5) * G1_poly + a2_max * G2_poly;
    H_G3 = G0_poly + a1_max * G1_poly + (a2_min+(a2_max-a2_min)/5) * G2_poly;
    H_G4 = G0_poly + a1_max * G1_poly + a2_min * G2_poly;
    cons0 = P >= Q;
    K = [K1,0];
    cons1 = (H_F1-H_G1*K)' * P * (H_F1-H_G1*K) - P <= -1*gamma*P;
    cons2 = (H_F2-H_G2*K)' * P * (H_F2-H_G2*K) - P <= -1*gamma*P;
    cons3 = (H_F4-H_G4*K)' * P * (H_F4-H_G4*K) - P <= -1*gamma*P;
    cons4 = (H_F3-H_G3*K)' * P * (H_F3-H_G3*K) - P <= -1*gamma*P;

    cons = [cons0, cons1,cons2,cons3,cons4];
    options = sdpsettings('verbose',0,'solver', "sedumi");
    obj = 0;
    result = optimize(cons, obj, options);
    result.info
    if strcmp(result.info, 'Successfully solved (SeDuMi-1.3)')
        disp("STABLE")
        rho4(1,i) = 0;
    else
        result.info
        rho4(1,i) = 1;
    end
end

scatter(hs,rho4, "ro")

%% Q5
sys_ol = ss(A,B,eye(2), [0;0]);
sys_cl = ss(A-B*K1,[0;0],eye(2), [0;0]);
t = 0:0.05:10;
lsim(sys_cl, zeros(1,length(t)),t, [1;0])

% [t,x] = ode45(@cont_sys, tspan, x0);
[t,x] = ode45(@(t,x) (A-B*K1)*x, tspan, x0);
plot(t,x,'-o')
%% Q5 - 1
P = sdpvar(2,2);
Q = 0.00001 * eye(2);
obj = 0;
options = sdpsettings('verbose',0,'solver', "sedumi");
cons = [P >= Q, (A-B*K1)'*P+P*(A-B*K1) <=0];
result = optimize(cons, obj, options);
result.info
P = value(P);
syms sigmaETC
Qn = -1*((A-B*K1)'*P+P*(A-B*K1));
ETC = [A'*P+P*A+sigmaETC*Qn, -P*B*K1; -(B*K1)'*P, zeros(2,2)];


%% Q5-2
figure(6)
time_interval = 0.1;
ts = 0:time_interval:50;
x0 = [1;0];
% x0 = [10;-5];
% x0 = [0;-1];
sigmas = 0.1;%[0.05,0.1:0.1:0.9];
for sigma=sigmas
    xi = [x0];
    xk = [1;0];
    sk = [];
    k = 1;
    u = 0;
    u_rec = [];
    samples = 0;
    for i=2:length(ts)
        tspan = ts(i-1):0.05:ts(i);
    %     if (xi(1,end) - xk(1,end))^2 > rho*xi(1,end)^2 % trigger condition
        if [xi(:,end)', xk(:,end)']*subs(ETC,sigmaETC,sigma)*[xi(:,end); xk(:,end)] >= 0 % GES trigger condition
            u = -K1*xi(:,end);
            xk = [xk, xi(:,end)];
            samples = samples + 1;
            sk = [sk, ts(i-1)];
        end
        xi0 = xi(:,end);
    %     [t,xi_new] = ode45(@(t,x) (A-B*K1)*x, tspan, xi0);
        [t,xi_new] = ode45(@(t,x) A*x+B*u(1,1), tspan, xi0);
        xi = [xi xi_new(2:end,:)'];
        u_rec = [u_rec,u];
    end
    tiledlayout(2,1);
    nexttile
    time_axis = ts(1):(ts(end)-ts(1))/length(xi):ts(end);
    plot(time_axis(1,1:end-1),xi(1,:))
    hold on
    plot(time_axis(1,1:end-1),xi(2,:))
    legend("\xi_1(t)", "\xi_2(t)")
    xlabel("time t[s]")
    ylabel("x_i[-]")
    nexttile
    time_axis = ts(1):(ts(end)-ts(1))/length(u_rec):ts(end);
    stairs(time_axis(1,1:end-1),u_rec)
    hold on
    plot(sk, zeros(1,length(sk)),"r|")
    legend("u","s_k")
    ylabel("u[-]")
    xlabel("time t[s]")
    % nexttile
    % time_axis = ts(1):(ts(end)-ts(1))/length(xk):ts(end);
    % stairs(time_axis(1,1:end-1),xk(2,:))
    % legend("x_{k1}")
    
    %calc average inter sample time
    dist = 0;
    for i=1:length(sk)-1
        dist = dist + sk(i+1)-sk(i);
    end
    h_avg=dist/samples;
    disp([sigma, "#samples=", samples, " | h_avg=", h_avg])
end
saveas(gcf,'figures4report/q5-ETC-simulation1.png')

%% Q5 - 3 checking whether time periodic sampling retains GES
% verify by checking spectral radius
FHH=double(subs(F_h,h,0.48));
GHH = double(subs(G_h,h,0.48));
max(abs(eig(FHH-GHH*K1)))
%% Q5-4 periodic event triggering

