% my student ID: 4834186, thus:
a = 4; b = 3; c = 6;

A = [a - b, 0.5 - c; 0, 1];
B = [0; 1];

%% controller from prev. assignment
p1 = [-1 + 2j, -1 - 2j];
K1 = place(A, B, p1);

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
Q = [1, 0;0,-5.5];
S = inv(Q)*exp(s) * [1, s; 0, 1]*Q;
intgr1 = int(S, s, [h-tau, h])*B;
intgr2 = int(S, s, [0, h-tau])*B;

alpha1 = exp(h-tau)*(tau - h + 1);
alpha2 = exp(h-tau);

F0_poly = [exp(h),-11/2*h*exp(h),-11/2*exp(h)*(h-1);0,exp(h),1;0,0,0];
F1_poly = [0, 0, -11/2; zeros(2,3)];
F2_poly = [zeros(1,3);0,0,-1;zeros(1,3)];
G0_poly = [11/2;-1;1];
G1_poly = [11/2;0;0];
G2_poly = [0;1;0];

Q = 1E-10*eye(3);

taus = 0:0.05:0.7;
hs = 0.01:0.1:0.7;

rho4 = zeros(length(hs), length(taus));
figure(6)
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
            
            a1s4 = a1s(1,1:j);
            a2s4 = a2s(1,1:j);
            %assess stability
            P = sdpvar(3);
            gamma = 0;
            F0_poly_i = double(subs(F0_poly, h,hs(i)));
            H_F1 = F0_poly_i + min(a1s4) * F1_poly + max(a2s4) * F2_poly;
            H_F2 = F0_poly_i + max(a1s4) * F1_poly + max(a2s4) * F2_poly;
            %H_F3 = F0_poly_i + min(a1s4) * F1_poly + min(a2s4) * F2_poly;
            H_F4 = F0_poly_i + max(a1s4) * F1_poly + min(a2s4) * F2_poly;
            H_G1 = G0_poly + min(a1s4) * G1_poly + max(a2s4) * G2_poly;
            H_G2 = G0_poly + max(a1s4) * G1_poly + max(a2s4) * G2_poly;
            %H_G3 = G0_poly + min(a1s4) * G1_poly + min(a2s4) * G2_poly;
            H_G4 = G0_poly + max(a1s4) * G1_poly + min(a2s4) * G2_poly;
            cons1 = P >= Q;
            K = [K1,0];
            cons2 = (H_F1-H_G1*K)' * P * (H_F1-H_G1*K) - P <= -1*gamma*P;
            cons3 = (H_F2-H_G2*K)' * P * (H_F2-H_G2*K) - P <= -1*gamma*P;
            %cons4 = (H_F3-H_G3*K)' * P * (H_F3-H_G3*K) - P <= -1*gamma*P;
            cons5 = (H_F4-H_G4*K)' * P * (H_F4-H_G4*K) - P <= -1*gamma*P;

            cons = [cons1,cons2,cons3,cons5];
            options = sdpsettings('verbose',0,'solver', "sedumi");
            obj = 0;
            result = optimize(cons, obj, options);
            result.info
            if strcmp(result.info, 'Infeasible problem (SeDuMi-1.3)')
                rho4(i,j) = 1;
            else
                rho4(i,j) = 0;
            end
        end
    end
end
spy(rho4,"red")
ylabel("sampling time [h]")
xlabel("delay time [\tau]")
hold on
spy(ones(size(rho4))-rho4, "green")