clear all
close all
% 中间点初值  初值赋值方式：[theta1'*m;theta2'*m;...;t']
% q0 = [pi/5,-3*pi/10;
%       2*pi/5,-1*pi/10;
%       3*pi/5,1*pi/10;
%       4*pi/5,3*pi/10];
piece_num = 5;
theta_num = 2;
T0 = 1;
t0 = realT2virtualT(T0);
pe0 = [pi,pi/2];
q0 = [pi/6;pi/3;pi/2;2*pi/3;
          -1*pi/3;0;pi/6;pi/3;
          pe0';
          t0;];
  
% 优化函数
fun = @MINCO;
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'Display','iter','MaxIterations',500,'MaxFunctionEvaluations',500 );
[q_ini, FVAL, EXITFLAG, OUTPUT] = fminunc(fun,q0,options);
OUTPUT.iterations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q = zeros(piece_num-1,theta_num);
    q(:,1) = q_ini(1:piece_num-1);
    q(:,2) = q_ini(piece_num:2*piece_num-2);

% 作图
% 杆长度
l = [1,0.5];
vir_t = q_ini(end);
% t=[1,0.8,0.5,0.7,1];
t = virtualT2realT(vir_t);
%初始状态
p0 = [0,-pi/2];
% pe = [pi,pi/2];
pe = q_ini(2*piece_num-1:2*piece_num)';
%更新M b
    M = zeros(8*piece_num,8*piece_num);
    b = zeros(8*piece_num,2);
    M(1:4,1:8) = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 2 0 0 0 0 0 ;0 0 0 6 0 0 0 0];
    b(1:4,:) = [p0;0,0;0,0;0,0];
    M_part = [1,t,t^2,t^3,t^4,t^5,t^6,t^7;
              1,t,t^2,t^3,t^4,t^5,t^6,t^7;
              0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6;
              0,0,2,6*t,12*t^2,20*t^3,30*t^4,42*t^5;
              0,0,0,6,24*t,60*t^2,120*t^3,210*t^4;
              0,0,0,0,24,120*t,360*t^2,840*t^3;
              0,0,0,0,0,120,720*t,2520*t^2;
              0,0,0,0,0,0,720,5040*t];
    for i=1:piece_num-1
        % E
        M(5+8*(i-1):12+8*(i-1),1+8*(i-1):8+8*(i-1)) = M_part;
        % F
        M(5+8*(i-1):12+8*(i-1),9+8*(i-1):16+8*(i-1)) = [0,0,0,0,0,0,0,0;
                                                      -1,0,0,0,0,0,0,0;
                                                      0,-1,0,0,0,0,0,0;
                                                      0,0,-2,0,0,0,0,0;
                                                      0,0,0,-6,0,0,0,0;
                                                      0,0,0,0,-24,0,0,0;
                                                      0,0,0,0,0,-120,0,0;
                                                      0,0,0,0,0,0,-720,0];
        b(5+8*(i-1):12+8*(i-1),:) = [q(i,:);0,0;0,0;0,0;0,0;0,0;0,0;0,0];
    end
    M(end-3:end,end-7:end) = [M_part(2,:);M_part(3,:);M_part(4,:);M_part(5,:)];
    b(end-3:end,:) = [pe;0,0;0,0;0,0];
c = M\b;


figure();
hold on;

%原点位置
x0 = [0,0,0,1]';

theta1 = [];
theta2 = [];
dtheta1 = [];
dtheta2 = [];
ddtheta1 = [];
dddtheta1 = [];
ddddtheta1 = [];

x3 = [];
x2 = [];
v = [];
plotT = [];
step = 0.01;
for i = 1:piece_num
    for T = 0:step:t
%         if abs(T-t(i))<=0.0001&&i~=piece_num
%             size(theta1);
%             continue;
%         end
        if i==1
            sumt = 0;
        else
            sumt = (i-1)*t;
        end
        plotT = [plotT,sumt+T];

        beta = [1,T,T^2,T^3,T^4,T^5,T^6,T^7]';
        dbeta = [0,1,2*T,3*T^2,4*T^3,5*T^4,6*T^5,7*T^6]';
        ddbeta = [0,0,2,6*T,12*T^2,20*T^3,30*T^4,42*T^5]';
        dddbeta = [0,0,0,6,24*T,60*T^2,120*T^3,210*T^4]';
        ddddbeta = [0,0,0,0,24,120*T,360*T^2,840*T^3]';
        
        theta1 = [theta1,c(1+8*(i-1):8*i,1)'*beta];
        theta2 = [theta2,c(1+8*(i-1):8*i,2)'*beta];

        dtheta1 = [dtheta1,c(1+8*(i-1):8*i,1)'*dbeta];
        dtheta2 = [dtheta2,c(1+8*(i-1):8*i,2)'*dbeta];
        
        ddtheta1 = [ddtheta1,c(1+8*(i-1):8*i,1)'*ddbeta];
        dddtheta1 = [dddtheta1,c(1+8*(i-1):8*i,1)'*dddbeta];
        ddddtheta1 = [ddddtheta1,c(1+8*(i-1):8*i,1)'*ddddbeta];
        
        T01 = [cos(theta1(end)),-sin(theta1(end)),0,0;sin(theta1(end)),cos(theta1(end)),0,0;0,0,1,0;0,0,0,1];
        T12 = [cos(theta2(end)),-sin(theta2(end)),0,l(1);sin(theta2(end)),cos(theta2(end)),0,0;0,0,1,0;0,0,0,1];
        T23 = [1,0,0,l(2);0,1,0,0;0,0,1,0;0,0,0,1];
        x3 = [x3,T01*T12*T23*x0];
        x2 = [x2,T01*T12*x0];
        ve = [-l(1)*sin(theta1(end))*dtheta1(end)-l(2)*sin(theta1(end)+theta2(end))*(dtheta1(end)+dtheta2(end));
              l(1)*cos(theta1(end))*dtheta1(end)+l(2)*cos(theta1(end)+theta2(end))*(dtheta1(end)+dtheta2(end));
              0];
        v = [v,norm(ve)];
    end
end


plot(plotT,theta1);
plot(plotT,dtheta1);
plot(plotT,ddtheta1);
plot(plotT,dddtheta1);
plot(plotT,ddddtheta1);
plot(plotT,theta2);
plot(plotT,v);
legend('theta1','v theta1','a theta1','j theta1','s theta1' ,'theta2','vofEE');


figure;

hold on
    O = getparam('obs_O');
    r = getparam('obs_r');
Pi = 0:pi/100:pi*2;
x = cos(Pi)*r+O(1);
y = sin(Pi)*r+O(2);
plot(x,y);

            O = [-1.34,0]';
            r = 0.2;
            Pi = 0:pi/100:pi*2;
x = cos(Pi)*r+O(1);
y = sin(Pi)*r+O(2);
plot(x,y,'-');

plot(x3(1,:),x3(2,:));
hold on
plot(x2(1,:),x2(2,:));
for i = 1:size(x2,2)
    x = [0,x2(1,i)];
    y = [0,x2(2,i)];
    line(x,y);
    drawnow

    x = [x2(1,i),x3(1,i)];
    y = [x2(2,i),x3(2,i)];
    line(x,y);
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g_ini] = MINCO(q_ini)
    
    l = [1,0.5];
    theta_num = 2;
    piece_num = 5;
    
    q = zeros(piece_num-1,theta_num);
    q(:,1) = q_ini(1:piece_num-1);
    q(:,2) = q_ini(piece_num:2*piece_num-2);
    t = virtualT2realT(q_ini(end));

    
    %初始状态
    p0 = [0,-pi/2];
%     pe = [pi,pi/2];
    pe = q_ini(2*piece_num-1:2*piece_num)';
    
    %更新M b
    M = zeros(8*piece_num,8*piece_num);
    b = zeros(8*piece_num,2);
    M(1:4,1:8) = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 2 0 0 0 0 0 ;0 0 0 6 0 0 0 0];
    b(1:4,:) = [p0;0,0;0,0;0,0];
    M_part = [1,t,t^2,t^3,t^4,t^5,t^6,t^7;
              1,t,t^2,t^3,t^4,t^5,t^6,t^7;
              0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6;
              0,0,2,6*t,12*t^2,20*t^3,30*t^4,42*t^5;
              0,0,0,6,24*t,60*t^2,120*t^3,210*t^4;
              0,0,0,0,24,120*t,360*t^2,840*t^3;
              0,0,0,0,0,120,720*t,2520*t^2;
              0,0,0,0,0,0,720,5040*t];
    for i=1:piece_num-1
        % E
        M(5+8*(i-1):12+8*(i-1),1+8*(i-1):8+8*(i-1)) = M_part;
        % F
        M(5+8*(i-1):12+8*(i-1),9+8*(i-1):16+8*(i-1)) = [0,0,0,0,0,0,0,0;
                                                      -1,0,0,0,0,0,0,0;
                                                      0,-1,0,0,0,0,0,0;
                                                      0,0,-2,0,0,0,0,0;
                                                      0,0,0,-6,0,0,0,0;
                                                      0,0,0,0,-24,0,0,0;
                                                      0,0,0,0,0,-120,0,0;
                                                      0,0,0,0,0,0,-720,0];
        b(5+8*(i-1):12+8*(i-1),:) = [q(i,:);0,0;0,0;0,0;0,0;0,0;0,0;0,0];
    end
    M(end-3:end,end-7:end) = [M_part(2,:);M_part(3,:);M_part(4,:);M_part(5,:)];
    b(end-3:end,:) = [pe;0,0;0,0;0,0];

    c = M\b;
%     c
    [f,Grad,GradT] = get_costandGandt_obs(c,t,l);
%     error("stop!!!!!!")
%     GradT
%     G
    G = M'\Grad;
%     gradp
    g = zeros(size(q));
    for i=1:piece_num-1
        g(i,:) = G(8*(i-1)+5,:);
    end
    g_end = G(8*(piece_num-1)+5,:);
    
    theta1 = pe(1);
    theta2 = pe(2);
            T01 = [cos(theta1),-sin(theta1),0,0;
                   sin(theta1),cos(theta1),0,0;
                   0,0,1,0;
                   0,0,0,1];
            T12 = [cos(theta2),-sin(theta2),0,l(1);
                   sin(theta2),cos(theta2),0,0;
                   0,0,1,0;
                   0,0,0,1];            
            dT01 = [-sin(theta1),-cos(theta1),0,0;
                   cos(theta1),-sin(theta1),0,0;
                   0,0,0,0;
                   0,0,0,0];
            dT12 = [-sin(theta2),-cos(theta2),0,0;
                   cos(theta2),-sin(theta2),0,0;
                   0,0,0,0;
                   0,0,0,0];
            T23 = [1,0,0,l(2);0,1,0,0;0,0,1,0;0,0,0,1];
            pp = [0,0,0,1]';
            p = T01*T12*T23*pp;
            
            O = [-1.34,0]';
            r = 0.2;
            
            weight = 10000;
            cost = abs(norm(p(1:2)-O)-r);
            
            fp_p = p(1:2)-O;
            if norm(p(1:2)-O)-r<0
                fp_p = - fp_p;
            end
            f = f + weight * cost;
            fp_p = weight * fp_p;
            g1 = [fp_p;0;0]' * dT01*T12*pp;
            g2 = [fp_p;0;0]' * T01*dT12*pp;
            
            
            
    
    
    g_ini = zeros(size(q_ini));
    g_ini(1:piece_num-1) = g(:,1);
    g_ini(piece_num:2*piece_num-2) = g(:,2);
    g_ini(2*piece_num-1:2*piece_num) = g_end'+[g1,g2]';

%     gradT
    % 时间正则项权重
    time_weight = 10;
    gradT = zeros(size(GradT));
    partialEi_ti_part = [ 0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6;
                        0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6;
                        0,0,2,6*t,12*t^2,20*t^3,30*t^4,42*t^5;
                        0,0,0,6,24*t,60*t^2,120*t^3,210*t^4;
                        0,0,0,0,24,120*t,360*t^2,840*t^3;
                        0,0,0,0,0,120,720*t,2520*t^2;
                        0,0,0,0,0,0,720,5040*t;
                        0,0,0,0,0,0,0,5040];
    for i=1:size(GradT,1)-1
       partialEi_ti = partialEi_ti_part;
        gradT(i) = GradT(i) - trace( (G(8*(i-1)+5:8*(i-1)+12,:))' * partialEi_ti * c(8*(i-1)+1:8*i,:) ); 
    end
    partialEi_ti = [partialEi_ti_part(2,:);partialEi_ti_part(3,:);partialEi_ti_part(4,:);partialEi_ti_part(5,:)];
    gradT(end) = GradT(end) - trace( G(end-3:end,:)' * partialEi_ti * c(end-7:end,:));
    
    gradsumT = sum(gradT) ;
    [gradt,cost_vt] = grad_real2vir(t,gradsumT,time_weight);
    g_ini(end) = gradt;
    f = f+cost_vt;
    
    
%     piece_t = t * ones(piece_num,1);
%     [gradt,cost_vt] = grad_real2vir(piece_t,gradT);
%     gradTT = sum(gradt);
%     g_ini(end) = gradTT;
%     f = f+cost_vt;
%     g_ini(end) = 0;
%     f = f;
    
end

function t = realT2virtualT(T)
    t = zeros(size(T));
    for i=1:size(T,1)
       if T>1
           t(i) = sqrt(2*T(i)-1)-1;
       else
           t(i) = 1 - sqrt(2/T(i)-1);
       end
    end
end

function T = virtualT2realT(t)
    T = zeros(size(t));
    for i=1:size(t,1)
        if t>0
            T(i) = t(i)*(t(i)/2+1)+1;
        else
            T(i) = 2/(t(i)*(t(i)-2)+2);
        end
    end
end

function [gradt,cost_vt] = grad_real2vir(RT,gradT,time_weight)

    
    vt = realT2virtualT(RT);
    gradt = zeros(size(vt));
    for i=1:size(vt,1)
       if vt(i)>0
           gradrt2vt = vt(i)+1;
       else
           denSqrt= (0.5*vt(i)-1)*vt(i)+1.0;
           gradrt2vt = (1-vt(i))/(denSqrt*denSqrt);
       end
       gradt(i) = (gradT(i)+time_weight)*gradrt2vt;
    end
    cost_vt = sum(RT)*time_weight;
end



