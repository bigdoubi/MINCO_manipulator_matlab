function [cost,G,Gt] = get_costandGandt_obs(c,t,l)
    
    smooth_weight = 0.0001;
    range_weight = 1000000;
    obs_weight   = 1000000;
    v_weight = 10;

    sampling_num = 20;
    step = t/sampling_num;
    t = ones(size(c,1)/8,1)*t;

    cost = 0;
    G = zeros(size(c)); 
    Gt = zeros(size(t));
    
    cost1 = 0;
    cost2 = 0;
    cost3 = 0;
    G1 = zeros(size(c));        
    
    cost_t1 = 0;
    G_t1 = zeros(size(t));
    % minimum snap
    
    
% 100800*c7^2*t^7 + 100800*c6*c7*t^6 + (25920*c6^2 + 40320*c5*c7)*t^5 + (10080*c4*c7 + 21600*c5*c6)*t^4 + (4800*c5^2 + 5760*c4*c6)*t^3 + 2880*c4*c5*t^2 + 576*c4^2*t
% 705600*c7^2*t^6 + 604800*c6*c7*t^5 + (129600*c6^2 + 201600*c5*c7)*t^4 + (40320*c4*c7 + 86400*c5*c6)*t^3 + (14400*c5^2 + 17280*c4*c6)*t^2 + 5760*c4*c5*t + 576*c4^2
    for i = 1:size(t,1)
        % 把这里写完，对比一下cost G有没有求错
        cost1 = cost1 + 100800*dot(c(8*(i-1)+8,:),c(8*(i-1)+8,:))*t(i)^7 ...
                      + 100800*dot(c(8*(i-1)+7,:),c(8*(i-1)+8,:))*t(i)^6 ...
                      + (25920*dot(c(8*(i-1)+7,:),c(8*(i-1)+7,:)) + 40320*dot(c(8*(i-1)+6,:),c(8*(i-1)+8,:)))*t(i)^5 ...
                      + (10080*dot(c(8*(i-1)+5,:),c(8*(i-1)+8,:)) + 21600*dot(c(8*(i-1)+6,:),c(8*(i-1)+7,:)))*t(i)^4 ...
                      + (4800*dot(c(8*(i-1)+6,:),c(8*(i-1)+6,:)) + 5760*dot(c(8*(i-1)+5,:),c(8*(i-1)+7,:)))*t(i)^3 ...
                      + 2880*dot(c(8*(i-1)+5,:),c(8*(i-1)+6,:))*t(i)^2 ...
                      + 576*dot(c(8*(i-1)+5,:),c(8*(i-1)+5,:))*t(i);
        G1(8*(i-1)+5:8*(i-1)+8,:) = G1(8*(i-1)+5:8*(i-1)+8,:) + ...
                                                             [2*576*c(8*(i-1)+5,:)*t(i) + 2880*c(8*(i-1)+6,:)*t(i)^2 + 5760*c(8*(i-1)+7,:)*t(i)^3 + 10080*c(8*(i-1)+8,:)*t(i)^4;
                                                              2880*c(8*(i-1)+5,:)*t(i)^2 + 2*4800*c(8*(i-1)+6,:)*t(i)^3 + 21600*c(8*(i-1)+7,:)*t(i)^4 + 40320*c(8*(i-1)+8,:)*t(i)^5;
                                                              5760*c(8*(i-1)+5,:)*t(i)^3 + 21600*c(8*(i-1)+6,:)*t(i)^4 + 2*25920*c(8*(i-1)+7,:)*t(i)^5 + 100800*c(8*(i-1)+8,:)*t(i)^6;
                                                              10080*c(8*(i-1)+5,:)*t(i)^4 + 40320*c(8*(i-1)+6,:)*t(i)^5 + 100800*c(8*(i-1)+7,:)*t(i)^6 + 2*100800*c(8*(i-1)+8,:)*t(i)^7];
        G_t1(i) = G_t1(i) + ...
                          + 705600*dot(c(8*(i-1)+8,:),c(8*(i-1)+8,:))*t(i)^6 ...
                          + 604800*dot(c(8*(i-1)+7,:),c(8*(i-1)+8,:))*t(i)^5 ...
                          + (129600*dot(c(8*(i-1)+7,:),c(8*(i-1)+7,:)) + 201600*dot(c(8*(i-1)+6,:),c(8*(i-1)+8,:)))*t(i)^4 ...
                          + (40320*dot(c(8*(i-1)+5,:),c(8*(i-1)+8,:)) + 86400*dot(c(8*(i-1)+6,:),c(8*(i-1)+7,:)))*t(i)^3 ...
                          + (14400*dot(c(8*(i-1)+6,:),c(8*(i-1)+6,:)) + 17280*dot(c(8*(i-1)+5,:),c(8*(i-1)+7,:)))*t(i)^2 ...
                          + 5760*dot(c(8*(i-1)+5,:),c(8*(i-1)+6,:))*t(i) ...
                          + 576*dot(c(8*(i-1)+5,:),c(8*(i-1)+5,:));                                                          
%         cost_try =  24*24*dot(c(8*(i-1)+5,:),c(8*(i-1)+5,:))*t(i) ...
%                     + 2*24*120*dot(c(8*(i-1)+5,:),c(8*(i-1)+6,:))*t(i)^2 / 2 ...
%                     + (120*120*dot(c(8*(i-1)+6,:),c(8*(i-1)+6,:)) + 2*24*360*dot(c(8*(i-1)+5,:),c(8*(i-1)+7,:)))*t(i)^3 / 3 ...
%                     + (120*2*360*dot(c(8*(i-1)+6,:),c(8*(i-1)+7,:)) + 24*840*2*dot(c(8*(i-1)+5,:),c(8*(i-1)+8,:)))*t(i)^4 / 4 ...
%                     + (360*360*dot(c(8*(i-1)+7,:),c(8*(i-1)+7,:)) + 120*840*2*dot(c(8*(i-1)+6,:),c(8*(i-1)+8,:)))*t(i)^5 / 5 ...
%                     + 2*360*840*dot(c(8*(i-1)+7,:),c(8*(i-1)+8,:))*t(i)^6 / 6 ...
%                     + 840*840*dot(c(8*(i-1)+8,:),c(8*(i-1)+8,:))*t(i)^7 / 7
%         G1(8*(i-1)+5:8*(i-1)+8,:) = G1(8*(i-1)+5:8*(i-1)+8,:) + ...
%                                                             [2*24*24*c(8*(i-1)+5)*t(i) + 2*24*120*c(8*(i-1)+6)*t(i)^2/2 + 2*24*360*c(8*(i-1)+7)*t(i)^3/3 + 24*840*2*c(8*(i-1)+8)*t(i)^4/4;
%                                                             2*24*120*c(8*(i-1)+5)*t(i)^2/2 + 2*120*120*c(8*(i-1)+6)*t(i)^3/3 + 120*2*360*c(8*(i-1)+7)*t(i)^4/4 + 120*840*2*c(8*(i-1)+8)*t(i)^5/5;
%                                                             2*24*360*c(8*(i-1)+5)*t(i)^3/3 + 120*2*360*c(8*(i-1)+6)*t(i)^4/4 + 2*360*360*c(8*(i-1)+7)*t(i)^5/5 + 2*360*840*c(8*(i-1)+8)*t(i)^6/6;
%                                                             24*840*2*c(8*(i-1)+5)*t(i)^4/4 + 120*840*2*c(8*(i-1)+6)*t(i)^5/5 + 2*360*840*c(8*(i-1)+7)*t(i)^6/6 + 2*840*840*c(8*(i-1)+8)*t(i)^7/7];         
%         G_t1(i) = G_t1(i) + ...
%                         24*24*dot(c(8*(i-1)+5,:),c(8*(i-1)+5,:)) ...
%                     + 2*24*120*dot(c(8*(i-1)+5,:),c(8*(i-1)+6,:))*t(i) ...
%                     + 120*120*dot(c(8*(i-1)+6,:),c(8*(i-1)+6,:))*t(i)^2 + 2*24*360*dot(c(8*(i-1)+5,:),c(8*(i-1)+7,:))*t(i)^2 ...
%                     + 120*2*360*dot(c(8*(i-1)+6,:),c(8*(i-1)+7,:))*t(i)^3 + 24*840*2*dot(c(8*(i-1)+5,:),c(8*(i-1)+8,:))*t(i)^3 ...
%                     + (360*360*dot(c(8*(i-1)+7,:),c(8*(i-1)+7,:)) + 120*840*2*dot(c(8*(i-1)+6,:),c(8*(i-1)+8,:)))*t(i)^4 ...
%                     + 2*360*840*dot(c(8*(i-1)+7,:),c(8*(i-1)+8,:))*t(i)^5 ...
%                     + 840*840*dot(c(8*(i-1)+8,:),c(8*(i-1)+8,:))*t(i)^6;
    end
    
    G1 = smooth_weight * G1;
    cost1 = smooth_weight * cost1;
    G_t1 = smooth_weight * G_t1;

    G2 = zeros(size(c));
    G3 = zeros(size(c));
    cost2 = 0;
    cost3 = 0;
    
    G_t2 = zeros(size(t));
    G_t3 = zeros(size(t));
    
    G4 = zeros(size(c));
    cost4 = 0;
    G_t4 = zeros(size(t));
    
    % 采样的方式计算range梯度和obs梯度
    for i = 1:size(t,1)
        for j = 0:1:sampling_num
            if(j==0||j==sampling_num)
                w = 1/2;
            else
                w = 1;
            end
            alpha = j/sampling_num;
            
            T = t(i)/sampling_num*j;
            C = c(8*(i-1)+1:8*i,:);
            beta = [1,T^1,T^2,T^3,T^4,T^5,T^6,T^7]';
            dbeta = [0,1,2*T,3*T^2,4*T^3,5*T^4,6*T^5,7*T^6]';
            ddbeta = [0,0,2,6*T,12*T^2,20*T^3,30*T^4,42*T^5]';
            theta1 = C(:,1)'*beta;
            theta2 = C(:,2)'*beta;
            dtheta1 = C(:,1)'*dbeta;
            dtheta2 = C(:,2)'*dbeta;
            ddtheta1 = C(:,1)'*ddbeta;
            ddtheta2 = C(:,2)'*ddbeta;
            % 将关节角度限制到一个范围内
            beta = [1,T^1,T^2,T^3,T^4,T^5,T^6,T^7]';
            costr = 0;
            Gr = zeros(size(C));           
            
            cost_th1 = (theta1-pi/2)^2 - (pi/2)^2;
            cost_th2 = (theta2)^2 - (pi/2)^2;            
            if(cost_th1 > 0)
                [f,df] = smoothL1(cost_th1);
                Gr(:,1) = 2*(theta1-pi/2)*beta*df;
                costr = costr + f;
                Gr_t = 2*(theta1-pi/2)*theta1;
                G_t2(i) = G_t2(i) + w*range_weight*(df*Gr_t*step*alpha + f/sampling_num);
            end
            if(cost_th2 > 0)
                [f,df] = smoothL1(cost_th2);
                Gr(:,2) = 2*(theta2)*beta*df;
                costr = costr + f; 
                Gr_t = 2*dtheta2*theta2;
                G_t2(i) = G_t2(i) + w*range_weight*(df*Gr_t*step*alpha + f/sampling_num);
            end
            cost2 = cost2 + costr*range_weight*w*step;
            G2(1+8*(i-1):8*i,:) = G2(1+8*(i-1):8*i,:) + Gr*range_weight*w*step;
            
            %速度约束
            vmax = 2.5;
            costv = 0;
            Gv = zeros(size(C));
            
            cost_v1 = dtheta1^2 - vmax^2;
            cost_v2 = dtheta2^2 - vmax^2;
            if(cost_v1 > 0)
                [f,df] = smoothL1(cost_v1);
                Gv(:,1) = 2*dtheta1*dbeta*df;
                costv = costv + f;
                Gv_t = 2*ddtheta1*dtheta1;
                G_t4(i) = G_t4(i) + w*v_weight*(df*Gv_t*step*alpha + f/sampling_num);
            end
            if(cost_v2 > 0)
                [f,df] = smoothL1(cost_v2);
                Gv(:,2) = 2*dtheta2*dbeta*df;
                costv = costv + f;
                Gv_t = 2*ddtheta2*dtheta2;
                G_t4(i) = G_t4(i) + w*v_weight*(df*Gv_t*step*alpha + f/sampling_num);
            end
            cost4 = cost4 + costv*v_weight*w*step;
            G4(1+8*(i-1):8*i,:) = G4(1+8*(i-1):8*i,:) + Gv*v_weight*w*step;
                
            
            %碰撞约束
            costo = 0;
            Go = zeros(size(C)); 
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
            dT12 = [-sin(theta2),-cos(theta2),0,0;%存疑 是否应该是l(1)
                   cos(theta2),-sin(theta2),0,0;
                   0,0,0,0;
                   0,0,0,0];     
            
                     
%             R01 = T01(1:3,1:3);
%             R12 = T12(1:3,1:3);
%             dR01 = dT01(1:3,1:3);
%             dR12 = dT12(1:3,1:3);
            
            interval = 0.1;
            for k=0:interval:1
                % 第一个机臂
%                 if k==0||k==1
%                     wk = 0.5;
%                 else
%                     wk = 1;
%                 end
                if k == 0
                    continue;
                end
                
                
                Pp = [k*l(1),0,0,1]';
                p = T01*Pp;
                [cost_o1,Gp1] = ESDF_costandG(p(1:2));
                Gp1 = [-Gp1;0;1];
                edge = 0.01;
                if cost_o1 < edge  
                    [f,df] = smoothL1(-(cost_o1-edge));
                    costo = costo + f;
                    Go(:,1) = Go(:,1) + df * Gp1' * dT01 * Pp * beta;                    
                    Go_t = Gp1' * dT01 * Pp * dtheta1;
                    G_t3(i) = G_t3(i) + w * obs_weight*(df*Go_t*step*alpha + f/sampling_num);
%                     costo = costo * wk + f;
%                     Go(:,1) = Go(:,1) + df * Gp1' * dT01 * Pp * beta * wk * w * step * obs_weight;
%                     Go_t = Gp1' * dT01 * Pp * dtheta1;
%                     G_t3(i) = G_t3(i) + wk * w * obs_weight * (df*Go_t*step + f/sampling_num);
                end
                % 第二个机臂
                Pp = [k*l(2),0,0,1]';
                p = T01*T12*Pp;
                [cost_o2,Gp2] = ESDF_costandG(p(1:2));
                Gp2 = [-Gp2;0;0];
                if cost_o2 < edge  
                    [f,df] = smoothL1(-(cost_o2-edge));
                    costo = costo + f;
                    Go(:,1) = Go(:,1) + df * Gp2' * dT01 * T12 * Pp * beta;
                    Go(:,2) = Go(:,2) + df * Gp2' * T01 * dT12 * Pp * beta;
                    Go_t = Gp2' * dT01 * T12 * Pp * dtheta1 + Gp2' * T01 * dT12 * Pp * dtheta2;
                    G_t3(i) = G_t3(i) + w * obs_weight * (df*Go_t*step*alpha + f/sampling_num);
                end
            end
            cost3 = cost3 + costo*obs_weight*step*w;
            G3(1+8*(i-1):8*i,:) = G3(1+8*(i-1):8*i,:) + Go * w * step * obs_weight;
%             G3(1+8*(i-1):8*i,:)
        end
    end
    G = G1+G2+G3+G4;
    cost = cost1+cost2+cost3+cost4;  
    Gt = G_t1+G_t2+G_t3+G_t4;
%     G = G1;
%     cost = cost1;
%     Gt = G_t1;
end




% function [costr,Gr_t,Gr] = get_pieceCG_range(C,T)
%     beta = [1,T^1,T^2,T^3,T^4,T^5,T^6,T^7]';
% %     dbeta = [0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6]';
%     theta1 = C(:,1)'*beta;
%     theta2 = C(:,2)'*beta;
%     costr = 0;
%     Gr = zeros(size(C));
%     Gr_t = 0;
%     
%     cost_th1 = (theta1-pi/2)^2 - (pi/2)^2;
%     cost_th2 = (theta2)^2 - (pi/2)^2;
% 
%     if(cost_th1 > 0)
%         [f,df] = smoothL1(cost_th1);
%         Gr(:,1) = Gr(:,1) + 2*(theta1-pi/2)*beta*df;
%         costr = costr + f;
%     end
%     if(cost_th2 > 0)
%         [f,df] = smoothL1(cost_th2);
%         Gr(:,2) = Gr(:,2) + 2*(theta2)*beta*df;
%         costr = costr + f; 
%     end  
% 
% end


% function [costo,Go] = get_pieceCG_obs(C,T,l)
%     beta = [1,T^1,T^2,T^3,T^4,T^5,T^6,T^7]';
% %     dbeta = [0,1,2*t,3*t^2,4*t^3,5*t^4,6*t^5,7*t^6]';
%     theta1 = C(:,1)'*beta;
%     theta2 = C(:,2)'*beta;
%     costo = 0;
%     Go = zeros(size(C));  
% 
% 
%     T01 = [cos(theta1),-sin(theta1),0,0;
%            sin(theta1),cos(theta1),0,0;
%            0,0,1,0;
%            0,0,0,1];
%     T12 = [cos(theta2),-sin(theta2),0,l(1);
%            sin(theta2),cos(theta2),0,0;
%            0,0,1,0;
%            0,0,0,1];
% %     T23 = [1,0,0,l(2);
% %            0,1,0,0;
% %            0,0,1,0;
% %            0,0,0,1];
% 
%     dT01 = [-sin(theta1),-cos(theta1),0,0;
%            cos(theta1),-sin(theta1),0,0;
%            0,0,1,0;
%            0,0,0,1];
%     dT12 = [-sin(theta2),-cos(theta2),0,l(1);
%            cos(theta2),-sin(theta2),0,0;
%            0,0,1,0;
%            0,0,0,1];
%     step = 0.1;
%     for i=0:step:1
%         % 第一个机臂
%         if i==0||i==1
%             w = 0.5;
%         else
%             w = 1;
%         end
%         Pp = [i*l(1),0,0,1]';
%         p = T01*Pp;
%         [cost1,Gp1] = ESDF_costandG(p(1:2));
%         Gp1 = [-Gp1;0;0];
% 
%         edge = 0.05;
% 
%         if cost1 < edge  %%正负号存疑
% %             cost1
%             [f,df] = smoothL1(-(cost1-edge));
%             costo = costo + f;
%             Go(:,1) = Go(:,1) + df * Gp1' * dT01 * Pp * beta * w * step;
%         end
%         % 第二个机臂
%         Pp = [i*l(2),0,0,1]';
%         p = T01*T12*Pp;
% %         plot(p(1),p(2),'x');
% %         hold on;
%         [cost2,Gp2] = ESDF_costandG(p(1:2));
%         Gp2 = [-Gp2;0;0];
% %         cost2
%         if cost2<edge  %%正负号存疑
% %             cost2
%             [f,df] = smoothL1(-(cost2-edge));
%             costo = costo + f;
%             Go(:,1) = Go(:,1) + df * Gp2' * dT01 * T12 * Pp * beta * w * step;
%             Go(:,2) = Go(:,2) + df * Gp2' * T01 * dT12 * Pp * beta * w * step;
%         end
%     end
% end

function [cost,fp_p] = ESDF_costandG(p)
    % 障碍物是一个圆，圆心坐标[1,1],半径0.3  没啥别的意思，纯粹为了ESDF好算
    O = getparam('obs_O');
    r = getparam('obs_r');
    cost = norm(p-O)-r;
    fp_p = p-O;
end


function [f,df] = smoothL1(x)
    pe = 0.01;
    half = 0.5 * pe;
    f3c = 1.0 / (pe * pe);
    f4c = -0.5 * f3c / pe;
    d2c = 3.0 * f3c;
    d3c = 4.0 * f4c;

    if (x < pe)
        f = (f4c * x + f3c) * x * x * x;
        df = (d3c * x + d2c) * x * x;
    else
        f = x - half;
        df = 1.0;
    end
end

    
    
