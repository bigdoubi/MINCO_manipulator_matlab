clear
close all

% 前端采样点
initroutewithtime = importdata("initroutewithtime.txt");
frondend_num = size(initroutewithtime, 1)-2;
% 前端传到后端的点
originroute = importdata("originroute.txt");
niheroute = importdata("niheroute.txt");

% 获得M矩阵
piece_num = size(originroute, 1)-1;
t = originroute(2,1);
M = zeros(8*piece_num,8*piece_num);
M(1:4,1:8) = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 ;0 0 2 0 0 0 0 0 ;0 0 0 6 0 0 0 0];

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
end
M(end-3:end,end-7:end) = [M_part(2,:);M_part(3,:);M_part(4,:);M_part(5,:)];

M_inv = inv(M);

B = zeros(8*piece_num, piece_num);
B(1, piece_num) = originroute(1,2);
for i=1:piece_num-1
   B(4+8*(i-1)+1, i) = 1; 
end
B(5+8*(piece_num-1), piece_num) = originroute(end,2);

beta = zeros(frondend_num, 8*piece_num);
for i=2:frondend_num-1
    index = floor(initroutewithtime(i,3)/t);
    beta(i, 8*index+1:8*index+8) = M_part(1,:);
end

% 定义最小二乘问题
A = beta * M_inv * B;
A11  = A(:,1:end-1);
A12  = A(:,end);
c = -A12 + initroutewithtime(2:end-1, 1);

overline_b = A11\c;

% 获得最新的曲线
b = zeros(8*piece_num,1);
b(1:4,:) = [originroute(1,2);0;0;0];
for i=1:piece_num-1
    b(5+8*(i-1):12+8*(i-1),:) = [overline_b(i);0;0;0;0;0;0;0];
end
b(end-3:end,:) = [originroute(end,2);0;0;0];
c = M_inv * b;

theta = [];
for ttt = 0:0.001:originroute(end, 1)
   index = min([floor(ttt/t), piece_num-1]);
   index_t = ttt-index*t;
   betat = [1,index_t,index_t^2,index_t^3,index_t^4,index_t^5,index_t^6,index_t^7];
   theta = [theta;betat * c(8*index+1:8*index+8)];
end

plot(0:0.001:originroute(end, 1), theta);
hold on 
plot(initroutewithtime(:,3), initroutewithtime(:,1),'x');
plot(niheroute(:,1),niheroute(:,2),'--');