function output = getparam(input_arg)
%getparam 获得参数
%       if strcmp(input_arg,'piece_num')
%       if strcmp(input_arg,'theta_num')
%       if strcmp(input_arg,'l')
%       if strcmp(input_arg,'p0')
%       if strcmp(input_arg,'pe')
%       if strcmp(input_arg,'q0')
%       if strcmp(input_arg,'T0')
% 
%     piece_num_ = 5;
%     theta_num_ = 5;
%     l_ = [0.1 , 1 , 0.8 , 0.4 , 0.2];
%     p0_ = [-pi/2,-pi/2,-pi/2,-pi/2,-pi/2];
%     pe_ = [pi/2,pi/2,pi/2,pi/2,pi/2];
%     T0_ = 1;
%     q0_ = [-1*pi/3;0;pi/6;pi/3;
%           -1*pi/3;0;pi/6;pi/3;
%           -1*pi/3;0;pi/6;pi/3;
%           -1*pi/3;0;pi/6;pi/3;
%           -1*pi/3;0;pi/6;pi/3];
    obs_O_ = [1.1,1.1]';
    obs_r_ = [0.2]';
    
    if strcmp(input_arg,'piece_num')
        output = piece_num_;
    elseif strcmp(input_arg,'theta_num')
        output = theta_num_;
    elseif strcmp(input_arg,'l')
        output = l_;
    elseif strcmp(input_arg,'p0')
        output = p0_;
    elseif strcmp(input_arg,'pe')
        output = pe_;
    elseif strcmp(input_arg,'T0')
        output = T0_;
    elseif strcmp(input_arg,'q0')
        output = q0_;
    elseif strcmp(input_arg,'obs_O')
        output = obs_O_;
    elseif strcmp(input_arg,'obs_r')
        output = obs_r_;
    else
        error('输入参数错误');
    end
end

