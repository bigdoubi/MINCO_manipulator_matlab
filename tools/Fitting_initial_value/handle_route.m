clear
close all

originroute = importdata("originroute.txt");
niheroute = importdata("niheroute.txt");
originroutexy = importdata("originroutexy.txt");
finalroute = importdata("finalroute.txt");

plot(originroute(:,1),originroute(:,2))
hold on
plot(niheroute(:,1), niheroute(:,2))
plot(originroute(:,1),originroute(:,2),'ob')
% line([0,12],[2*pi, 2*pi])
% line([0,12],[3*pi, 3*pi])
% line([0,12],[pi, pi])
max_yaw = max([niheroute(:,2);originroute(:,2)]);
min_yaw = min([niheroute(:,2);originroute(:,2)]);
for i=floor(min_yaw/pi):ceil(max_yaw/pi)
    line([niheroute(1,1),niheroute(end,1)],[i*pi, i*pi])
end

figure();
plot(originroute(:,1),originroute(:,3))
hold on
plot(niheroute(:,1), niheroute(:,3))
plot(originroute(:,1),originroute(:,3),'ob')

figure();
subplot(1,2,1);
plot(finalroute(:,1),finalroute(:,3));
subplot(1,2,2);
plot(finalroute(:,1),finalroute(:,2));

