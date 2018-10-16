%----------------------Figures---------------------------------
figure(1)
for count=1:length(tkeep),
    surf(x,y,squeeze(ukeep(:,:,count))'); view(-30,70); shading('interp')
    xlabel('x');ylabel('y');zlabel('u')
    axis([-L L -L L -1 1]);
    Fmovie(count)=getframe;
end
movie(Fmovie)
figure(2)
for count=1:length(tkeep),
    surf(x,y,squeeze(vkeep(:,:,count))'); view(-30,70); shading('interp');
    xlabel('x');ylabel('y');zlabel('v')
    axis([-L L -L L -1 1]);
    Gmovie(count)=getframe;
end
movie(Gmovie)
figure(3)
subplot(2,1,1)
plot(x,squeeze(ukeep(:,find(y==0),:))')
xlabel('x'); ylabel('u');
subplot(2,1,2)
plot(y,squeeze(vkeep(:,find(y==0),:))')
xlabel('x'); ylabel('v');
figure(4)
tplot=[200, 400, 600, 1000];
for count=1:4,
    subplot(2,2,count)
    surf(x,y,squeeze(ukeep(:,:,find(tkeep==tplot(count))))'); view(-30,70); shading('interp')
    xlabel('x');ylabel('y');zlabel('u')
    axis([-L L -L L -1 1]);
end
