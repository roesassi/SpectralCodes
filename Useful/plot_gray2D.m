%----------------------Figures---------------------------------
figure(1)
for count=1:length(tkeep),
    clf
    surf(x,y,squeeze(ukeep(:,:,count))'); view(-30,70); shading('interp')
    xlabel('x');ylabel('y');zlabel('u')
    axis([-L L -L L 0 1]);
    Fmovie(count)=getframe;
end
movie(Fmovie)
figure(2)
for count=1:length(tkeep),
    clf
    surf(x,y,squeeze(vkeep(:,:,count))'); view(-30,70); shading('interp');
    xlabel('x');ylabel('y');zlabel('v')
    axis([-L L -L L 0 4]);
    Gmovie(count)=getframe;
end
movie(Gmovie)
figure(3)
subplot(2,1,1)
plot(x,squeeze(ukeep(:,find(y==0),:))')
xlabel('x'); ylabel('u'); axis([-L L 0 1])
subplot(2,1,2)
plot(y,squeeze(vkeep(:,find(y==0),:))')
xlabel('x'); ylabel('v'); axis([-L L 0 4])
figure(4)
subplot(2,2,1)
surf(x,y,squeeze(ukeep(:,:,end))'); view(-30,70); shading('interp')
xlabel('x');ylabel('y');zlabel('u')
axis([-L L -L L 0 1]);
subplot(2,2,2)
surf(x,y,squeeze(vkeep(:,:,end))'); view(-30,70); shading('interp');
xlabel('x');ylabel('y');zlabel('v')
axis([-L L -L L 0 4]);
