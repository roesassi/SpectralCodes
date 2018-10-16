%----------------------Figures---------------------------------
figure(1)
for count=1:length(tkeep),
    surf(x,y,squeeze(ukeep(:,:,count))'); view(-30,70); shading('interp')
    xlabel('x');ylabel('y');zlabel('u')
    axis([-L L -L L 0 1]);
    Fmovie(count)=getframe;
end
movie(Fmovie)
figure(2)
subplot(2,1,1)
plot(x,squeeze(ukeep(:,find(y==0),:))')
xlabel('x'); ylabel('u'); axis([-L L 0 1])
subplot(2,1,2)
plot(y,squeeze(ukeep(find(x==0),:,:))')
xlabel('y'); ylabel('u'); axis([-L L 0 1])
