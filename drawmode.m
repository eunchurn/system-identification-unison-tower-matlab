function me=drawmode(M30,M58,M80,len,figname,vis,filename)
me=1;
try
    %     len=10;
    inix=[-len/2 len/2 len/2 -len/2 -len/2];
    iniy=[-len/2 -len/2 len/2 len/2 -len/2];
    iniz1=[0 0 0 0 0];
    iniz2=[30 30 30 30 30];
    iniz3=[58 58 58 58 58];
    iniz4=[80 80 80 80 80];
    mode30x=inix+M30(2);
    mode30y=[iniy(1:2)+M30(1) iniy(3:4)+M30(3) iniy(5)+M30(1)];
    mode58x=inix+M58(2);
    mode58y=[iniy(1:2)+M58(1) iniy(3:4)+M58(3) iniy(5)+M58(1)];
    mode80x=inix+M80(2);
    mode80y=[iniy(1:2)+M80(1) iniy(3:4)+M80(3) iniy(5)+M80(1)];
    inivz=[0 80 0 80 0 80 0 80];
    figure('Color',[1 1 1],'Position',[10 10 480 800],'Name',figname,'Visible',vis)
    line([inix inix inix inix],[iniy iniy iniy iniy],[iniz1 iniz2 iniz3 iniz4],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':')
    hold on
    for kk=1:4
        line([inix(kk) inix(kk)],[iniy(kk) iniy(kk)],[iniz1(kk) iniz4(kk)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
    end
    line([inix mode30x mode58x mode80x],[iniy mode30y mode58y mode80y],[iniz1 iniz2 iniz3 iniz4],'Color','r','LineWidth',1.5)
    for kk=1:4
        line([inix(kk) mode30x(kk) mode58x(kk) mode80x(kk)],[iniy(kk) mode30y(kk) mode58y(kk) mode80y(kk)],[iniz1(kk) iniz2(kk) iniz3(kk) iniz4(kk)],'Color','r','LineWidth',1.5);
    end
    title(figname);
    view(gca,[-27.5 10])
    xlim([-(len/2+1) len/2+1]),ylim([-(len/2+1) len/2+1]),zlim([0 85])
    text(0,0,100,figname)
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    set(gca,'ZTick',[])
    set(gca,'Visible','off')
    set(gcf,'PaperPositionMode','auto')
    print(gcf,'-dtiff',filename)
catch me
end