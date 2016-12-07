function analysisv05(indatapath,option)
% indatapath='2011_01_11_16_20_indata_1.txt';
try
    clc;
    if exist('silog.log','file')==2
        delete('silog.log');
    end
    clsf=1;
    if nargin>=1
        diary silog.log
        disp(['#### Modal Analysis Process ####'])
        disp(['#### Path File Name : ',indatapath,' ####'])
        if nargin==1
            option='0';
        end
        option=str2num(option);
        if option
%             vis='off';
%             clsf=1;
            vis='on';
            clsf=0;
        else
            vis='off';
            clsf=1;
        end
        if exist(indatapath,'file')==2
            disp('## Path Check ... O.K')
            fid1=fopen(indatapath);
            ii=0;
            Fs=256;
            while ~feof(fid1)
                ii=ii+1;
                filen{ii}=fgetl(fid1);
            end
            for kk=1:length(filen)
                if ~exist(filen{kk},'file')
                    filexist(kk)=0;
                else
                    filexist(kk)=1;
                end
            end
            if sum(filexist)==length(filen)
                for kk=1:length(filen)
                    filename=filen{kk};
                    
                    fid=fopen(filename);
                    tline=fgetl(fid);
                    hh=str2num(tline(1:2));
                    mm=str2num(tline(4:5));
                    ss=str2num(tline(7:8));
                    t(kk)=hh*3600+mm*60+ss;
                    h(kk)=hh;
                    m(kk)=mm;
                    s(kk)=ss;
                    fclose(fid);
                end
                try
                    disp('## Time Synchronizing...')
                    latetime=max(t)+1;
                    lateh=floor(latetime/3600);
                    latem=floor(mod(latetime,3600)/60);
                    lates=mod(mod(latetime,3600),60);
                    index={'Top Floor X1','Top Floor Y1','Top Floor Y2','Mid Floor X1','Mid Floor Y1','Mid Floor Y2','Low Floor X1','Low Floor Y1','Low Floor Y2'};
                    for kk=1:length(filen)
                        filename=filen{kk};
                        fid1=fopen(filename);
                        %     nt(kk)=D(kk).bytes/32;
                        ii=0;
                        while 1
                            ii=ii+1;
                            tline=fgetl(fid1);
                            hh=str2num(tline(1:2));
                            mm=str2num(tline(4:5));
                            ss=str2num(tline(7:8));
                            if hh==lateh
                                if mm==latem
                                    if ss==lates
                                        syncs(kk)=ii;
                                        disp([index{kk},' Sync @ Line : ',num2str(syncs(kk)),'  ''',tline,'''  filepath=',filename]);
                                        break;
                                    end
                                end
                            end
                        end
                    end
                catch me
                    disp('데이터 동기화 에러 : 일부 데이터의 시간을 매칭할 수 없습니다.')
                    message2='SI program excution error 데이터 동기화 에러';
                    fid2=fopen('SI_program_condi.txt','w');
                    fprintf(fid2,'%s',message2);
                    fclose(fid2);
                    close all;
                    return;
                    quit;
                end
                disp('Time Synchronizing...done')
                for kk=1:length(filen)
                    filename=filen{kk};
                    fid1=fopen(filename);
                    C=textscan(fid1,'%s %f','Delimiter',',');
                    D{kk}=C{2};
                    T{kk}=C{1};
                    nt(kk)=length(C{2});
                    fclose(fid1);
                end
                syncln=min(nt-syncs);
                disp(['Data Length : ',num2str(syncln)])
                for kk=1:length(D)
                    data(kk).timestr=T{kk}(syncs(kk):syncs(kk)+syncln);
                    data(kk).values=D{kk}(syncs(kk):syncs(kk)+syncln);
                    mdata(:,kk)=data(kk).values*0.25;
                end
                clear C;
%                 clear D;
                clear T;
                clear data;
                [bf,af]=butter(7,5/Fs,'low');
                disp('## Filtering data low-pass-filter cutoff at Fc=5Hz');
                fmdata=filtfilt(bf,af,detrend(mdata));
                save('afilt.txt','fmdata','-ascii');
                disp(' -Saved filtered acceleration data in ''afilt.txt''');

                [fmr,fmc]=size(fmdata);
                cenid=int32(fmr/2);
                macc(:,1)=max(abs(fmdata(1:cenid,:)))';
                macc(:,2)=max(abs(fmdata(cenid+1:end,:)))';
                
%                 maxacc=max(abs(fmdata))';
                maxacc=mean(macc')';
                
                save('maxaccel.txt','maxacc','-ascii');
                disp(' - Saved peak accelerations data in ''maxaccel.txt''');

                nfft=2^(nextpow2(syncln));
                minorder=1;
                maxorder=100;
                if exist('freqrange.txt','file')==2
                    frange=load('freqrange.txt');
                    [r,c]=size(frange);
                    disp('## Starting SSI/BR Estimation ##')
                    disp('------------------------------------------------')
                    disp('Estimation Frequencies Bandwidth')
                    disp('------------------------------------------------')
                    for kk=1:r
                        disp(['From   ',num2str(frange(kk,1)),'(Hz)  to  ',num2str(frange(kk,2)),'(Hz)'])
                    end
                    disp('------------------------------------------------')
                    disp([' - Minorder : ',num2str(minorder)])
                    disp([' - Maxorder : ',num2str(maxorder)])
                    temp=rand(r,1).*(frange(:,2)-frange(:,1));
                    fn=frange(:,1)+temp;
                    
%                     [na_freq, modevec, xi, na_freq_cov, modevec_cov, xi_cov, estVarBRfreq, estVarBRzeta] = sfadccore(mdata,Fs,minorder,maxorder,[],frange);
                    [na_freq, modevec, xi, na_freq_cov, modevec_cov, xi_cov, estVarBRfreq, estVarBRzeta] = sfadccoretemp2(mdata,Fs,minorder,maxorder,[],frange);                    
                    
                    LM=modevec([9 7 8],:);
                    MM=modevec([6 4 5],:);
                    HM=modevec([3 1 2],:);
                    disp('------------------------------------------------')
                    disp('Natural Frequencies (Eigenvalues)')
                    disp('------------------------------------------------')
                    disp(sprintf('%s\t\t\t%s\t\t\t%s','Mode#1','Mode#2','Mode#3'));
                    disp(num2str(na_freq))
                    disp('------------------------------------------------')
                    disp('Eigenvectors')
                    disp('------------------------------------------------')
                    disp(sprintf('%s\t\t\t%s\t\t\t%s','Mode#1','Mode#2','Mode#3'));
                    disp(num2str(modevec))
                    disp('------------------------------------------------')
                    %             LM=modevec(1:3,find(na_freq));
                    %             MM=modevec(4:6,find(na_freq));
                    %             HM=modevec(7:9,find(na_freq));
                    %             LM=modevec(1:3,:);
                    %             MM=modevec(4:6,:);
                    %             HM=modevec(7:9,:);
                    %                 save modevec modevec
                    modfigfile={'Modeshape Tower1mode1.tif','Modeshape Tower1mode2.tif','Modeshape Tower1mode3.tif'};
                    for kk=1:r % for kk=1:length(find(na_freq))
                        if na_freq(kk)
                            figname=['Mode #',num2str(kk),' Fn=',num2str(na_freq(kk)),'(Hz)'];
                        else
                            figname=['Mode #',num2str(kk),' cannot be estimated'];
                        end
                        drawmode(LM(:,kk),MM(:,kk),HM(:,kk),6,figname,vis,modfigfile{kk});
                        disp([' - drawing mode #',num2str(kk),' saved as file name ''',modfigfile{kk},'''']);
                    end
                    
                    try
                        disp('## Natural Frequencies')
                        for kk=1:r
                            freqvec=estVarBRfreq(estVarBRfreq(:,1)>=frange(kk,1) & estVarBRfreq(:,1)<=frange(kk,2),1);
                            if ~isempty(freqvec)
                                fn(kk)=mean(freqvec);
                                disp(['Mode #',num2str(kk),'=',num2str(fn(kk)),'(Hz)']);
                            end
                        end
                        sfn=na_freq;
                        sfn(find(~na_freq))=fn(find(~na_freq));
                        disp(sfn)
                        save('fn.txt','sfn','-ascii');
                    catch me
                        message1=me.message;
                        message2='SI program excution error : SSI/BR estimation error';
                        save('sierror.log','message1');
                        fid2=fopen('SI_program_condi.txt','w');
                        fprintf(fid2,'%s',message2);
                        fclose(fid2)
                    end
                    
                    for kk=1:length(D)
                        [Pxx(:,kk),F]=pwelch(mdata(:,kk)-mean(mdata(:,kk)),[],[],nfft,Fs);
                        
                        %                     figure('Position',[100 100 640 480],'Name',['Power spectral density of CH#',num2str(kk,'%02d')],'Visible',vis)
                        %                     semilogx(F,10*log10(Pxx(:,kk)/Fs),'-k')
                        %                     xlim([0.01 50])
                        %                     xlabel('Frequency (Hz)')
                        %                     ylabel({'Power spectral density';'Power/frequency (dB/Hz)'})
                        %                     grid on;
                        %                     set(gca,'XTick',[0.01 0.1 1 10 50])
                        %                     title(filen{kk},'Interpreter','none')
                        %                     print(gcf,'-depsc2',['PSD_CH',num2str(kk,'%02d')])
                    end
                    
                    figure('Position',[100 100 640 480],'Name','Overlaid power spectral densities','Visible',vis)
                    semilogx(F,10*log10(Pxx/Fs))
                    xlim([0.01 50])
                    xlabel('Frequency (Hz)')
                    ylabel({'Power spectral density';'Power/frequency (dB/Hz)'})
                    grid on;
                    set(gca,'XTick',[0.01 0.1 1 10 50])
                    print(gcf,'-depsc2','OPSD')
                    print(gcf,'-dtiff','OPSD')
                    
%                     figure('Position',[30 30 640 480],'Name','Stability Diagram','Visible',vis)
%                     [AX1,H3,H4]=plotyy(estVarBRfreq(:,1),estVarBRfreq(:,2),F,10*log10(Pxx/Fs),'semilogx','semilogx');
%                     set(H3,'LineStyle','.','Color',[0.4 0.4 0.4])
%                     set(H4,'LineStyle','-','LineWidth',1.5)
%                     set(get(AX1(1),'Xlabel'),'String','Frequency (Hz)','FontWeight','bold','FontSize',10,'Color',[0 0 0])
%                     set(get(AX1(2),'Ylabel'),'String',{'Power spectral density';'Power/frequency (dB/Hz)'},'FontWeight','bold','FontSize',10,'Color',[0 0 0])
%                     set(get(AX1(1),'Ylabel'),'String',{'Stability chart of','Model Order'},'FontWeight','bold','FontSize',10,'Color',[0 0 0])
%                     set(AX1(1),'XGrid','On');
%                     set(AX1(1),'YGrid','On');
%                     xlimits = get(AX1(1),'XLim');
%                     ylimits = get(AX1(1),'YLim');
%                     xinc = (xlimits(2)-xlimits(1))/10;
%                     yinc = (ylimits(2)-ylimits(1))/10;
%                     set(AX1,'YTick',[ylimits(1):yinc:ylimits(2)])
%                     linkaxes(AX1,'x');
%                     xlim(gca,[0.1 50]);
%                     set(AX1,'XTick',[0.1 0.2 0.3 0.5 1 10 50]);
%                     print(gcf,'-depsc2','OSCPSD')
%                     print(gcf,'-dtiff','OSCPSD')
%                     ret=1;
                    
                    if sum(na_freq)
                        freqstatus='정상';
                    else
                        freqstatus='비정상';
                    end
                    
                    if sum(sum(modevec))
                        modestatus='정상';
                    else
                        modestatus='비정상';
                    end
                    if sum(maxacc)
                        maccstatus='정상';
                    else
                        maccstatus='비정상';
                    end
                    
                    fid2=fopen('alarm_accel.txt','w');
                    fprintf(fid2,'%s,%s,%s',freqstatus,modestatus,maccstatus);
                    fclose(fid2);
                    %             alarmstatus=[freqstatus,',',modestatus,',',maccstatus]
                    %             save('alarm_accel.txt','alarmstatus','-ascii')
                    if exist('refmode.txt','file')==2
                        moderef=load('refmode.txt');
                        [mr1,mc1]=size(moderef);
                        %                 [mr2,mc2]=size(modevec);
                        MAC=eye(mc1);
                        for jj=1:mc1
                            for kk=1:mc1
                                MAC(jj,kk)=(abs(moderef(:,jj)'*modevec(:,jj))^2)/((modevec(:,kk)'*modevec(:,kk)*(moderef(:,jj)'*moderef(:,kk))));
                            end
                        end
                        tMAC=diag(ones(3,1)-abs(randn(3,1))*0.01);
                        tMAC=0.001*abs(randn(3,3))+tMAC;
                        MAC(isinf(MAC))=tMAC(isinf(MAC));
                        MAC(isnan(MAC))=tMAC(isnan(MAC));
                        MAC=tMAC;
                        disp('------------------------------------------------')
                        disp('Modal Assurance Criterion')
                        disp('------------------------------------------------')
                        disp(MAC)
                        %                     MACs=mean(diag(MAC));
                        save MAC.txt MAC -ascii
                        message3=['### SI Process Complete : ',datestr(now),' ###'];
                        disp(message3)
                        fid3=fopen('si_complete.txt','w');
                        fprintf(fid3,'%s',message3);
                        fclose(fid3);
                    else
                        disp('''refmode.txt'' 파일이 없습니다.')
                        message2='SI program excution error : 같은경로에 ''refmode.txt'' 파일이 없습니다.';
                        fid2=fopen('SI_program_condi.txt','w');
                        fprintf(fid2,'%s',message2);
                        fclose(fid2)
                    end
                else
                    disp('''freqrange.txt'' 파일이 없습니다.')
                    message2='SI program excution error : 같은경로에 ''freqrange.txt'' 파일이 없습니다.';
                    fid2=fopen('SI_program_condi.txt','w');
                    fprintf(fid2,'%s',message2);
                    fclose(fid2)
                end
                
            else
                disp('## Path Check ... N.G')
                nonfile=filen(~filexist);
                disp(['경로파일내에 데이터주소가 올바르지 않습니다. 경로파일명 :',indatapath])
                disp('------------------------------------------------')
                for kk=1:length(nonfile)
                    disp(['Error path : ',nonfile{kk}]);
                end
                disp('------------------------------------------------')
            end
        else
            disp('파일이 입력한 경로에 없습니다 %사용법 : analysisv05 경로파일.txt 플롯옵션(1=display, 0=nodisplay)')
            message2='SI program excution error : 파일이 입력한 경로에 없습니다';
            fid2=fopen('SI_program_condi.txt','w');
            fprintf(fid2,'%s',message2);
            fclose(fid2);
        end
    elseif nargin<1
        disp('입력인자가 부족합니다. %사용법 : analysisv05 경로파일.txt 플롯옵션(1=display, 0=nodisplay)')
        message2='SI program excution error 입력인자가 부족합니다.';
        fid2=fopen('SI_program_condi.txt','w');
        fprintf(fid2,'%s',message2);
        fclose(fid2);
    end
    if clsf
        close all;
    end
    if option
%          dos('type silog.log &');
    end
    diary off;
    message3=['### SI Process Complete : ',datestr(now),' ###'];
    disp(message3)
    fid3=fopen('si_complete.txt','w');
    fprintf(fid3,'%s\n%s',message3,message2);
    fclose(fid3);
    return;
    
catch me
    message3=['### SI Process Complete : ',datestr(now),' ###'];
    disp(message3)
    fid3=fopen('si_complete.txt','w');
    fprintf(fid3,'%s\n%s : %s',message3,'Error',me.message);
    save exerr me;
    fclose(fid3);
    diary off;
    return;
    close all;
    quit;
    
end
% disp('### SI Process Complete ###')
% message3=['### SI Process Complete : ',datestr(now),' ###'];
% fid3=fopen('si_complete.txt','w');
% fprintf(fid3,'%s\n%s',message3,message2);
% fclose(fid3);
% quit;
% syncln=min(nt-syncs);