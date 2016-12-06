function [na_freq, modevec, xi, na_freq_cov, modevec_cov, xi_cov, estVarBRfreq, estVarBRzeta] = sfadccore(data,Fs,minorder,maxorder,filename,freqband)

% lambda=0;
% phi=0;
% xi=0;
fwr=0;

if nargin==4
    fwr=0;
end

ecvar=[];
Nsensor=size(data,2);
GD.dt=1/Fs;
FD.fcut=Fs/2;
[GD.Ndata,GD.Nsensor]=size(data);
disp(['## Calculating the Correlation Function Matrices...']);

TD.NdataCorr=256;
TD.Max_order=maxorder;
for isensor=1:GD.Nsensor
    for jsensor=1:GD.Nsensor
        Rtemp=xcorr(data(:,isensor),data(:,jsensor),TD.NdataCorr);
        Loc(isensor,jsensor).Ryy(1:TD.NdataCorr)=Rtemp(TD.NdataCorr:2*TD.NdataCorr-1);
        clear Rtemp
    end
    disp([num2str(isensor/GD.Nsensor*100,3),'% processed'])
end

disp(['## Constructing the Generalized Hankel Matrices for SSI/BR...'])
TD.n1=128;
TD.n2=128;

for isensor=1:GD.Nsensor
    for jsensor=1:GD.Nsensor
        for i=1:max(TD.n1,TD.n2)
            R(isensor,jsensor,:)=Loc(isensor,jsensor).Ryy;
        end
    end
end

H0.data=[];
H1.data=[];
for i=1:TD.n1
    temp1=[];
    temp2=[];
    for j=1:TD.n2
        temp1=[temp1 R(:,:,i+j-1)];
    end
    H0.data=[H0.data;temp1];
end
disp(['## Singular Value Decomposition'])
[H0.u, H0.s, H0.v]=svds(H0.data,TD.Max_order);

%%
disp(['## Building Modal Stability Chart...'])
tempVar=[];

Min_order=minorder;
Max_order=maxorder;

for Norder=Min_order:Max_order
    tempVar.U1=H0.u(1:TD.n1*GD.Nsensor,1:Norder);
    tempVar.S1=H0.s(1:Norder,1:Norder);
    tempVar.V1=H0.v(1:TD.n2*GD.Nsensor,1:Norder);
    tempVar.Op=tempVar.U1*sqrt(tempVar.S1);
    
    tempVar.Op1=tempVar.Op(GD.Nsensor+1:TD.n1*GD.Nsensor,:);
    tempVar.Op2=tempVar.Op(1:(TD.n1-1)*GD.Nsensor,:);
    
    tempVar.A=pinv(tempVar.Op2)*tempVar.Op1;
    tempVar.C=tempVar.Op(1:GD.Nsensor,:);
    
    [tempVar.psi,tempVar.z]=eig(tempVar.A);
    tempVar.phi=tempVar.C*tempVar.psi;
    
    %   modal frequency and modal damping
    
    tempVar.icheck=[];
    for imode=1:Norder,
        if (imag(tempVar.z(imode,imode))<0) tempVar.icheck=[tempVar.icheck;imode];,end
    end
    
    if (isempty(tempVar.icheck))
        estVarBR{Norder}.freq=[];
    end
    
    if (~isempty(tempVar.icheck))
        tempVar.Nmodes=length(tempVar.icheck);
        tempVar.temp=diag(tempVar.z);
        tempVar.lamda=log(tempVar.temp(tempVar.icheck))/GD.dt;
        
        tempVar.wr=imag(tempVar.lamda);
        tempVar.sr=real(tempVar.lamda);
        tempVar.e_zeta =-tempVar.sr./sqrt(tempVar.wr.^2+tempVar.sr.^2);
        tempVar.e_omega=-tempVar.wr./sqrt(1-tempVar.e_zeta.^2);
        tempVar.phi=tempVar.phi(:,tempVar.icheck);
        
        [tempVar.e_omega,tempVar.iorder]=sort(tempVar.e_omega);
        tempVar.e_zeta=tempVar.e_zeta(  tempVar.iorder);
        tempVar.phi   =tempVar.phi   (:,tempVar.iorder);
        for imode=1:tempVar.Nmodes
            [maxAbs,maxLoc]=max(abs(tempVar.phi(:,imode)));
            tempVar.phi2(:,imode)=real(tempVar.phi(:,imode)*exp(-phase(tempVar.phi(maxLoc,imode))*sqrt(-1)));
        end
        estVarBR{Norder}.freq=tempVar.e_omega/2/pi;
        %         ecvar=[ecvar;Norder,tempVar.e_omega/2/pi]
        estVarBR{Norder}.zeta=tempVar.e_zeta;
        estVarBR{Norder}.mode=[];
        for imode=1:tempVar.Nmodes
            [maxAbs,maxLoc]=max(abs(tempVar.phi(:,imode)));
            estVarBR{Norder}.mode(:,imode)=real(tempVar.phi(:,imode)*exp(-phase(tempVar.phi(maxLoc,imode))*sqrt(-1)));
        end
    end
end

%   critical values for determining the stability of each poles

critical_frequency_difference=0.1;
critical_zeta_ratio=1;
critical_damping=1;
critical_MAC_value=0.85;

for Norder=Min_order:Max_order,
    estVarBR{Norder}.stability=[];
end

%
%   Criteria for noise modes
%

for Norder=1:Max_order,
    if (~isempty(estVarBR{Norder}.freq)) Min_order=Norder; break; end
end


for Norder=Min_order:Max_order,
    for i=1:length(estVarBR{Norder}.freq)
        if (estVarBR{Norder}.zeta(i)<=0 | estVarBR{Norder}.zeta(i)>=critical_damping)
            estVarBR{Norder}.stability(i)=1;   % stability=1: noise mode
        else
            estVarBR{Norder}.stability(i)=0;   % stability=0: unclassified, now.
        end
    end
end

%
%   Criteria for stable and unstable modes
%

for Norder=Min_order+1:Max_order,
    for i=1:length(estVarBR{Norder}.freq)
        for j=1:length(estVarBR{Norder-1}.freq)
            if ((estVarBR{Norder-1}.stability(j)~=1) && (estVarBR{Norder}.stability(i)~=1))
                frequency_difference=abs(estVarBR{Norder-1}.freq(j)-estVarBR{Norder}.freq(i));
                if (frequency_difference<critical_frequency_difference)
                    zeta_ratio=(estVarBR{Norder-1}.zeta(j)-estVarBR{Norder}.zeta(i))/estVarBR{Norder}.zeta(i);
                    MAC_value=((estVarBR{Norder-1}.mode(:,j)'*estVarBR{Norder}.mode(:,i))/(norm(estVarBR{Norder-1}.mode(:,j))*norm(estVarBR{Norder}.mode(:,i))))^2;
                    if (zeta_ratio<critical_zeta_ratio && MAC_value>critical_MAC_value)
                        estVarBR{Norder-1}.stability(j)=2;
                        estVarBR{Norder}  .stability(i)=2;
                    end
                end
            end
        end
    end
end

% eval(['icheck=exist(''',infile.root,'.schart'',''file'');'])
%
% if (icheck==2)
%     eval(['save ',infile.root,'.schart estVarBR Min_order Max_order -mat -append'])
% else
%     eval(['save ',infile.root,'.schart estVarBR Min_order Max_order -mat'])
% end
% getmBR

%%
disp(['## Stabilization chart has constructed'])
if fwr
    fid1=fopen(filename,'wt');
end
kk=0;estVarBRfreq=[];estVarBRzeta=[];
for Norder=1:maxorder
    if isfield(estVarBR{Norder},{'mode'})
        [rs,cs]=size(estVarBR{Norder}.mode);
        if fwr
            wstring=['order#',num2str(Norder),',',num2str(rs),',',num2str(cs),',freq,',num2str(estVarBR{Norder}.freq','%.15e'),',zeta,',num2str(estVarBR{Norder}.zeta','%.15e')];
            fprintf(fid1,'order#%i,%i,%i,freq,',Norder,rs,cs);
        end
        for jj=1:cs
            if fwr
                fprintf(fid1,'%.15e,',estVarBR{Norder}.freq(jj));
            end
            kk=kk+1;
            estVarBRfreq(kk,:)=[estVarBR{Norder}.freq(jj),Norder];
            estVarBRzeta(kk,:)=[estVarBR{Norder}.zeta(jj),Norder];
        end
        if fwr
            fprintf(fid1,'zeta,');
            for jj=1:cs
                fprintf(fid1,'%.15e,',estVarBR{Norder}.zeta(jj));
            end
            fprintf(fid1,'mode,');
            for ii=1:rs
                for jj=1:cs
                    fprintf(fid1,'%.15e,',estVarBR{Norder}.mode(ii,jj));
                end
            end
            fprintf(fid1,'stability,');
            for jj=1:cs
                fprintf(fid1,'%.15e,',estVarBR{Norder}.stability(jj));
            end
            fprintf(fid1,'\n');
        end
    else
        if fwr
            [rs,cs]=size(estVarBR{Norder}.freq);
            wstring=['order#',num2str(Norder),',',num2str(rs),',',num2str(cs),',freq,',num2str(estVarBR{Norder}.freq,'%.15e'),',zeta,','0'];
            fprintf(fid1,'order#%i,%i,%i,freq,0,zeta,0,mode,0,stability,0\n',Norder,rs,cs);
        end
    end
    if fwr
        disp(wstring)
    end
    %     fprintf(fid1,'%s\n',wstring);
end
if fwr
    fclose(fid1);
end

ifreq=size(freqband,1);
f=[];
z=[];
m=[];
for kk=1:ifreq
    for Norder=minorder:maxorder
        for i=1:length(estVarBR{Norder}.freq)
            if (estVarBR{Norder}.stability(i)>=1 && estVarBR{Norder}.freq(i)>freqband(kk,1) && estVarBR{Norder}.freq(i)<freqband(kk,2))
                f=[f;estVarBR{Norder}.freq(i)];
                z=[z;estVarBR{Norder}.zeta(i)];
                m=[m;estVarBR{Norder}.mode(:,i)'];
%                 plot(estVarBR{Norder}.freq(i),Norder,'r.')
            end
        end
    end

    N=length(f);

    if (N>1)
        for i=1:N
            m(i,:)=m(i,:)/norm(m(i,:));
        end

        mac=[];
        for i=1:N
            for j=1:N
                mac(i,j)=(m(i,:)*m(j,:)')/norm(m(i,:))/norm(m(j,:));
            end
        end

        for i=2:N
            if (mac(1,i)<0) m(i,:)=-m(i,:);
            end
        end

        mac=abs(mac);
        
        [i,j,s]=find(mac>0.999);
        S=sparse(i,j,s,N,N);
        M=full(S);
        nz=[];
        for i=1:N
            nz(i)=nnz(M(i,:));
        end

        [temp,iord]=sort(-nz);
        mac=[];
        for i=1:N
            for j=1:N
                mac(iord(i),iord(j))=(m(iord(i),:)*m(iord(j),:)')/norm(m(iord(i),:))/norm(m(iord(j),:));
            end
        end

        for i=2:N
            if (mac(iord(1),iord(i))<0)
                m(iord(i),:)=-m(iord(i),:);
                mac(iord(1),iord(i))=-mac(iord(1),iord(i));
            end
        end
%         save mac mac
        cov_f=[];
        cov_z=[];
        for i=1+1:N
            cov_f(i)=std(f(iord(1:i)))/mean(f(iord(1:i)));
            cov_z(i)=std(z(iord(1:i)))/mean(z(iord(1:i)));
        end

        icheck=3;
        size(iord)
        %     while (icheck<0 | icheck>N)
        %         icheck=input('Number of considering modes:');
        %     end
%         size(f)
        freq=mean(f(iord(1:icheck)));
        zeta=mean(z(iord(1:icheck)));
        mode=mean(m(iord(1:icheck),:));
        freq_cov=std(f(iord(1:icheck)))/freq;
        zeta_cov=std(z(iord(1:icheck)))/zeta;
        mode_cov=std(m(iord(1:icheck),:))./mode;
%         mode
        na_freq(kk)=freq;
        modevec(:,kk)=mode';
        xi(kk)=zeta;
        na_freq_cov(kk)=freq_cov;
        modevec_cov(:,kk)=mode_cov';
        xi_cov(kk)=zeta_cov;
        disp([' - frequency    : f=',num2str(freq,5)    ,'  f_std=',num2str(freq_cov*100,5)])
        disp([' - damping ratio: z=',num2str(zeta*100,5),'  z_std=',num2str(zeta_cov*100,5)])
    else
        na_freq(kk)=0;
        modevec(:,kk)=zeros(Nsensor,1);
        xi(kk)=0;
        na_freq_cov(kk)=0;
        modevec_cov(:,kk)=zeros(Nsensor,1);
        xi_cov(kk)=0;
    end
f=[];
z=[];
m=[];
end