% changing the STO thickness


clc;
clf;
clear all;
hold all;

c=3e8;

count=1;
w=55:0.005:90; % step size for angles (X - axis)
% try and change angles if FOM calc me error
kkk_range = 1:1:2; % 1.33,1.34 - analyte RI varies above water
iii_range = w; % jitni baar w hoga chalega
total_iterations = length(iii_range);
td_values = [2500e-9]; % Array of different td values ; thickness of analyte(2500nm)
IR12_nd_133 = [];

for td_index = 1:length(td_values) %
   
    td = td_values(td_index); % Current td value
    figure; % Creating a new figure for the current td value
    title(sprintf('td = %e', td));
    hold all;

    for kkk=kkk_range
        for iii=iii_range
            lambda0=1.550e-6; %1550nm wavelength irradiation
            lamdac=2.4511e-5; %collision wavelength
            lamdap=1.0657e-7; %plasmon wavelength
            %below 2 lines are expression for dielectric (drude model ig)
            epsilonreal=1-((lambda0.^2.*(lamdac)^2)./(lamdap^2.*(lambda0.^2+lamdac^2)));
            epsilonim=((lambda0.^3.*(lamdac))./(lamdap^2.*(lambda0.^2+lamdac^2)));
            % sqrt of dielectric= RI
            n2=sqrt((sqrt(epsilonreal.^2+epsilonim.^2)+epsilonreal)./2); %
            k2=sqrt((sqrt(epsilonreal.^2+epsilonim.^2)-epsilonreal)./2);
   
            Numords=1;  % number of diffractive orders maintained, always odd
                            % order of diff [-50,50]

%             region 1 cover refractive index
%             nc=1.5007;  %BK7
            %nc=1.426;   %CaF2
%             nc=1.642    %SF5
%             nc=1.69     %SF10
%             nc=1.743    %SF11
              % nc=1.3245;  %NaF Sodium Floride  at 633nm
              nc=1.3194; %NaF at 1550nm
         
            % different Substrate:Caf2,BK7,SF10
           
           
            ns=1;  % region 3 substrate refractive index (air)
            Ngrat=4;   % number of grating slices; (no of layers)
            period=100e-9;  % grating period in microns
   
            nd=1.33+(kkk-1)*.01; %varies the RI of analyte for current iteration
            %nd = 1.33;

            nm = n2-1i*k2; % basically real +im RI
            A = 3.44904; % line no 67 -71 are constants nsi formula
            A1 = 2271.88813;
            A2 = 3.39538;
            t1 = 0.058304;
            t2 = 0.30384;
            nsi = A+A1*exp(-1.55/t1)+A2*exp(-1.55/t2); %Si ka RI
           
            %below are thickness and RI of Alo and eth
            neth = 1.36;
            nALO = 1.746; %both not being used rn
            teth = 0e-9;
            tALO = kkk*5e-9; %both not being used rn
           
           nSTO = 2.284;
           
            %Second metal
            nm1=0.57470-1i*9.6643; %RI of Gold at 1550
           
             %DielectricMetal -MgF2
             D=10;  %no of layers of MgF2
             nd1=1.3705;
             
            %MoS2
            T=1 ;   % layers(single layer MoS2 0.65nm)
            nT=3.647 ;  %%%% RI of MoS2
           
            % BP
            BP=1;   % layers(single layer BP 0.53nm)
            nbp=3.0945-1i*1.58369;   %%% RI of  BP (black phosphorous)
               
            ngr=2.974-2.803*1i;

            td = td_values(td_index); %baadme

           
           % IMP positioning of layers starts from here
           
            %%%%    [Al-30nm,Au-5nm,MgF2-MoS2,BP,Analyte-2500]
            % depth = [30e-9,0e-9,D*0e-9,T*0e-9,BP*0e-9,td];  %
           
            %ridge is bulky part; groove is the bottom part of
            %nanostructure

            depth=[30e-9,6e-9,0.34e-9,td];
            nr = [nm,nSTO,ngr,nd];  % Ridge refractive index for each grating
            ng = [nm,nSTO,ngr,nd];  % Groove refractive index for each grating
            %ngrat = 4  excludes air and glass; ngrat is total layers

            Filfac = [.5 .5 .5 .5];  % fill factor for ridges comes when dealing eith nanostruc
                                           % basically equal width of
                                           % ridges and grooves

            Disp = [0 0 0 0];   % ridge displacement in a frac  when layer non uniform                                                                                                                                                                                                                                            tion of period
   
            %%%% rest is Common
           
            theta0=iii;  %   angle of incidence for this iteration
            phi0=0;    % azimuthal angle of incidence
            deg=pi/180;
            Nmax=(Numords-1)/2; % highest order number retained (50
            I=(-Nmax:Nmax)';  % I is the order indexes
            p=(Numords+1)/2;  % index of zeroth order
            theta0=theta0*deg;
            phi0=phi0*deg;  % converting in radians
            epsc=nc^2;    % relative permittivity
            epss=ns^2;
   
            k0=2*pi/lambda0;  % free space vector
            K=2*pi/period;    % grating vector
   
            kc=k0*nc;
            kx=kc*sin(theta0)*cos(phi0)-I*K;  % region1 wave vector components
   
            ky=kc*sin(theta0)*sin(phi0)*ones(size(kx));
            kzc=sqrt(kc^2-kx.^2-ky.^2);
            bad_indices=find((real(kzc)-imag(kzc))<0);  
            kzc(bad_indices)=-kzc(bad_indices);              
            ks=k0*ns;
            kzs=sqrt(ks^2-kx.^2-ky.^2);
            bad_indices=find((real(kzs)-imag(kzs))<0);  % region3 wavevector
            kzs(bad_indices)=-kzs(bad_indices);
   
            %%%%%% define some  matrices and vectors %%%%%%%
            Zv=zeros(Numords,1);
            Dv=Zv;
            Dv(p)=1;
            Zv2=[Zv;Zv];
            Eye=eye(Numords);  % identity matrix
            Kx=diag(kx/k0);
            Kxsq=Kx.^2;
            Kzc=diag(kzc/k0);
            Kzcsq= Kzc.^2;
            Kzs=diag(kzs/k0);
            Kzssq= Kzs.^2;
            M=Numords-1;
            temp1=Eye/ns;
            fmat=Eye;
            gmat=j*Kzs/ns^2;
   
            for ii=Ngrat:-1:1
                epsg=ng(ii).^2;  % groove permittivity
                epsr=nr(ii).^2;   % ridge permittivity
                epsG=(1-Filfac(ii))*epsg+Filfac(ii)*epsr;  % average grating
                iepsG=(1-Filfac(ii))/epsg+Filfac(ii)/epsr;
                Sinc=sin(pi*Filfac(ii)*(1:M))./(pi*(1:M));
                vm=(epsr-epsg)*fliplr(Sinc);
                v0=epsG;
                vp=(epsr-epsg)*Sinc;
                v=[vm v0 vp].*exp(+j*2*pi*Disp(ii)*(-M:M));
                ivm=(1/epsr-1/epsg)*fliplr(Sinc);
                iv0=iepsG;
                ivp=(1/epsr-1/epsg)*Sinc;
                iv=[ivm iv0 ivp].*exp(+j*2*pi*Disp(ii)*(-M:M));
                Epsilon=toeplitz(fliplr(v(1:Numords)),v(Numords:2*Numords-1));
                Alpha=toeplitz(fliplr(iv(1:Numords)),iv(Numords:2*Numords-1));
   
                clear Sinc v  vm  v0  vp
                B=Kx*(Epsilon\Kx)-Eye;  % cofficient matrix
                [W,V]=eig(Alpha\B);    % W is the eigen vector and V are the eigen values
                Q=sqrt(V);
                M0=Alpha*W*Q;
                E=expm(-k0*Q*depth(ii));
                v=[W W;M0,-M0]\[fmat;gmat];
                temp2=v(1:Numords,:)\E;
                temp3=E*v(Numords+1:2*Numords,:)*temp2;
                temp1=temp1*temp2;
                fmat=W+W*temp3;
                gmat=M0-M0*temp3;
   
            end
   
            gfi=gmat/fmat;
            RHS=-gfi(:,p);
            RHS(p)=RHS(p)+j*kzc(p)/k0/epsc;
            Rs=(gfi+j*Kzc/nc^2)\RHS;
            Ts=(temp1/fmat)*(Rs+Dv)*nc;
         
   
           IR1=(abs(Rs).^2).*real(kzc./kzc(p));
           IT1=(abs(Ts).^2).*real(kzs./kzc(p));
   
           e=sum(IT1);
           f=sum(IR1);
           g=1-e-f;
           IT12(count)=e;
           IR12(count)=f;
           loss(count)=g;

           if nd == 1.33
               IR12_nd_133 = IR12;
           end
 
           count=count+1;
     
           clc;
           progress = (count/total_iterations)*100;
           %fprintf('td = %.2e\t', td);
           
           % fprintf('RI of BP = %d\t',nbp);
           % fprintf('No. of BP Layers(.53)=%d\n',BP);
           % fprintf('RI of MoS2 = %d\t',nT);
           % fprintf('No. of MoS2 Layers(0.65nm) =%d\n',T);
           
           fprintf('nd = %.2f\t', nd);
           fprintf('Progress: [%s%s] %.2f%%\r', repmat('=',1,floor(progress/2)), repmat(' ',1,50-floor(progress/2)), floor(progress));
   
        end
   
        %%%%% Plotting results %%%%%
       
        % Finding the minimum y-value and its corresponding x-value
        [minVal, minIdx] = min(IR12);
        minW = w(minIdx);
       
        % Storing min values and respective x values for each kkk
        minVals(kkk) = minVal;
        minWs(kkk) = minW;
       
        % Plotting the curve for current nd val
        plot(w, IR12, 'DisplayName', ['nd = ', num2str(nd)]);
        box ON;
       
        % Resetting count for next iteration
        count=1;
       
        % Storing IR12 values for sensitivity calculation in next step
        IR12_store{kkk} = IR12;
       
    end
   
    if nd == 1.33+(kkk-1)*.01
        % Calculating sensitivity
        sensitivity = abs(minWs(1)-minWs(2)) / 0.01;
   
         % Finding all the local minima (dips) in the IR12 curve
        [mins, minIndices] = findpeaks(-IR12);
       
        % If there are multiple dips, find the one with the largest magnitude (lowest y value)
%         if length(mins) > 1
%             [minValue1, minLoc] = min(-mins);
%             minIndex1 = minIndices(minLoc);
%         else
%             [minValue1, minIndex1] = min(IR12);
%         end
   
        % Calculating the FWHM for the SPR dip
        % half_max_value = (max(IR12) - min(IR12))/2;
       
        % Finding the main dip in the IR12 curve for nd=1.33
        [minValue1, minIndex1] = min(IR12_nd_133);
       
        % Calculating the FWHM for the SPR dip
        half_max_value = (max(IR12_nd_133) + min(IR12_nd_133))/2;
       
        left_cross_index = find(IR12_nd_133(1:minIndex1) >= half_max_value, 1, 'last');
        right_cross_index = find(IR12_nd_133(minIndex1:end) >= half_max_value, 1, 'first') + minIndex1 - 1;
        FWHM = w(right_cross_index) - w(left_cross_index);
       
        % Debug statements
disp(['Sensitivity: ', num2str(sensitivity)]);
disp(['FWHM: ', num2str(FWHM)]);

        % Calculating FOM (Sensitivity / FWHM)
        FOM = sensitivity / FWHM;
   
        scatter(w(left_cross_index), IR12_nd_133(left_cross_index), 'rx'); % Marking the left cross point with a red 'x'
        scatter(w(right_cross_index), IR12_nd_133(right_cross_index), 'bx'); % Marking the right cross point with a blue 'x'
       
         Rmin=minValue1;
         
        % Displaying sensitivity on the plot
        dim = [.7 .5 .3 .3];
        str = {['Sensitivity: ', num2str(sensitivity)], ['FWHM: ', num2str(FWHM)], ['FOM: ', num2str(FOM)],['Rmin: ', num2str(Rmin)]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','white');
    end

    % Setting plot title and labels
    %title(['td: ', num2str(td), '--- No of MoS2 layers:',num2str(M_L),'No.of  TiSi2 layers:',num2str(T)])
    title('NaF-Al-Au-MgF2(1000)-MoS2-BP-Analyte- 1.33-1.34')
    xlabel('Angle');
    ylabel('Reflectivity');
    legend;
   
    % Resetting minVals and minWs for next td value
    minVals = [];
    minWs = [];
end
   
%     format long ;
%     Rmin