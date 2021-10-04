function [V] = layer_SH_analysis(nmax,geoid,Re,rhoE,max_bin,fupper,flower,fdens,falpha)
%
% This function converts a variable density layer into SH-coefficients. The
% following input variables are used:
%
%   - nmax:             maximum degree SH-coefficient (nmax=length(lat)-1)
%   - geoid:            string describing geoid used
%                           - 'none'    (spherical geoid)
%                           - 'WGS84'
%   - Re:               value of radius [m]
%   - rhoE:             value of mean density [kg/m3] (linking to specific GM)
%   - max_bin:          maximum binomial terms that can be used (now max=3)
%
%%%%%%% Layer input values (format: matrix m*n [m=lat,n=lon (block grid definition)]):
%
%   - fupper:           value for upper boundary varying from Re [m]
%   - flower:           value for lower boundary varying from Re [m]
%   - fdens:            value for density [kg/m3]
%   - falpha:           value for linear density gradient [1/m]
%                       default=0
%
%%%%%% output values:
%
%   - V:                Spherical harmonic coefficients (4pi-normalization)
%     format:           [l,m,Clm,Slm]     l = 0,1,2,3,... m = 0,0,0,0,... etc.
%
% Made by: Bart Root, TUDelft, Astrodynamics and Space Missions
% First version 1.0: date: 08-04-2013
%
% 06-feb-2014: Put in the anti-aliasing part, but don't use this with high
% degree of SH-coefficients.
%
% Matlab functions used:
%
% - GSHA
% - cs2sc
% - sc2vecml
% - geocradius
% - size
% - zeros
%
%%%%%%%%%%%%%%% Start function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initializing certain values

grid = 'block';
method = 'ls';
rhoRatio = 1./rhoE;

% maximum binomial terms check (for now max = 3)

if max_bin>8
    error('Maximum number of binomial terms is overriden')
end

% No alpha stated, don't do extra the GSHA loops
if nargin<9
    radOn = 0;
    falpha = 0;
else
    radOn = 1;
end

% Do all the matrices have the same format
if size(fupper)~=size(flower)
    error('Boundary matrices are not the same size')
elseif size(fdens)~=size(fupper)
    error('Density matrix is not the same size as the boundary matrices')
elseif ( radOn==1 )
    if size(fdens)~=size(falpha)
        error('Linear density gradient matrix is not the same size as the boundary matrices')
    end
else
    % all is correct: do nothing
end

% Check if nmax value is correct
if nmax+1>size(fupper,1)
    error('Size input matrix is incorrect: Latitude length is to small for nmax')
elseif nmax+1*2>size(fupper,2)
    error('Size input matrix is incorrect: Longitude length is to small for nmax')
end

% Geoidetic radius reference correction
if strcmp(geoid,'none')
    % do nothing: spherical case
elseif strcmp(geoid,'WGS84')        
        res = 180/size(fupper,1);        
        latV = ((-90 +res/2:res:90-res/2).*-1)';
        lat = repmat(latV,1,size(fupper,2));
        
        correctionR = geocradius(lat,'WGS84') - Re;
        fupper = fupper + correctionR;
        correctionR = geocradius(lat,'WGS84') - Re;
        flower = flower + correctionR;
        clear lat
        
elseif strcmp(geoid,'benchmark_Mikhail')
        % not needed anymore!
        disp('Benchmark Mikhail reference')       
        
        f  = 1/298.257223563;
        aRe = 6371000;
        
        res = 180/size(fupper,1);        
        latV = ((-90 +res/2:res:90-res/2).*-1)';
        lat = repmat(latV,1,size(fupper,2));
        
        correctionR = geocradius(lat,f,aRe) - Re;
        fupper = fupper + correctionR;
        correctionR = geocradius(lat,f,aRe) - Re;
        flower = flower + correctionR;        
        clear lat 
        
else
    % An incorrect geoid is chosen
    error('An incorrect geoid reference is chosen')
end

%%%%                                                       %%%%
%%%% Starting with the binomial terms loop in the Analysis %%%%
%%%%                                                       %%%%

for Hi = 1:max_bin
    
    %disp(['binary number: ' num2str(Hi)])

    % calibrate density if alpha is used
%     if radOn ==1
%         fdens = fdens./(1+falpha .* ((fupper+flower)./2));
%     end

    % K = Hi
    switch Hi
        case 1
            m = (fdens).*(fupper-flower);
        case 2
            m = (fdens).*(fupper.*fupper-flower.*flower);
        case 3
            m = (fdens).*(fupper.*fupper.*fupper-flower.*flower.*flower);
        case 4
            m = (fdens).*(fupper.*fupper.*fupper.*fupper-flower.*flower.*flower.*flower);
        case 5
            m = (fdens).*(fupper.*fupper.*fupper.*fupper.*fupper-flower.*flower.*flower.*flower.*flower);
        case 6
            m = (fdens).*(fupper.*fupper.*fupper.*fupper.*fupper.*fupper-flower.*flower.*flower.*flower.*flower.*flower);
        case 7
            m = (fdens).*(fupper.*fupper.*fupper.*fupper.*fupper.*fupper.*fupper-flower.*flower.*flower.*flower.*flower.*flower.*flower);
        case 8
            m = (fdens).*(fupper.*fupper.*fupper.*fupper.*fupper.*fupper.*fupper.*fupper-flower.*flower.*flower.*flower.*flower.*flower.*flower.*flower);
        otherwise
            error('Incorrect number of binomial terms')
    end
    
%     %%%%%%% test anti-aliasing
%     m_a = zeros(size(m)*4);
% 
%     for ii = 1:size(m,1)
%         for jj = 1:size(m,2)
% 
%             i = 1 + (ii-1) * 4;
%             j = 1 + (jj-1) * 4;
% 
%             D = m(ii,jj)*ones(4,4);
% 
%             m_a(i:i+(4-1),j:j+(4-1)) = D;
%         end
%     end
    
    
    % Do the analysis
    cs = GSHA(m,method,grid,nmax); sc = cs2sc(cs);
    %figure;plot(1:size(sc,1),log10(sqrt(sum(sc.^2,2))./(1:size(sc,1))'));
    [Clm,Slm,llvec,mmvec] = sc2vecml(sc,nmax);
    
    %%%%%%%%%%%%%%%%%%
    Vpart = [llvec' mmvec' Clm Slm];

    % Stating the factors
    fac =  (3./(((2.*Vpart(:,1))+1)));
    fac2 = (((Vpart(:,1)+2)/2));
    fac3 = (((Vpart(:,1)+2).*(Vpart(:,1)+1)/6));
    fac4 = (((Vpart(:,1)+2).*(Vpart(:,1)+1).*(Vpart(:,1))/24));
    fac5 = (((Vpart(:,1)+2).*(Vpart(:,1)+1).*(Vpart(:,1)).*(Vpart(:,1)-1)/120));
    fac6 = (((Vpart(:,1)+2).*(Vpart(:,1)+1).*(Vpart(:,1)).*(Vpart(:,1)-1).*(Vpart(:,1)-2)/720));
    fac7 = (((Vpart(:,1)+2).*(Vpart(:,1)+1).*(Vpart(:,1)).*(Vpart(:,1)-1).*(Vpart(:,1)-2).*(Vpart(:,1)-3)/5040));
    fac8 = (((Vpart(:,1)+2).*(Vpart(:,1)+1).*(Vpart(:,1)).*(Vpart(:,1)-1).*(Vpart(:,1)-2).*(Vpart(:,1)-3).*(Vpart(:,1)-4)/40320));
    
    % If needed also analyse the linear density gradient term
    if radOn ==1   
        if Hi==2
            % Do the SHA on the radial dens component.
            % K = 2
            m = (falpha.*Re).*(fdens).*(fupper.*fupper-flower.*flower);
%             %%%%%%% test anti-aliasing
%             m_a = zeros(size(m)*4);
% 
%             for ii = 1:size(m,1)
%                 for jj = 1:size(m,2)
% 
%                     i = 1 + (ii-1) * 4;
%                     j = 1 + (jj-1) * 4;
% 
%                     D = m(ii,jj)*ones(4,4);
% 
%                     m_a(i:i+(4-1),j:j+(4-1)) = D;
%                 end
%             end
            
            cs = GSHA(m,'wls',grid,nmax); sc = cs2sc(cs);
            [Clm,Slm,llvec,mmvec] = sc2vecml(sc,nmax);
            
            %%%%%%%%%%%%%%%%%% 
            Vpart_rad = [llvec' mmvec' Clm Slm];
            
            % Adding the radial component to the other coefficients
            %Vpart(:,3:4) = Vpart(:,3:4) + Vpart_rad(:,3:4);
        elseif Hi==3
            % K = 3
            m = (2.*falpha.*Re).*(fdens).*(fupper.*fupper.*fupper-flower.*flower.*flower);
%             %%%%%%% test anti-aliasing
%             m_a = zeros(size(m)*4);
% 
%             for ii = 1:size(m,1)
%                 for jj = 1:size(m,2)
% 
%                     i = 1 + (ii-1) * 4;
%                     j = 1 + (jj-1) * 4;
% 
%                     D = m(ii,jj)*ones(4,4);
% 
%                     m_a(i:i+(4-1),j:j+(4-1)) = D;
%                 end
%             end
            
            cs = GSHA(m,'wls',grid,nmax); sc = cs2sc(cs);
            [Clm,Slm,llvec,mmvec] = sc2vecml(sc,nmax);
            
            %%%%%%%%%%%%%%%%%%
            Vpart_rad = [llvec' mmvec' (Vpart(:,1)+2).*Clm (Vpart(:,1)+2).*Slm];
            
            % Adding the radial component to the other coefficients
            %Vpart(:,3:4) = Vpart(:,3:4) + Vpart_rad(:,3:4);
        end 
        
    end
    
    % Constructing the correct coefficients
    if Hi==1
        V = zeros(size(Vpart));
        V(:,1) = Vpart(:,1);
        V(:,2) = Vpart(:,2);
    end

    switch Hi
        case 1
            VC = (Vpart(:,3)).*fac.*rhoRatio./Re;
            VS = (Vpart(:,4)).*fac.*rhoRatio./Re;
        case 2
            VC = fac2.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re;
            VS = fac2.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re;
            
            % with lineair density variation
            if radOn ==1
                VC = VC + (Vpart_rad(:,3)).*fac.*rhoRatio./Re./Re/2;
                VS = VS + (Vpart_rad(:,4)).*fac.*rhoRatio./Re./Re/2;
            end
            
        case 3
            VC = fac3.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re;
            VS = fac3.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re;
            
            % with lineair density variation
            if radOn ==1
                VC = VC + (Vpart_rad(:,3)).*fac.*rhoRatio./Re./Re./Re/6;
                VS = VS + (Vpart_rad(:,4)).*fac.*rhoRatio./Re./Re./Re/6;
            end
            
        case 4
            VC = fac4.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re./Re;
            VS = fac4.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re./Re;
        case 5
            VC = fac5.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re./Re./Re;
            VS = fac5.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re./Re./Re;
        case 6
            VC = fac6.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re;
            VS = fac6.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re;
        case 7
            VC = fac7.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re./Re;
            VS = fac7.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re./Re;
        case 8
            VC = fac8.*(Vpart(:,3)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re./Re./Re;
            VS = fac8.*(Vpart(:,4)).*fac.*rhoRatio./Re./Re./Re./Re./Re./Re./Re./Re;
        otherwise
            error('Incorrect number of binomial terms')
    end
    
    % Summate the coefficients
    switch Hi
        case 1
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 2
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 3
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 4
            VC(Vpart(:,1)==0) = 0;
            VS(Vpart(:,1)==0) = 0;
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 5
            VC(Vpart(:,1)==0|Vpart(:,1)==1) = 0;
            VS(Vpart(:,1)==0|Vpart(:,1)==1) = 0;
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 6
            VC(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2) = 0;
            VS(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2) = 0;
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 7
            VC(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2|Vpart(:,1)==3) = 0;
            VS(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2|Vpart(:,1)==3) = 0;
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        case 8
            VC(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2|Vpart(:,1)==3|Vpart(:,1)==4) = 0;
            VS(Vpart(:,1)==0|Vpart(:,1)==1|Vpart(:,1)==2|Vpart(:,1)==3|Vpart(:,1)==4) = 0;
            V(:,3) = V(:,3) + VC;
            V(:,4) = V(:,4) + VS;
        otherwise
            error('Incorrect number of binomial terms')
    end

%     figure
%     plot(log10(abs(VC(1:180))))
%     hold on
%     %plot(log10(abs(VS(1:180,4))),'r')
%     hold off
end

%     figure
%     plot(log10(fac(1:180).*rhoRatio./Re))
%     hold on
%     plot(log10(fac2(1:180).*fac(1:180).*rhoRatio./Re./Re),'r')
%     plot(log10(fac3(1:180).*fac(1:180).*rhoRatio./Re./Re./Re),'g')
%     plot(log10(fac4(1:180).*fac(1:180).*rhoRatio./Re./Re./Re./Re),'k')
%     plot(log10(fac5(1:180).*fac(1:180).*rhoRatio./Re./Re./Re./Re./Re),'m')
%     hold off
