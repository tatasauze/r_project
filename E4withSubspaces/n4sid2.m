function [TH,bestchoice,nchoice,failflag] = ...
          n4sid(z,order,l,auxord,dkx,maxsize,Tsamp,refine,arg,trace)
%N4SID   Estimates a state-space model using a sub-space method.
%       TH=N4SID(Z) or [TH,AO]=N4SID(Z,ORDER,NY,AUXORD,DKX,MAXSIZE,TSAMP)
%
%	TH: Returned as the estimated state-space model in the THETA format.
%           No model covariances are given. 
%	Z : The output input data [y u], with y and u as column vectors
%           For multi-variable systems, Z=[y1 y2 ... yp u1 u2 ... un]
%	ORDER: The order of the model (Dimension of state vector). If entered
%	    as a vector (e.g. 3:10) information about all these orders will be
%	    given in a plot, Default; ORDER=1:10;
%	    If ORDER is entered as 'best', the default order among 1:10 is
%	    chosen.
%	NY: The number of outputs in the data matrix. Default NY =1.
%	AUXORD: An auxiliary order, that is used for the selection of state
%           variables. Default 1.2*ORDER+3. If AUXORD is entered as a row vector
%	    the best value (min pred error) in this vector will be selected.
%	DKX: This is a three-dimensional vector: DKX =[D,K,X]
%           D=1 indicates that a direct term from input to output will be 
%               estimated, while D=0 means that a delay from input to output
%	        is postulated.
%	    K=1 indicates that the K-matrix is estimated, while K=0 means that 
%               K will be fixed to zero. 
%           X=1 indicates that the initial state is estimated, X=0 that the
%	        initial state is set to zero.
%           Default: DKX = [0, 1, 0]
%       TRACE: Letting the last given argument be 'trace' gives info to screen
%	    about fit and choice of AUXORD
%       MAXSIZE: See also AUXVAR.
%
%	AO: The chosen value of AUXORD.
%
% 	The algorithm implements Van Overschee's and De Moor's method for 
%	identification of general multivariable linear systems in state space.
%	See also CANSTART, PEM.

% 	M. Viberg, 8-13-1992, T. McKelvey, L. Ljung 9-26-1993.
%	Copyright (c) 1986-95 by the MathWorks, Inc.
%	$Revision: 1.11 $  $Date: 1995/06/01 10:19:25 $

global XIDplotw
if nargin<1
  disp('Usage: TH = N4SID(Z);')
  disp('       TH = N4SID(Z,ORDER,NY,AUXORD,DKX,MAXSIZE,TSAMP,''trace'');')
  return
end
[Ncap,NN] = size(z);
if NN>Ncap,error('Data should be organized in column vectors.'),end
failflag=0;
trace=0;
if nargin==10,if strcmp(lower(trace),'trace'),trace=1;end,end
if nargin==9,if strcmp(lower(arg),'trace'),trace=1;end,end
if nargin==8,if strcmp(lower(refine),'trace'),trace=1;end,end
if nargin==7,if strcmp(lower(Tsamp),'trace'),trace=1;end,end
if nargin==6,if strcmp(lower(maxsize),'trace'),trace=1;end,end
if nargin==5,if strcmp(lower(dkx),'trace'),trace=1;end,end
if nargin==4,if strcmp(lower(auxord),'trace'),trace=1;end,end
if nargin==3,if strcmp(lower(l),'trace'),trace=1;end,end
if nargin==2,if strcmp(lower(order),'trace'),trace=1;end,end
if nargin < 9, arg='nogui';end
%if strcmp(arg,'gui')
%   eval('zv=iduigetd(''v'');')
%else
   zv=z;
%end

[Ncapv,nz_zv] = size(zv);
if nargin < 8, refine=[];end,if isempty(refine),refine=1;end
if nargin < 7, Tsamp=[];end
if nargin < 6, maxsize=[];end
if nargin < 5, dkx=[];end,if isempty(dkx),dkx=[0,1,0];end
if nargin < 4, auxord=[];end
if nargin < 3, l=[];end, if isempty(l),l=1;end
if nargin < 2, order=[];end,if isempty(order),order=1:10;end
m=NN-l;
if isstr(order)
    if Ncap<19*l+30*(m+1)+1
       error('There are too few data points to support this order choice.')
    end
    def_order=1;order=1:10;
else 
   def_order=0;
end
if dkx(2),oe='no';else oe='oe';end
if dkx(1),dmat='d';else dmat='no';end
if isempty(Tsamp),Tsamp=1;end
if length(auxord>1),tryaux=1;else tryaux=0;end
if isempty(trace)
   if tryaux
      trace=1;
   else
      trace=0;
   end
end
if m==0,timeseries=1;else timeseries=0;end
if timeseries&dkx(2)==0
   error('For a time-series model, the K-matrix has to be estimated.')
end
auxord=auxord(:)';
if l>NN
    error('You have specified more outputs than the number of data columns.')
end
if length(order)>1,
   n=max(order);
   ordchoice=1;
   if tryaux
      error('The chosen ORDER and the chosen AUXORD cannot both be vectors.') 
   end
else
   n=order;
   ordchoice=0;
   laux=length(auxord);
   if laux>1
      if strcmp(lower(oe),'oe'),ncheck=n+1;else ncheck=n+2;end
      auxord=auxord(find(auxord>=ncheck&auxord<=(Ncap+1)/6));
      if laux>length(auxord)&strcmp(arg,'nogui'), disp(['AUXORD changed to ',...
       int2str(min(auxord)),':',int2str(max(auxord))]),end
   end
end

if isempty(auxord)
      auxord=fix(1.2*n)+3;
end

if isempty(maxsize),maxsize=idmsize(Ncap,2*auxord(1)*(l+m));end

y = z(:,1:l);
u = z(:,l+1:l+m);
bestloss=[];
if length(auxord)>1&strcmp(arg,'gui')
   handax=get(XIDplotw(11,1),'userdata');handax=handax(3);
   axes(handax),cla,xlabel('Auxiliary order'),ylabel('Fit to working data')
   set(gca,'xlim',[min(auxord)-1,max(auxord)+1])
   set(XIDplotw(11,1),'vis','on')
end
for icount=auxord %
if strcmp(arg,'guichoice')
   eval('hh=findobj(XIDplotw(10,1),''label'',menulabel(''&Help''));')
   R=get(hh,'userdata');
   [nR,cR]=size(R);
   l=R(1,1);i=R(1,2);m=R(1,3);Ncap=R(1,4);
   j = Ncap+1-2*i;
   if m==0,timeseries=1;else timeseries=0;end
   R=R(2:nR,:);
   n=order;
else
   i=icount;
   if strcmp(lower(oe),'oe'),ncheck=n+1;else ncheck=n+2;end
   j = Ncap+1-2*i;
   if j<2*(m+1)*i
      itemp = fix((Ncap-l)/(1.5*(l+m)+m));
      if itemp<ncheck
         error('Too few data points for this choice of orders.')
      end
      if itemp<i
         i = itemp;
         disp(['AUXORD has been changed to ' num2str(i)])
         j = Ncap+1-2*i;
      end
   end
   if i<ncheck
      i=ncheck;
      disp(['AUXORD has been changed to ',int2str(i)])
   end
   nrr=2*i*(l+m);
   if nrr*1.2*nrr>maxsize,
      maxsize=nrr*1.2*nrr;
      disp(['MAXSIZE has been changed to ',int2str(maxsize)]);
   end
   %%% Hankel matrices:
   nele=floor(maxsize/(2*i*(l+m)));
   if nele==0,error(['MAXSIZE must be larger than ',int2str(maxsize) '.']),end
   nloop=floor((Ncap-2*i)/nele-eps)+1;
   R=[];H1=zeros(nrr,nrr+nele+2*i-1);
   for kk=1:nloop
        jj=[1+(kk-1)*nele:min(Ncap,kk*nele+2*i-1)];
        Y = ssssaux('idblockh',y(jj,:)/sqrt(j),2*i);
        if ~timeseries
            U = ssssaux('idblockh',u(jj,:)/sqrt(j),2*i);
        else
            U=[];
        end
        H = [U;Y];
        if ~isempty(R)
          H1(:,1:min(ncR,nrr)+length(jj)-2*i+1) =...
              [H,R(1:min(nrR,2*i*(l+m)),1:min(ncR,2*i*(l+m)))];
        else
          H1(:,1:length(jj)-2*i+1)=H;ncR=0;
        end
        R=triu(qr(H1(:,1:min(ncR,nrr)+length(jj)-2*i+1)'))';[nrR,ncR]=size(R);
   end  %for kk
   if ncR<i*(2*m+l)+l
         failflag=2;
         disp('Error: Too few data points for this choice of orders.')
         return
   end
end % if strcmp(arg,'guichoice')
%%% Z_i:
       if ordchoice&~def_order
        if strcmp(arg,'gui')
           handax=get(XIDplotw(10,1),'userdata');handax=handax(3);
           set(handax,'vis','off');axes(handax),cla
           [nR,nC]=size(R);
           eval('hh=findobj(XIDplotw(10,1),''label'',menulabel(''&Help''));')
           set(hh,'userdata',[[l,i,m,Ncap,zeros(1,nC-4)];R]);
        else
          hfig=figure;
        end
       end
if timeseries,
     %%% Split R
     R41t3 = R(l*(i+1)+1:l*2*i,1:l*(i+1));
     R44   = R(l*(i+1)+1:l*2*i,l*(i+1)+1:l*2*i);
     R3t41t2 = R(l*i+1:l*2*i,1:l*i);
     R3t43t4 = R(l*i+1:l*2*i,l*i+1:l*2*i);
     R2t41   = R(l*(i-1)+1:l*2*i,1:l*(i-1));
     R2t42t4 = R(l*(i-1)+1:l*2*i,l*(i-1)+1:l*2*i);

     %%% Calculate QSVD's
     [Um1,Sm1,Xm1,Vm1,Tm1] = ssssaux('gsvd',R41t3',R44');
     [U,S,X,V,T]           = ssssaux('gsvd',R3t41t2',R3t43t4');
     [Up1,Sp1,Xp1,Vp1,Tp1] = ssssaux('gsvd',R2t41',R2t42t4');

Sm1
S
Sp1
   if ordchoice
      dS=diag(S);
      testo=log(diag(S));
      ndef=max(find(testo>(max(testo)+min(testo))/2));
      if ~def_order
      [xx,yy]=bar(1:length(dS),log(dS)');
      mdS=floor(min(log(dS)));
      zer=find(yy==0);yy(zer)=mdS*ones(size(zer));
      if ~strcmp(arg,'gui')
          handax=gca;
      else
          hfig=XIDplotw(10,1);
      end
      axes(handax)
      line(xx,yy,'color','y');%%
      ylabel('Log of Singular values');xlabel('Model order')
      title('Model singular values vs order')
      text(0.97,0.95,'Red: Default Choice','units','norm','fontsize',10,...
          'HorizontalAlignment','right');
      set(hfig,'vis','on')
      patch(xx(1:ndef*5-3),yy(1:ndef*5-3),'y');%%
      patch(xx(5*ndef-3:5*ndef),yy(5*ndef-3:5*ndef),'r');%%%
      patch(xx(5*ndef:n*5+1),yy(5*ndef:n*5+1),'y');%%
      set(handax,'vis','on')
      if strcmp(arg,'gui')
         set(XIDplotw(10,3),'userdata',[[ndef,1:length(dS)];[dS(ndef),dS']]);
         return
      end
      title('Select model order in Command Window.')
      n=input('Select model order:(''Return'' gives default) ');
      if isempty(n)
         n=ndef;
         disp(['Order chosen to ',int2str(n)]);
      end
      else
         n=ndef;
      end % if def_order
      nchoice=n;
   end

     %%% Find matrices

     A = inv(sqrt(S(1:n,1:n)))*pinv(X(1:l*(i-1),1:n))*Xm1(:,1:n)*...
         Sm1(1:n,1:n)*Um1(1:l*i,1:n)'*U(:,1:n)*inv(sqrt(S(1:n,1:n)))

     CC = X(:,1:n)*sqrt(S(1:n,1:n));
     C = CC(1:l,:)

     GG = R(1:l*i,1:l*i)*U(:,1:n)*sqrt(S(1:n,1:n));

     G = GG(l*(i-1)+1:l*i,:)';

 G = G*C
C = 1;

     Lam0 = R(l*i+1:l*(i+1),1:l*(i+1))*R(l*i+1:l*(i+1),1:l*(i+1))';

     %% Solve the Riccati equation (8)

     iLam = inv(Lam0);
     [W,d] = eig([ A'-C'*iLam*G' zeros(n,n);...
          -G*iLam*G'   eye(n) ], ...
        [eye(n) -C'*iLam*C; ...
         zeros(n,n) A-G*iLam*C ]);
     d = diag(d);
     [e,index] = sort(abs(d));        % sort on magnitude of eigenvalues
     if (~((e(n) < 1) & (e(n+1)>1)))
           disp('Can''t order eigenvalues');
     end

    % select vectors with eigenvalues inside unit circle
    WW = W(:,index(1:n));

    %%% State covariance:
    P = WW(n+1:2*n,:)*inv(WW(1:n,:));

    Q = P-A*P*A';
    S = G-A*P*C';
    R = Lam0-C*P*C';
S*inv(R)

    oe='not_oe';
    B=[];
    D=[];

else  % i.e if not timeseries

   invma= R((1:2*i*m+i*l),(1:2*i*m+i*l));
   if rank(invma)<i
    disp('WARNING: Your input signal is apparently not persistently exciting.')
    disp('         Use other input or lower model order.')
    disp('         No model returned.')
    TH=[];failflag=1;return
   end

   L = (R(2*m*i+l*i+1:2*m*i+2*l*i,1:2*i*m+i*l)/invma);
   Li1 = L(:,1:m*i);
   Li3 = L(:,2*m*i+1:2*m*i+l*i);

%%% Gamma_i:
   [UU,SS,VV] = svd([Li1 zeros(l*i,m*i) Li3]*R((1:2*i*m+i*l),(1:2*i*m+i*l)));
   if ordchoice
      dS=diag(SS);
      testo=log(diag(SS));
      ndef=max(find(testo>(max(testo)+min(testo))/2));
      if ~def_order
      [xx,yy]=bar(1:length(dS),log(dS)');
      mdS=floor(min(log(dS)));
      zer=find(yy==0);yy(zer)=mdS*ones(size(zer));
      if ~strcmp(arg,'gui')
          handax=gca;
      else
          hfig=XIDplotw(10,1);
      end
      axes(handax)

      line(xx,yy,'color','y');%%
      ylabel('Log of Singular values');xlabel('Model order')
      title('Model sigular values vs order')
      text(0.97,0.95,'Red: Default Choice','units','norm','fontsize',10,...
          'HorizontalAlignment','right');
      set(hfig,'vis','on')
      patch(xx(1:ndef*5-3),yy(1:ndef*5-3),'y');%%
      patch(xx(5*ndef-3:5*ndef),yy(5*ndef-3:5*ndef),'r');%%%
      patch(xx(5*ndef:n*5+1),yy(5*ndef:n*5+1),'y');%%
      set(handax,'vis','on')

     if strcmp(arg,'gui')
         set(XIDplotw(10,3),'userdata',[[ndef,1:length(dS)];[dS(ndef),dS']]);
         return
      end
      title('Select model order in command window.')
      n=input('Select model order:(''Return'' gives default) ');
      if isempty(n)
         %testo=log(diag(SS));
         n=ndef;
         disp(['Order chosen to ',int2str(n)]);
      end
      else
         n=ndef;
      end  % if def_order
      nchoice=n;
   end
   U1 = UU(:,1:n); Sigma1 = SS(1:n,1:n);

%%% LS problem:

   BB0 = [U1(1:l*(i-1),:)\R(2*m*i+l*i+l+1:2*m*i+2*l*i,1:2*i*m+(i+1)*l);
	R(2*m*i+l*i+1:2*m*i+l*i+l,1:2*m*i+l*(i+1))];
   BB = BB0(:,1:2*i*m+i*l);
   AA = [U1'*R(2*m*i+l*i+1:2*m*i+2*l*i,1:2*i*m+i*l);
	R(m*i+1:2*m*i,1:2*m*i+l*i)];
   K = BB/AA;

%%% Extraction of system matrices:

   A = K(1:n,1:n);


   C = K(n+1:n+l,1:n);
   G1 = U1(1:l,:)'; G2 = U1(l+1:l*i,:)'; U2 = U1(1:l*(i-1),:);
%   Gamma = obs_mat(A,C,i); G = pinv(Gamma); G1 = G(:,1:l); 
%   G2 = G(:,l+1:l*i); U2 = Gamma(1:l*(i-1),:)
   S1 = [-A*G1 eye(n)-A*G2*U2 ; eye(l)-C*G1 -C*G2*U2];
   S2 = [ssssaux('idspech',pinv(U2)-A*G2,l) ; ssssaux('idspech',-C*G2,l)];
   S2 = S2*[eye(l) zeros(l,n); zeros((i-2)*l,l) U2(1:l*(i-2),:)];
 k = [K(:,n+1:n+m); ssssaux('idblockc',K(1:n,n+m+1:n+m*i),m);...
                    ssssaux('idblockc',K(n+1:n+l,n+m+1:n+m*i),m)];
   DB = [S1 ; S2]\k;
   D = DB(1:l,:);
   B = DB(l+1:l+n,:);
end


 if ~timeseries
   res = BB - K*AA;
   BBrest=BB0(:,2*i*m+i*l+1:2*i*m+(i+1)*l);
   QRS=res*res'+BBrest*BBrest';
   Q = QRS(1:n,1:n);
   S = QRS(1:n,n+1:n+l);
   R = QRS(n+1:n+l,n+1:n+l);
 end
 errflag=0;
 eval('K1=ssssaux(''kric'',A,C,Q,R,S);','errflag=1;')
 if errflag,failflag=1;return,end
if ~strcmp(lower(oe),'oe')		
   K=K1;
else
   K=zeros(n,l);
end

THprel = ms2th(real(modstruc(A,B,C,D,K)));
if tryaux
  e=pe(zv,THprel);
  Vloss=det(e'*e/Ncapv);
  if trace,
     mess=['auxord = ',int2str(i),' gives a fit of ',num2str(Vloss)];
     disp(mess)
  end
  if strcmp(arg,'gui')
        [xx,yy]=bar(i,Vloss);
        patch(xx,yy,'y');drawnow
  end %if 'gui'

  if isempty(bestloss),
     update=1;
  elseif Vloss<bestloss
     update=1;
  else
     update=0;
  end
  if update
     TH=THprel;bestloss=Vloss;bestchoice=i;
  end
else
     TH=THprel;bestchoice=i;
end % if tryaux
end % for icount=auxord
if tryaux&strcmp(arg,'gui')
   [xx,yy]=bar(bestchoice,bestloss);
   patch(xx,yy,'red')
end
e=pe(z,TH);

refine = 0;

if refine&~timeseries % refine the estimate of B and D
'hola'
   [A,B,C,D,K,X0]=th2ss(TH);
   lam=e'*e/Ncap;sqrlam=inv(sqrtm(lam));
   [ny,nx]=size(C);
   [nx,nu]=size(B);nz=ny+nu;
   n=nu*nx+nx;if strcmp(dmat,'d'),n=n+ny*nu;end
   % *** Compute the gradient PSI. If N>M do it in portions ***
   rowmax=max(n*ny,nx+nz);
   M=floor(maxsize/rowmax);
   R=zeros(n,n);Fcap=zeros(n,1);
   for kc=1:M:Ncap-1
      jj=(kc:min(Ncap,kc-1+M));
      if jj(length(jj))<Ncap,jjz=[jj,jj(length(jj))+1];else jjz=jj;end
      psitemp=zeros(length(jj),ny);
      psi=zeros(ny*length(jj),n);
      x=ltitr(A-K1*C,K1,z(jjz,1:ny),X0); 
           %We use the good K even for an OE model
      yh=(C*x(1:length(jj),:)')';
      e=(z(jj,1:ny)-yh)*sqrlam;
      [nxr,nxc]=size(x);X0=x(nxr,:)';
      evec=e(:);
      kl=1;
      for kx=1:nx,for ku=1:nu
          dB=zeros(nx,1);dB(kx,1)=1;
          if kc==1,dX=zeros(nx,1);else dX=dXk(:,kl);end
          psix=ltitr(A-K1*C,dB,z(jjz,ny+ku),dX);[rp,cp]=size(psix);
	  dXk(:,kl)=psix(rp,:)';
	  psitemp=(C*psix(1:length(jj),:)')'*sqrlam;
          psi(:,kl)=psitemp(:);kl=kl+1;	
      end,end
      if strcmp(dmat,'d'),
          for ky=1:ny,for ku=1:nu
          psitemp=...
           [zeros(length(jjz),ky-1),z(jjz,ny+ku),zeros(length(jjz),ny-ky)]...
           *sqrlam;psitemp=psitemp(1:length(jj),:);
           psi(:,kl)=psitemp(:);kl=kl+1;
          end,end
       end
       %% x0
       for kx=1:nx
          if kc==1
             x0dum=zeros(nx,1);x0dum(kx,1)=1;
          else
             x0dum=X00(:,kl);
          end
          psix=ltitr(A-K1*C,zeros(nx,1),zeros(length(jjz),1),x0dum);
          [rp,cp]=size(psix);
	  X00(:,kl)=psix(rp,:)';
	  psitemp=(C*psix(1:length(jj),:)')'*sqrlam;
          psi(:,kl)=psitemp(:);kl=kl+1;	
       end
      R=R+psi'*psi;  Fcap=Fcap+ psi'*evec;
   end

  % *** Compute the estimate of B and D ***

   if Ncap>M, g=R\Fcap; else g=psi\evec;end
   B=reshape(g(1:nx*nu),nu,nx)';
   if strcmp(dmat,'d'),
      D=reshape(g(nx*nu+1:n-nx),nu,ny)';
   else 
      D=zeros(ny,nu);
   end
   if dkx(3),x0=g(n-nx+1:n);else x0=zeros(nx,1);end
   TH = ms2th(real(modstruc(A,B,C,D,K,x0)));
   e=pe(z,TH);
end % if refine
lambda=e'*e/Ncap;Vloss=det(lambda);
npar=l*n+m*n;if ~strcmp(lower(oe),'oe'),npar=npar+l*n;end
Vfpe=Vloss*(1+npar/Ncap)/(1-npar/Ncap);
TH(1,1)=Vloss;TH(2,1)=Vfpe;
TH(2,7)=50;TH=sett(TH,Tsamp);
rarg=TH(1,6);
TH(4+rarg:3+l+rarg,1:l)=lambda;

