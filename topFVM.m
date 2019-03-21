function topFVM(N,volfrac,penal,rmin)
% INITIALIZE
x(1:N,1:N) = volfrac; 
loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=heattop(N,reshape(x,N^2,1),penal);       
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 2:N-1
    for elx = 2:N-1
      n1 = N*(elx-2)+ely; 
      n2 = N* (elx-1) +ely;
      n3 = N* elx   +ely;
      Ue = U([n1; n2-1; n2+1; n3],1);
      c = c + (0.001+0.999*x(ely,elx)^penal)*Ue'*KE*Ue;
      dc(ely,elx) = -0.999*penal*x(ely,elx)^(penal- 1)*Ue'*KE*Ue;
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(N,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(N,x,volfrac,dc); 
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(N*N)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
  colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end 

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(N,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  if sum(sum(xnew)) - volfrac*N*N > 0
    l1 = lmid;
  else
    l2 = lmid;
  end
end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(N,rmin,x,dc)
dcn=zeros(N,N);
for i = 2:N-1
  for j = 2:N-1
    sum=0.0; 
    for k = max(i-floor(rmin),2):min(i+floor(rmin),N-1)
      for l = max(j-floor(rmin),2):min(j+floor(rmin),N-1)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*dc(l,k)*x(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
KE = [ 2/3 -1/6 -1/3 -1/6
      -1/6 2/3 -1/6 -1/3
      -1/3 -1/6 2/3 -1/6
      -1/6 -1/3 -1/6 2/3];

  
  

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(N,x,penal)
[KE] = lk; 
K = sparse((N +1) *(N +1), (N +1) *(N +1));
F = sparse((N +1) *(N +1),1); U = sparse((N +1) *(N +1),1);
for elx = 1:N
  for ely = 1:N
    n1 = (N+1)*(elx-1)+ely; 
    n2 = (N+1)* elx   +ely;
    edof = [n1; n2; n2+1; n1+1];
    K(edof,edof) = K(edof,edof) + (0.001+0.999*x(ely,elx)^penal)*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (SQUARE PLATE WITH HEAT SINK)
F(:,1) = 0.01;
%fixeddofs = [N/2+1-(N/20):2:N/2+1+(N/20)];
fixeddofs = [N/2+1-(N/20):2:N/2+1+(N/20);(N+1)*(N+1)+1-(N/2+1-(N/20):2:N/2+1+(N/20))];
alldofs = [1:(N+1)*(N+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;

  
  
  
  
%   
%           nlist = [];
%         if elx > 1
%             nlist = [nlist; (N)*(elx-2)+ely];
%         end
%         if ely > 1      
%             nlist = [nlist; (N)*(elx-1)+ely-1];
%         end
%         if ely < N
%             nlist = [nlist; (N)*(elx-1)+ely+1];
%         end
%         if elx < N
%             nlist = [nlist; (N)*(elx)+ely];
%         end
%       Ue = U(nlist,1);
%       if length(nlist)==3
%           K = K3;
%       elseif length(nlist)==2
%           K = K2;
%       else
%           K = K4;
%       end