function [c] = costfunc(temp,bigk, N)

% temp = solution matrix

dp = N+1;
 
c = 0;

hx = 1/N;
hy = 1/N;




%inwendige punten
for i = 2:(dp-1)
    for j = 2:(dp-1)   
        kv =0.25*(bigk(2*j-1-1,2*i-1-1) + bigk(2*j-1+1,2*i-1-1) + bigk(2*j-1+1,2*i-1+1) + bigk(2*j-1-1,2*i-1+1);
        xdir = ((temp(j,(i+1))-temp(j,(i-1)))/(2*hx))^2;
        ydir = ((temp((j+1),i)-temp((j-1),i))/(2*hy))^2;
        c = c + (xdir+ydir)*hx*hy*kv;  
    end
end


%neumann boven & onder

for j = 2:(dp-1)
    %boven
    kv = 0.5*(bigk(1+1,2*j-1-1) + bigk(1+1,2*j-1+1));
    xdir = ((temp(1,(j+1)) - temp(1,(j-1)))/(2*hx))^2;
    ydir = (temp(2,j)/(hy/2))^2;
    c = c + kv*(xdir+ydir)*hx*hy/2;
    
    %onder
    kv = 0.25*(bigk(2*dp-1-1,2*j-1-1) + bigk(2*dp-1-1,2*j-1+1));
    xdir = ((temp(dp,(j+1)) - temp(dp,(j-1)))/(2*hx))^2;
    ydir = (temp((dp-1),j)/(hy/2))^2;
    c = c + kv*(xdir+ydir)*hx*hy/2;
end

%hoekpunten
%linksboven 
kv = bigk(1+1,1+1);
xdir = (temp(1,2)/(0.5*hx))^2;
ydir = (temp(2,1)/(0.5*hy))^2;
c = c + kv*(xdir+ydir)*hx*hy/4;
%linksonder
kv = bigk(2*dp-1-1,1+1);
xdir = (temp(dp,2)/(0.5*hx))^2;
ydir = (temp(dp-1,1)/(0.5*hy))^2;
c = c + kv*(xdir+ydir)*hx*hy/4;
%rechtsboven
kv = bigk(1+1,2*dp-1-1);
xdir = (temp(1,dp-1)/(0.5*hx))^2;
ydir = (temp(2,dp)/(0.5*hy))^2;
c = c + kv*(xdir+ydir)*hx*hy/4;
%rechtsonder
kv = bigk(2*dp-1-1, 2*dp-1-1);
xdir = (temp(dp,dp-1)/(0.5*hx))^2;
ydir = (temp(dp-1,dp)/(0.5*hy))^2;
c = c + kv*(xdir+ydir)*hx*hy/4;


%links&rechts
for j = 2:(dp-1)
    %links
    kv = 0.5*(bigk(2*j-1-1,1+1) + bigk(2*j-1+1,1+1));
    xdir = (temp(j,1+1)/(hx/2))^2;
    ydir = ((temp(j-1,1)+temp(j+1,1))/hy)^2;
    c = c + kv*(xdir+ydir)*hx*hy/2;
    
    %rechts
    kv = 0.5*(bigk(2*j-1-1,2*dp-1-1) + bigk(2*j-1+1,2*dp-1-1));
    xdir = (temp(j,dp-1)/(hx/2))^2;
    ydir = ((temp(j-1,dp)+temp(j+1,dp))/hy)^2;
    c = c + kv*(xdir+ydir)*hx*hy/2;
end

