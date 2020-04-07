function bigkmat = bigKmat(kmat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[s1,s2] = size(kmat);


bigkmat = zeros(s1*2+1,s2*2+1);
% bigkmat opvullen met kmat
for i = 1:(s1*2+1)
    for j = 1:(s2*2+1)
        if (mod(i,2) == 0)&&(mod(j,2)==0)
            bigkmat(i,j) = kmat(i/2,j/2);
        end  
    end
end
% gemiddeldes nemen in bigkmat
% speciale gevallen:
% 1. inwendige punten (punten die tussen de waardes van de originele kmat
% liggen)
% 2. randpunten (die geen hoekpunt zijn)


% 1.

for i = 2:2:(s1*2) %(tijdens elke iteratie het gemiddelde nemen tussen element [i,j] en [i+2,j] en tussen [i,j] en [i,j+2] (element rechts en onder)
    for j = 2:2:(s2*2)
        if (i==s1*2)&&(j~=s2*2)
            bigkmat(i,j+1) = (bigkmat(i,j)+bigkmat(i,j+2))/2;
        elseif(j==s2*2)&&(i~=s1*2)
            bigkmat(i+1,j) = (bigkmat(i,j)+bigkmat(i+2,j))/2;
        elseif(i<=s1*2-2)&&(j<=s2*2-2)
            bigkmat(i+1,j) = (bigkmat(i,j)+bigkmat(i+2,j))/2;
            bigkmat(i,j+1) = (bigkmat(i,j)+bigkmat(i,j+2))/2;
        end
        
    end
end

% 2.
% eerste en laatste kolom
for i = 2:2:2*s1
    bigkmat(i,1) = bigkmat(i,2);
    bigkmat(i,(2*s2+1)) = bigkmat(i,(2*s2));
end

% eerste en laatste rij
for j = 2:2:2*s2
    bigkmat(1,j) = bigkmat(2,j);
    bigkmat((2*s1+1),j) = bigkmat((2*s1),j);
end
end

