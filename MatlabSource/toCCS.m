function [CCSformat,d] =  toCCS(A)

n = size(A,1);
d = n;
%should be changed to an approx value


s= size(nonzeros(A),1);

%----------------
%    value      |
%----------------
%  row index    |
%----------------
%  col ptr      |
%----------------
CCSformat = zeros(3,s);




index = 1;

for i = 1:n  % col
CCSformat(3,i) = index;
    for j = 1:n %row
        if(A(j,i) ~= 0)
            
            CCSformat(1,index) = A(j,i);
            CCSformat(2,index) = j;

            index = index+1;
        end
    end
    if(i == n)
        CCSformat(3,i+1) = index;
    end
    
end

end
