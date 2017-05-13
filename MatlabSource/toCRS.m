function [CRSformat] =  toCRS(A)

n = size(A,1);
%should be changed to an approx value
s= size(nonzeros(A),1);

%----------------
%    value      |
%----------------
%  col index    |
%----------------
%  row pointer  |
%----------------
CRSformat = zeros(3,s);




index = 1;

for i = 1:n  % row
    CRSformat(3,i) = index;
    for j = 1:n
        if(A(i,j) ~= 0)
            
            CRSformat(1,index) = A(i,j);
            CRSformat(2,index) = j;
           
            index = index+1;
        end
    end
    if(i == n)
        CRSformat(3,i+1) = index;
    end
end

end
