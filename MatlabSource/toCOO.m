function [COOformat] =  toCOO(A)

n = size(A,1);
%should be changed to an approx value


s= size(nonzeros(A),1);

%----------------
%    value      |
%----------------
%  col index    |
%----------------
%  row index    |
%----------------
COOformat = zeros(3,s);




index = 1;

for i = 1:n  % row

    for j = 1:n
        if(A(i,j) ~= 0)
            
            COOformat(1,index) = A(i,j);
            COOformat(2,index) = j;
            COOformat(3,index) = i;
            index = index+1;
        end
    end
end

end
