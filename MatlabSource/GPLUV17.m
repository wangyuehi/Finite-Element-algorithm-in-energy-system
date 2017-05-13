% V14: bug fixed, using As save L&U, removed reach
% V15: added col dependency part
% V16: removed fliplr
% V17: fill-in changed to sparse traverse,removed L and U
% V17: warning: non-Zeros of As is nonZeros of L and U, merged As and Atest. 
%      but other entries in As is original value of As
function [x1] = GPLUV17(A1,BB)
% return value: x1 : solution of A1*x = BB
%               LevelOf: the level of every col. index is col number,
%                        entries are level of this column.
%               ancester: the ith entry is the column nunber that i column depend on, with 
%                         highest level 
  B = BB;
 n = size(A1,1);
 levelOf = zeros(n,1)+1;
 ancester = zeros(n,1);

 As = A1;





nonZeroIndexLUtest = zeros(n,460);
s = 1;
for i = 1:n
    
    if(A1(i,1)~=0)
        nonZeroIndexLUtest(1,s) = i;
        s = s + 1;
    end
end

 
    for i = 2:n %col
        s = 1;
        for j = 1:n % row
            
            if (As(j,i)~=0) % go to col J
                As(j,i) = 1;
                nonZeroIndexLUtest(i,s) = j;
                s = s + 1;
                if(j<i)
                    %compute dependency
                    if(levelOf(j)+1>levelOf(i))
                        ancester(i) = j;
                        levelOf(i) = levelOf(j)+1;
                        
                    end
                    
           %find fill in         
                    
                    for m = 1:460
                        in = nonZeroIndexLUtest(j,m);
                        if(in~=0)
                            As(in,i) = 1;
                        else
                            break;
                        end
                    end
                        
                  
                end
            end
        end
        
    end
    
    
for k = 1:n
    b = A1(:,k);

   %-------------- Solve  Lx = b---------------;
   index = nonZeroIndexLUtest(k,:);
   x1 = b;

   for m = 1:size(index,2)
%        size(index,2)
         j = index(m);
  %     j=m;

    if(j>0)


%if j>=k, only one nonzero entry is (j,j), do nothing.
if (j>=k)
    continue;
else
    %if j<k, the nonzero entries are stored in nonzeroindexLU
    for ii = 1:460
        nonZIndex = nonZeroIndexLUtest(j,ii);
        if(nonZIndex~=0)
            if(nonZIndex<=j)
                continue;
            else % nonZIndex>j
                x1(nonZIndex) = x1(nonZIndex) - As(nonZIndex,j)*x1(j);
            end
        else
            break;
        end
    end
end
        
    
    else
        break;
    end
   end
 
      
for i = 1:size(index,2)
    j = index(i);
    if(j>0)
    if (j<=k)
        As(j,k) = x1(j);
    else
        As(j,k) = x1(j)/As(k,k); 
    end
    else
        break;
    end
end
 
end


% % % 
% % % 
%----backward substitution-------
    %----first solve Ly = b--------
    y1 = zeros(1,n);
    for k = 1:n
        y1(k) = B(k);
        
        % L use indices larger than k
        nonZeroIndexL = nonZeroIndexLUtest(k,:);
%         nonZeroIndexL = nonZeroIndexL(nonZeroIndexL>0);
% % %         %---debug
% % %         nonZeroIndexLindex = find(nonZeroIndexLU(k,:)>k);
% % %         nonZeroIndexL = fliplr(nonZeroIndexLU(k,nonZeroIndexLindex));
        %----
        for i = 1:size(nonZeroIndexL,2)% here should be the non-zero ones   

               
                j = nonZeroIndexL(i);
                if(j==0)
                    break;
                end
                if(j<k)
                    continue;
                else 
                    if( j == k)
                         B(k)=B(k)-(y1(k));
                    else
                        B(j)=B(j)-(y1(k)*As(j,k));
                    end
                end
                
                
              
              
                
        end
       
    end
% %   
    x1 = zeros(1,n);
    %---- then solve Ux = y--------
    for k = n:-1:1
        x1(k) = y1(k)/As(k,k);
        
        nonZeroIndexU = nonZeroIndexLUtest(k,:);% WARNING: this is in descending order!
%         nonZeroIndexU = nonZeroIndexU(nonZeroIndexU>0);
        %has removed indices larger than k.
       
        for i = 1:size(nonZeroIndexU,2)
            j = nonZeroIndexU(i);
            if(j==0)
                break;% end of indices. 
            else
                
                y1(j) = y1(j) - (x1(k)*As(j,k));
               
            end
            
        
        end
    end
  end