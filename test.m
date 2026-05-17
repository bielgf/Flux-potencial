function test()
    
%     Nw = 512;
%     Nh = 512/4;
% 
%     a1 = strings(Nw+Nh,Nw+Nh);
%     b1 = strings(Nw+Nh,1);

    N1 = 4;
    N2 = 2;

    a2 = strings(N1+N2,N1+N2);
    b2 = strings(N1+N2,1);

%     tic
% 
%     for ii = 1:Nw
%         b1(ii,1) = 'W';
%         for jj = 1:Nw
%             if ii == jj
%                 a1(ii,jj) = 'WW_ij';
%             else
%                 a1(ii,jj) = 'WW';
%             end
%         end
% 
%         for jj = 1:Nh
%             a1(ii,Nw+jj) = 'WC';
%         end
%     end
% 
%     for ii = 1:Nh
%         b1(Nw+ii,1) = 'C';
%         for jj = 1:Nh
%             if ii == jj
%                 a1(Nw+ii,Nw+jj) = 'CC_ij';
%             else
%                 a1(Nw+ii,Nw+jj) = 'CC';
%             end
%         end
% 
%         for jj = 1:Nw
%             a1(Nw+ii,jj) = 'CW';
%         end
%     end
% 
%     t1 = toc;

    tic

    for i=1:N1
        b2(i,1) = 'W';
        b2(N1+i,1) = 'C';
        for j=1:N1
            
            if j==i
                %Wing-Wing
                a2(i,j) = 'WW-ij';
                
            else
                %Wing-Wing
                a2(i,j) = 'WW';
            end 
            %Wing-Estabilizador
            a2(i,N2+j) = 'WC';
        end
        for j=1:N1
    
            if j==i
                %Estabilizador-Estabilizador
                a2(N1+i,N2+j) = 'CC-ij';
                
            else
                %Estabilizador-Estabilizador
                a2(N1+i,N2+j) = 'CC';
            end
            %Estabilizador-Wing
            a2(N1+i,j) = 'CW';
        end
    end

    t2 = toc;

end