function [ EDs, ECs, maps ] = getOneMap( EDs, ECs, maps, map , pt)
% pt  =  [i,j]

% map =  { [] [] []
%          [] [] []
%          [] [] [] }

% idx = [ 7 8 9
%         4 5 6
%         1 2 3 ]

% p = [(0,0) (0,1) (0,2)
%      (1,0) (1,1) (1,2)
%      (2,0) (2,1) (2,2)];

idx = (pt(1) + 3 * pt(2)) + 1;


for i = pt(1)-1:pt(1)
    for j = pt(2)-1:pt(2)
        if i >=0 && j>=0
            map{idx} = [i,j];
            
            if idx == 9
                ec = [ 0 0 2
                    0 0 1
                    2 1 1 ];
                ed = [0 0 3
                    0 0 2
                    3 2 2];
                for m = 1:3
                    for n = 1:3
                        ec(map{m,n}(1)+1,map{m,n}(2)+1) = ec(map{m,n}(1)+1,map{m,n}(2)+1)+1;
                        ed(map{m,n}(1)+1,map{m,n}(2)+1) = ed(map{m,n}(1)+1,map{m,n}(2)+1)+...
                            sum(abs([map{m,n}(1)+1,map{m,n}(2)+1] - [m,n]));
                    end
                end
                
                if sum(sum(ec==0)) == 0
                    maps{end+1} = map;
                    ECs{end+1} = log2(ec);
                    EDs{end+1} = ed./ec;
                end
                continue;
            end
            
            new_idx = idx + 1;
            pt_new = [ mod(new_idx-1,3), floor((new_idx-1)/3)];
            [ EDs, ECs, maps ] = getOneMap(EDs, ECs, maps, map, pt_new);
        end
    end
end


end

