clc;clear;
% an example of how to use the GetOneMap() function

maps = {};
map = cell(3);
ECs = {};
EDs = {};

[EDs, ECs, maps] = getOneMap(EDs, ECs, maps,map,[0,0]);

% 测试有没有重复映射
equalIndex = {};

for i = 1:1996
    for j = i+1:1996
        tmpEC = ECs{i} -ECs{j};
        tmpED = EDs{i} -EDs{j};
        if sum(sum(tmpEC~=0)) == 0 && sum(sum(tmpED~=0)) == 0
           equalIndex{end+1} = [i,j]; 
        end
    end
end