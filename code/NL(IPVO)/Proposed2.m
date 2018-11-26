clear all; clc;
addpath(genpath('Origin Images')); addpath(genpath('result')); addpath(genpath('tools'));

mask2 = [false, true; true, true];
mask3 = [false, true, true; true, true, true; true, true, true];

Imgs = {'Lena', 'Baboon', 'Barbara', 'Airplane', 'Lake', 'Peppers', 'Boat', 'Elaine'};
%%
for tt = 1:1:8
    Iname = Imgs{tt};
    istr = ['Proposed_2019_',Iname,'.mat']
    I = double(imread([Iname,'.bmp']));
    [A, B] = size(I);
    
    %%
    EdgInfo = 40;
    a = 2; b = 2;
    R = zeros(2,1000);
    cnt = 0;
    for Payload =10000+EdgInfo:1000:100000+EdgInfo
        %         Payload
        %%
        H1 = cell(1,2048);
        for i = 1 : 2048
            H1{i} = zeros(1,256);
            H2{i} = zeros(1,256);
        end
        
        EC1 = zeros(A-2-a,B-2-b);
        ED1 = zeros(A-2-a,B-2-b);
        NL = zeros(A-2-a,B-2-b);
        for i = 1:(A-2)-a
            for j = 1:(B-2)-b
                for ii = 1:a+1
                    for jj = 1:b+2
                        if ii == a+1 || jj == b+1 || jj == b+2
                            NL(i,j) = NL(i,j) + abs(I(1*(i-1)+ii,1*(j-1)+jj) - I(1*(i-1)+ii+1,1*(j-1)+jj));
                        end
                    end
                end
                for ii = 1:a+2
                    for jj = 1:b+1
                        if ii == a+1 || ii == a+2 || jj == b+1
                            NL(i,j) = NL(i,j) + abs(I(1*(i-1)+ii,1*(j-1)+jj) - I(1*(i-1)+ii,1*(j-1)+jj+1));
                        end
                    end
                end
                T = NL(i,j);
                
                % 2x2 block
                X = I(i:i+1,j:j+1);
                X = X(mask2);
                [Y, In] = sort(X);
                if Y(end) ~= Y(1)
                    % dmax
                    if I(i,j) >= Y(end)
                        dmax = I(i,j) - Y(end); % >= 0
                        H1{T+1}(dmax+1) = H1{T+1}(dmax+1) + 1;
                        if dmax == 0
                            EC1(i,j) = EC1(i,j) + 1;
                            ED1(i,j) = ED1(i,j) + 0.5;
                        else
                            ED1(i,j) = ED1(i,j) + 1;
                        end
                    end
                    % dmin
                    if I(i,j) <= Y(1)
                        dmin = Y(1) - I(i,j); % >= 0
                        H1{T+1}(dmin+1) = H1{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC1(i,j) = EC1(i,j) + 1;
                            ED1(i,j) = ED1(i,j) + 0.5;
                        else
                            ED1(i,j) = ED1(i,j) + 1;
                        end
                    end
                else % In(end) == In(1)
                    if I(i,j) <= Y(end)
                        dmin = Y(1) - I(i,j); % >= 0
                        H1{T+1}(dmin+1) = H1{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC1(i,j) = EC1(i,j) + 1;
                            ED1(i,j) = ED1(i,j) + 0.5;
                        else
                            ED1(i,j) = ED1(i,j) + 1;
                        end
                    end
                end
                
            end
        end
        
        Ratio = 100000;
        nls = unique(NL(:));
        t1 = 0;
        for i = 1 : numel(nls)
            T1 = nls(i) + 1;
            h1 = zeros(1,256);
            for k = 1 : T1
                h1 = h1 + H1{k};
            end
            ec1 = h1(1);
            
            if ec1 > Payload
                ed1 = 0.5*h1(1) + sum(h1(2:end));
                ratio = ed1/ec1;
                if Ratio > ratio
                    Ratio = ratio;
                    EC = ec1;
                    ED = ed1;
                    t1 = T1;
                end
                break;
            end
        end
        
        if t1 == 0
            fprintf('Max Payload : %d\n', Payload-1000);
            break;
        end
        
        ec = 0;
        ed = 0;
        for i = 1:(A-2)-a
            for j = 1:(B-2)-b
                if NL(i,j) < t1
                    ec = ec + EC1(i,j);
                    ed = ed + ED1(i,j);
                end
                if ec >= Payload
                    break;
                end
            end
        end
        
        PSNR = 10*log10(A*B*255^2 / ed);
        cnt = cnt + 1;
        R(:,cnt) = [Payload, PSNR];
    end
    res = R(:,1:cnt);
    save(istr, 'res');
end