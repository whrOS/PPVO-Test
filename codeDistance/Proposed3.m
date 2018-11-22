clear all; clc;
addpath(genpath('Origin Images')); addpath(genpath('result')); addpath(genpath('tools'));

mask1 = [false, true; true, true];
mask2 = [false, false, true; true, true, true];
mask3 = [false, true, true; true, true, true; true, true, true];
mask4 = [false, false, true, true; true, true, true, true; true, true, true, true];
mask5 = [false, false, false, true, true; true, true, true, true, true; true, true, true, true, true];

Imgs = {'Lena', 'Baboon', 'Airplane', 'Barbara', 'Lake', 'Peppers', 'Boat', 'Elaine'};
%%
for tt = 1:1:8
    Iname = Imgs{tt}; 
    istr = ['Proposed_2019_',Iname,'.mat']
    I = double(imread([Iname,'.bmp'])); 
    [A, B] = size(I);

    %%
    EdgInfo = 40;
    a = 3; b = 3;
    R = zeros(2,1000);
    cnt = 0;
    for Payload = 20000+EdgInfo:1000:100000+EdgInfo
%         Payload
        %%
        H1 = cell(1,2048);
        H2 = cell(1,2048);
        for i = 1 : 2048
            H1{i} = zeros(1,256);
            H2{i} = zeros(1,256);
        end

        EC1 = zeros(A-2,B-2);
        EC2 = zeros(A-2,B-2);
        ED1 = zeros(A-2,B-2);
        ED2 = zeros(A-2,B-2); 
        NL = zeros(A-2,B-2);
        for i = 1:(A-2)
            for j = 1:(B-2)
                flag1 = 0;
                flag2 = 0;
                for ii = 1:2
                    for jj = -1:3
                        if ii == 2 || jj == 2 || jj == 3
                            if 1*(j-1)+jj == 0
                                flag1 = 1;
                            elseif 1*(j-1)+jj == -1
                                flag2 = 1;
                            else
%                                 [1*(i-1)+ii,1*(j-1)+jj, 1*(i-1)+ii+1,1*(j-1)+jj]
                                NL(i,j) = NL(i,j) + abs(I(1*(i-1)+ii,1*(j-1)+jj) - I(1*(i-1)+ii+1,1*(j-1)+jj));
                            end
                        end
                    end
                end
                for ii = 1:3
                    for jj = -1:2
                        if ii == 2 || ii == 3 || jj == 2
                            if 1*(j-1)+jj <=0
                            else
%                                 [1*(i-1)+ii,1*(j-1)+jj, 1*(i-1)+ii,1*(j-1)+jj+1]
                                NL(i,j) = NL(i,j) + abs(I(1*(i-1)+ii,1*(j-1)+jj) - I(1*(i-1)+ii,1*(j-1)+jj+1));
                            end
                        end
                    end
                end
                if flag1 == 1 && flag2 == 1
%                     flag
                    NL(i,j) = floor(NL(i,j) * 10/6);
                end
                if flag1 == 1 && flag2 == 0
%                     flag
                    NL(i,j) = floor(NL(i,j) * 10/8);
                end
                T = NL(i,j);

                % 2x2 block
                if flag1 == 1
                    X = I(i:i+1,j:j+1);
                    X = X(mask1);
                else
                    X = I(i:i+1,j-1:j+1);
                    X = X(mask2);
                end
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

                % 3x3 block
                if flag1 == 1 && flag2 == 1
                    X = I(i:i+2,j:j+2);
                    X = X(mask3);
                end
                if flag1 == 1 && flag2 == 0
                    X = I(i:i+2,j-1:j+2);
                    X = X(mask4);
                end
                if flag1 == 0 && flag2 == 0
                    X = I(i:i+2,j-2:j+2);
                    X = X(mask5);
                end
                [Y, In] = sort(X);
                if Y(end) ~= Y(1)
                    % dmax
                    if I(i,j) >= Y(end)
                        dmax = I(i,j) - Y(end); % >= 0
                        H2{T+1}(dmax+1) = H2{T+1}(dmax+1) + 1;
                        if dmax == 0
                            EC2(i,j) = EC2(i,j) + 1;
                            ED2(i,j) = ED2(i,j) + 0.5;
                        else
                            ED2(i,j) = ED2(i,j) + 1;
                        end
                    end
                    % dmin
                    if I(i,j) <= Y(1)
                        dmin = Y(1) - I(i,j); % >= 0
                        H2{T+1}(dmin+1) = H2{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC2(i,j) = EC2(i,j) + 1;
                            ED2(i,j) = ED2(i,j) + 0.5;
                        else
                            ED2(i,j) = ED2(i,j) + 1;
                        end
                    end
                else % In(end) == In(1)
                    if I(i,j) <= Y(end)
                        dmin = Y(1) - I(i,j); % >= 0
                        H2{T+1}(dmin+1) = H2{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC2(i,j) = EC2(i,j) + 1;
                            ED2(i,j) = ED2(i,j) + 0.5;
                        else
                            ED2(i,j) = ED2(i,j) + 1;
                        end
                    end
                end
            end
        end

        Ratio = 100000;
        nls = unique(NL(:));
        t1 = 0; t2 = 0;
        for i = 1 : numel(nls)
            T1 = nls(i) + 1;
            h1 = zeros(1,256);
            for k = 1 : T1
                h1 = h1 + H1{k};
            end
            ec1 = h1(1);
            if ec1 > Payload
                break;
            end
            h2 = zeros(1,256);
            for j = i+1 : numel(nls)
                T2 = nls(j) + 1;
                h2 = h2 + H2{T2};
                ec2 = h2(1);
                if ec1 + ec2 > Payload
                    ed1 = 0.5*h1(1) + sum(h1(2:end));
                    ed2 = 0.5*h2(1) + sum(h2(2:end));
                    ec = ec1 + ec2;
                    ed = ed1 + ed2;
                    ratio = ed/ec;
                    if Ratio > ratio
                        Ratio = ratio;
                        EC = ec;
                        ED = ed;
                        t1 = T1;
                        t2 = T2;
                    end
                    break;
                end
            end
        end

        if t1 == 0
            fprintf('Max Payload : %d\n', Payload-1000);
            break;
        end

        ec = 0;
        ed = 0;
        flag = 0;
        for i = 1:(A-2)
            if flag == 1
                break;
            end
            for j = 1:(B-2)
                if NL(i,j) < t1
                    ec = ec + EC1(i,j);
                    ed = ed + ED1(i,j);
                end
                if NL(i,j) < t2 && NL(i,j) >= t1
                    ec = ec + EC2(i,j);
                    ed = ed + ED2(i,j);
                end
                if ec >= Payload
                    kend = [i,j];
                    flag = 1;
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