clear all; clc;
addpath(genpath('Origin Images')); addpath(genpath('result')); addpath(genpath('tools'));

mask2 = [false, true; true, true];
mask3 = [false, true, true; true, true, true; true, true, true];
mask4 = [false, true, true, true; true, true, true, true; true, true, true, true; true, true, true, true];
mask5 = [false, true, true, true, true; true, true, true, true, true; true, true, true, true, true; true, true, true, true, true; true, true, true, true, true];

Imgs = {'Lena', 'Baboon', 'Airplane', 'Barbara', 'Lake', 'Peppers', 'Boat', 'Elaine'};
%[28000, ]
%%
for tt = 1:1:8
    Iname = Imgs{tt}; 
    istr = ['5x5/Proposed_2019_',Iname,'.mat']
    I = double(imread([Iname,'.bmp'])); 
    [A, B] = size(I);

    %%
    EdgInfo = 40;
    a = 5; b = 5;
    R = zeros(2,1000);
    cnt = 0;
    Tmax = 3000;
    for Payload =5000+EdgInfo:1000:100000+EdgInfo
        tic
%         Payload
        %%
        H1 = cell(1,Tmax);
        H2 = cell(1,Tmax);
        H3 = cell(1,Tmax);
        H4 = cell(1,Tmax);
        for i = 1 : Tmax
            H1{i} = zeros(1,256);
            H2{i} = zeros(1,256);
            H3{i} = zeros(1,256);
            H4{i} = zeros(1,256);
        end

        EC1 = zeros(A-4,B-4);
        EC2 = zeros(A-4,B-4);
        EC3 = zeros(A-4,B-4);
        EC4 = zeros(A-4,B-4);
        ED1 = zeros(A-4,B-4);
        ED2 = zeros(A-4,B-4); 
        ED3 = zeros(A-4,B-4);
        ED4 = zeros(A-4,B-4);
        NL = zeros(A-4,B-4);
        for i = 1:(A-4)
            for j = 1:(B-4)
                for ii = 1:4
                    for jj = 1:5
                        if ii == 2 || ii == 3 || ii == 4 || jj == 2 || jj == 3 || jj == 4 || jj == 5
%                             [1*(i-1)+ii,1*(j-1)+jj, 1*(i-1)+ii+1,1*(j-1)+jj]
                            NL(i,j) = NL(i,j) + abs(I(1*(i-1)+ii,1*(j-1)+jj) - I(1*(i-1)+ii+1,1*(j-1)+jj));
                        end
                    end
                end
                for ii = 1:5
                    for jj = 1:4
                        if ii == 2 || ii == 3 || ii == 4 || ii == 5 || jj == 2 || jj == 3 || jj == 4
%                             [1*(i-1)+ii,1*(j-1)+jj, 1*(i-1)+ii,1*(j-1)+jj+1]
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

                % 3x3 block
                X = I(i:i+2,j:j+2);
                X = X(mask3);
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
                
                % 4x4 block
                X = I(i:i+3,j:j+3);
                X = X(mask4);
                [Y, In] = sort(X);
                if Y(end) ~= Y(1)
                    % dmax
                    if I(i,j) >= Y(end)
                        dmax = I(i,j) - Y(end); % >= 0
                        H3{T+1}(dmax+1) = H3{T+1}(dmax+1) + 1;
                        if dmax == 0
                            EC3(i,j) = EC3(i,j) + 1;
                            ED3(i,j) = ED3(i,j) + 0.5;
                        else
                            ED3(i,j) = ED3(i,j) + 1;
                        end
                    end
                    % dmin
                    if I(i,j) <= Y(1)
                        dmin = Y(1) - I(i,j); % >= 0
                        H3{T+1}(dmin+1) = H3{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC3(i,j) = EC3(i,j) + 1;
                            ED3(i,j) = ED3(i,j) + 0.5;
                        else
                            ED3(i,j) = ED3(i,j) + 1;
                        end
                    end
                else % In(end) == In(1)
                    if I(i,j) <= Y(end)
                        dmin = Y(1) - I(i,j); % >= 0
                        H3{T+1}(dmin+1) = H3{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC3(i,j) = EC3(i,j) + 1;
                            ED3(i,j) = ED3(i,j) + 0.5;
                        else
                            ED3(i,j) = ED3(i,j) + 1;
                        end
                    end
                end
                
                    % 5x5 block
                X = I(i:i+4,j:j+4);
                X = X(mask5);
                [Y, In] = sort(X);
                if Y(end) ~= Y(1)
                    % dmax
                    if I(i,j) >= Y(end)
                        dmax = I(i,j) - Y(end); % >= 0
                        H4{T+1}(dmax+1) = H4{T+1}(dmax+1) + 1;
                        if dmax == 0
                            EC4(i,j) = EC4(i,j) + 1;
                            ED4(i,j) = ED4(i,j) + 0.5;
                        else
                            ED4(i,j) = ED4(i,j) + 1;
                        end
                    end
                    % dmin
                    if I(i,j) <= Y(1)
                        dmin = Y(1) - I(i,j); % >= 0
                        H4{T+1}(dmin+1) = H4{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC4(i,j) = EC4(i,j) + 1;
                            ED4(i,j) = ED4(i,j) + 0.5;
                        else
                            ED4(i,j) = ED4(i,j) + 1;
                        end
                    end
                else % In(end) == In(1)
                    if I(i,j) <= Y(end)
                        dmin = Y(1) - I(i,j); % >= 0
                        H4{T+1}(dmin+1) = H4{T+1}(dmin+1) + 1;
                        if dmin == 0
                            EC4(i,j) = EC4(i,j) + 1;
                            ED4(i,j) = ED4(i,j) + 0.5;
                        else
                            ED4(i,j) = ED4(i,j) + 1;
                        end
                    end
                end
                
            end
        end

        for i = 2 : 1 : Tmax
            H1{i} = H1{i} + H1{i-1};
            H2{i} = H2{i} + H2{i-1};
            H3{i} = H3{i} + H3{i-1};
            H4{i} = H4{i} + H4{i-1};
        end
        ECs1 = zeros(1,Tmax);
        ECs2 = zeros(1,Tmax);
        ECs3 = zeros(1,Tmax);
        ECs4 = zeros(1,Tmax);
        EDs1 = zeros(1,Tmax);
        EDs2 = zeros(1,Tmax);
        EDs3 = zeros(1,Tmax);
        EDs4 = zeros(1,Tmax);
        for i = 1 : 1 : Tmax
            ECs1(i) = H1{i}(1);
            ECs2(i) = H2{i}(1);
            ECs3(i) = H3{i}(1);
            ECs4(i) = H4{i}(1);
            
            EDs1(i) = 0.5*H1{i}(1) + sum(H1{i}(2:end));
            EDs2(i) = 0.5*H2{i}(1) + sum(H2{i}(2:end));
            EDs3(i) = 0.5*H3{i}(1) + sum(H3{i}(2:end));
            EDs4(i) = 0.5*H4{i}(1) + sum(H4{i}(2:end));
        end
        
        Ratio = 100000;
        nls = unique(NL(:));
        t1 = 0; t2 = 0; t3 = 0; t4 = 0;
        for i = 1 : numel(nls)
            T1 = nls(i) + 1;
            ec1 = ECs1(T1);
            if ec1 > Payload
                break;
            end
            for j = i+1 : numel(nls)
                T2 = nls(j) + 1;
                ec2 = ECs2(T2)-ECs2(T1);
                if ec1 + ec2 > Payload
                    break;
                end
                for l = j+1 : numel(nls)
                    T3 = nls(l) + 1;
                    ec3 = ECs3(T3)-ECs3(T2);
                    if ec1 + ec2 + ec3 > Payload
                        break;
                    end
                    
                    h4 = zeros(1,256);
                    for k = l+1 : numel(nls)
                        T4 = nls(k) + 1;
                        ec4 = ECs4(T4)-ECs4(T3);

                        if ec1 + ec2 + ec3 + ec4 > Payload
                            ec = ec1 + ec2 + ec3 + ec4;
                            ed = EDs1(T1) + EDs2(T2)-EDs2(T1) + EDs3(T3)-EDs3(T2) + EDs4(T4)-EDs4(T3);
                            ratio = ed/ec;
                            if Ratio > ratio
                                Ratio = ratio;
                                EC = ec;
                                ED = ed;
                                t1 = T1;
                                t2 = T2;
                                t3 = T3;
                                t4 = T4;
                            end
                            break;
                        end
                    end
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
        for i = 1:(A-4)
            if flag == 1
                break;
            end
            for j = 1:(B-4)
                if NL(i,j) < t1
                    ec = ec + EC1(i,j);
                    ed = ed + ED1(i,j);
                end
                if NL(i,j) >= t1 && NL(i,j) < t2
                    ec = ec + EC2(i,j);
                    ed = ed + ED2(i,j);
                end
                if NL(i,j) >= t2 && NL(i,j) < t3
                    ec = ec + EC3(i,j);
                    ed = ed + ED3(i,j);
                end
                if NL(i,j) >= t3 && NL(i,j) < t4
                    ec = ec + EC4(i,j);
                    ed = ed + ED4(i,j);
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
        [Payload, PSNR, toc]
    end
    res = R(:,1:cnt);
    save(istr, 'res');
end