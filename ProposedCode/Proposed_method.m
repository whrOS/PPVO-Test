tic
% clear all;
clc

addpath(genpath('Origin Images'));
addpath(genpath('pairwise_IPVO'));
addpath(genpath('tools'));
addpath(genpath('result'));
Imgs = {'Lena', 'Baboon', 'Airplane', 'Lake', 'Peppers', 'Boat', 'Barbara', 'Elaine'};

%%
% n1xn2      4  bits
% T          8  bits
% kend      12  bits
% LM        12  bits
% In total: 36  bits
% Optimal Map Index : 11 bits
edge_info = 36 + 11;
maps = {};
map = cell(3);
mapECs = {};
mapEDs = {};

[mapEDs, mapECs, maps] = getOneMap(mapEDs, mapECs, maps,map,[0,0]);
mapNum = numel(maps);

%%
for tt = 1:8
    Iname = Imgs{tt};
    
    istr = ['Proposed_2019_',Iname,'.mat']
    I = double(imread([Iname,'.bmp']));
    [A, B] = size(I);
    %%
    S = zeros(8,100);
    for payload = 5000 + edge_info : 1000 : 100000 + edge_info
        payload - edge_info
        Tmax = 500;
        R = zeros(8,1);
        rc = 1;
        for a = 2:5
            for b = 2:5
                Hs = cell(1,Tmax); % T = 0 - Tmax
                for i = 1:Tmax
                    Hs{i} = zeros(512);
                end
                NL = zeros(floor((A-2)/a),floor((B-2)/b));
                for i = 1:floor((A-2)/a)
                    for j = 1:floor((B-2)/b)
                        for ii = 1:a+1
                            for jj = 1:b+2
                                if ii == a+1 || jj == b+1 || jj == b+2
                                    NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii+1,b*(j-1)+jj));
                                end
                            end
                        end
                        for ii = 1:a+2
                            for jj = 1:b+1
                                if ii == a+1 || ii == a+2 || jj == b+1
                                    NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii,b*(j-1)+jj+1));
                                end
                            end
                        end
                        if NL(i,j) < Tmax
                            X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
                            X = X(:);
                            [Y, In] = sort(X);
                            % max
                            if In(a*b) < In(a*b-1)
                                dmax = Y(a*b-1) - Y(a*b);
                            else
                                dmax = Y(a*b) - Y(a*b-1);
                            end
                            % min
                            if In(2) < In(1)
                                dmin = Y(1) - Y(2);
                            else
                                dmin = Y(2) - Y(1);
                            end
                            Hs{NL(i,j)+1}(256-dmax, dmin+256) = Hs{NL(i,j)+1}(256-dmax, dmin+256) + 1;
                        end
                    end
                end
                
                for i = 2:Tmax
                    Hs{i} = Hs{i} + Hs{i-1};
                end
                
                ec0 = zeros(1,Tmax);
                ed0 = zeros(1,Tmax);
                optMapEC = zeros(512,512); optMapEC(256:256+1,:) = 1;  optMapEC(:,256-1:256) = 1;  optMapEC(256-2:256+3,256-3:256+2) = 0;
                optMapED = ones(512,512)*2;optMapED(256:256+1,:) = 3/2;optMapED(:,256-1:256) = 3/2;optMapED(256-2:256+3,256-3:256+2) = 0;
                for T = 1:Tmax
                    ec0(T) = sum(sum(Hs{T} .* optMapEC));
                    ed0(T) = sum(sum(Hs{T} .* optMapED));
                end
                
                %%
                optMapEC = zeros(6,6);
                optMapED = zeros(6,6);
                for m = 1 : mapNum
                    ec = mapECs(m);
                    ed = mapEDs(m);
                    ec = ec{1}; ec1 = rot90(ec'); ec2 = fliplr(rot90(ec')); ec3 = rot90(rot90(ec'),2); ec4 = flipud(rot90(ec'));
                    ed = ed{1}; ed1 = rot90(ed'); ed2 = fliplr(rot90(ed')); ed3 = rot90(rot90(ed'),2); ed4 = flipud(rot90(ed'));
                    optMapEC = [ec2 ec1; ec3 ec4];
                    optMapED = [ed2 ed1; ed3 ed4];
                    for T = 1:Tmax
                        H = Hs{T}(256-2:256+3, 256-3:256+2);
                        ECMT = sum(sum(H .* optMapEC)) + ec0(T);
                        if ECMT >= payload
                            EDMT = sum(sum(H .* optMapED)) + ed0(T);
                            R(:,rc) = [payload,ECMT,EDMT,EDMT/ECMT,T,a,b,m]'; % T -> T+1
                            %                             EDMT = 10*log10(A*B*255^2/EDMT);
                            rc = rc+1;
                            break;
                        end
                    end
                end
            end
        end
        if sum(R) == 0
            break;
        end
        %%
        [~,idx] = min(R(4,:));
        M = R(end,idx);
        Topt = R(5,idx);
        a = R(6,idx);
        b = R(7,idx);

        ec = mapECs(M);
        ed = mapEDs(M);
        ec = ec{1}; ec1 = rot90(ec'); ec2 = fliplr(rot90(ec')); ec3 = rot90(rot90(ec'),2); ec4 = flipud(rot90(ec'));
        ed = ed{1}; ed1 = rot90(ed'); ed2 = fliplr(rot90(ed')); ed3 = rot90(rot90(ed'),2); ed4 = flipud(rot90(ed'));
        optMapEC = zeros(512,512);
        optMapEC(256:256+1,:) = 1;
        optMapEC(:,256-1:256) = 1;
        optMapEC(256-2:256+3,256-3:256+2) = [ec2 ec1;
            ec3 ec4];
        optMapED = ones(512,512)*2;
        optMapED(256:256+1,:) = 3/2;
        optMapED(:,256-1:256) = 3/2;
        optMapED(256-2:256+3,256-3:256+2) = [ed2 ed1;
            ed3 ed4];
        
        NL = zeros(floor((A-2)/a),floor((B-2)/b));
        EC = zeros(floor((A-2)/a),floor((B-2)/b));
        ED = zeros(floor((A-2)/a),floor((B-2)/b));
        for i = 1:floor((A-2)/a)
            for j = 1:floor((B-2)/b)
                for ii = 1:a+1
                    for jj = 1:b+2
                        if ii == a+1 || jj == b+1 || jj == b+2
                            NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii+1,b*(j-1)+jj));
                        end
                    end
                end
                for ii = 1:a+2
                    for jj = 1:b+1
                        if ii == a+1 || ii == a+2 || jj == b+1
                            NL(i,j) = NL(i,j) + abs(I(a*(i-1)+ii,b*(j-1)+jj) - I(a*(i-1)+ii,b*(j-1)+jj+1));
                        end
                    end
                end
                X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
                X = X(:);
                [Y, In] = sort(X);
                % max
                if In(a*b) < In(a*b-1)
                    dmax = Y(a*b-1) - Y(a*b);
                else
                    dmax = Y(a*b) - Y(a*b-1);
                end
                % min
                if In(2) < In(1)
                    dmin = Y(1) - Y(2);
                else
                    dmin = Y(2) - Y(1);
                end
                EC(i,j) = optMapEC(256-dmax,256+dmin);
                ED(i,j) = optMapED(256-dmax,256+dmin);
            end
        end
        
        optEC = 0;
        optED = 0;
        for i = 1:floor((A-2)/a)
            for j = 1:floor((B-2)/b)
                if optEC < payload && NL(i,j) < Topt
                    optEC = optEC+EC(i,j);
                    optED = optED+ED(i,j);
                    if optEC >= payload
                        Tend = [i,j];
                    end
                end
            end
        end
        %%
        m_index = (payload-edge_info-5000)/1000 + 1;
        R(3,idx) = 10*log10(A*B*255^2 / optED);
        S(:,m_index) = R(:,idx);
    end
    S = S(:,1:m_index);
    res = [S(1,:);S(3:end,:)];
        save(istr, 'res');
end
