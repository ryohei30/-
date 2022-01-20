%function Yd = tsvd_dct_version(Y)
% t-SVDによるガウシアンテンソルデノイズ
% 疑似乱数のシード
%rng(1)

% パラメータ
%   原テンソルのパラメータ
%N = 256;
%M = 256;
%K = 3;

% 観測テンソルのパラメータ
%noise_sigma = 0.2;

% t-SVDのパラメータ
%lambda = 4.4;

%最適なlambdaを探す場合
%lambda_list = [];
%psnr_list = [];
%for i=0:0.1:10
    %lambda = 0.1+i;
    %lambda_list(end+1) = lambda;


% 適当な低ランクテンソルを作る
%a = randn(N,1,1);
%b = randn(1,M,1);
%c = randn(1,1,K);
%X = a.*b.*c;

%画像を使う場合
%pic_name = './Mandrill.bmp';
%X = double(imread(pic_name));
%X = X/255;


%Y = X + randn(size(X))*noise_sigma;


function [SURE_dct,psnr_dct,MSE_dct] = tsvd_dct_version(Y,K,noise_sigma,X)

    
% 第３軸方向にDCT
Yb = dct(Y,[],3); %dctしたテンソルをYbとする
Ybd = zeros(size(Yb)); %Ybの同じサイズの空のテンソルを定義（ここに処理後の値を入れる）
SURE_list = [];
lambda = 4.4;

% 特異値分解してソフト閾値処理
for i = 1:K
    Ybd(:,:,i) = SVDandThresh(Yb(:,:,i), lambda);
    %SUREの計算
    [~, SURE] = calcSURE(Ybd(:,:,i), lambda, noise_sigma);
    SURE_list(end+1) = SURE;
end
Yd = idct(Ybd,[],3); %閾値処理したテンソルを逆dctしたものをYdとする


% MSEをはかる
MSE_dct = mean(sum((X - Yd).^2),'all'); 

%PSNR

psnr_dct = PSNR(X,Yd,1.0);
SURE_dct = sum(SURE_list);

%psnr_list(end+1) = PSNR(X,Yd,1.0);
%end


% 結果表示
%fprintf("result = %f\n", result) %MSEの結果
%fprintf("psnr = %f\n", psnr)
%result2 = max(psnr_list);
%fprintf("result2 = %f\n", result2)
%max(lambda_list)


%プロット
%figure(1)
%plot(lambda_list,psnr_list)
%figure(2)
%subplot(1,3,1)
%imshow(X)
%subplot(1,3,2)
%imshow(Y)
%subplot(1,5,3)
%imshow(Xhat)
%subplot(1,5,4)
%imshow(Xhat2)
%subplot(1,3,3)
%imshow(Yd)



       
function X = SVDandThresh(Y, lambda)
    [U,S,V] = svd(Y); %特異値分解を実行
    S=max(S-lambda, 0); %Sを閾値処理し、Sを再定義（閾値はlambda)
    X = U*S*V'; %閾値処理した特異値でテンソルXを定義
end
function psnr = PSNR(Xfull,Xrecover,maxP)

% Written by Canyi Lu (canyilu@gmail.com)

Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
[n1,n2,n3] = size(Xrecover);
MSE = norm(Xfull(:)-Xrecover(:))^2/(n1*n2*n3);
psnr = 10*log10(maxP^2/MSE);
end

function [C_df, SURE] = calcSURE(Y, lambda, nsig)
    [U,SigmaY,V]=svd(Y);
    sY  = diag(SigmaY);
    [n, m] = size(Y); 
    lambda;
    minnm = min(n,m);
    %idx_w = (sY-lambda >0);%しきい値処理で0にされない範囲（f_iの微分値があるのはここだけ）
    term1 = abs(m-n) * sum(max((1- lambda./sY),0)); %Eq. (9)第1項
    term2 = sum(sY > lambda); %Eq. (9)第2項 びっくり
    numer = repmat(sY.*max((sY-lambda),0), [1,minnm]); %Eq. (9)第3項の分子
    denom = repmat(sY.^2, [1,minnm]) - repmat((sY.^2)', [minnm,1]); %Eq. (9)第3項の分母
    DD = numer./denom; % i=jのときNaNが出るが…
    DD( isnan(DD) | isinf(DD) | abs(DD) > 1E6 )  =   0;%ここで消すので問題なし←割と邪道なのでちゃんとやったほうがいい
    term3 = 2*sum(DD(:)); %Eq. (5)第3項
    C_df = term1+term2+term3;
    SURE = - m*n*nsig^2 + sum(min(lambda^2,sY.^2)) + 2*nsig^2*C_df;
end
end