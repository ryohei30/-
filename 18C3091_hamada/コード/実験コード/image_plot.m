clear all;
close all;
format long

% 疑似乱数のシー�?
seed = 1;
rng(seed, 'twister')

% 観測�?ンソルのパラメータ
noise_sigma = 0.5;

X = double(imread('Balloon.bmp'))/255;
Y = X + randn(size(X))*noise_sigma;

%lambdaラインサーチ�?��?割数
%N=40;
%ラインサーチするlambdaの上限と下限
%{
lambda_min = 10;
lambda_max = 15;

lambda_list = linspace(lambda_min,lambda_max,N);
MSE_dct_list = zeros(size(lambda_list));
MSE_dft_list = zeros(size(lambda_list));
MSE_pca_list = zeros(size(lambda_list));
SURE_dct_list = zeros(size(lambda_list));
SURE_dft_list = zeros(size(lambda_list));
SURE_pca_list = zeros(size(lambda_list));
PSNR_dft_list = zeros(size(lambda_list));
PSNR_dct_list = zeros(size(lambda_list));
PSNR_pca_list = zeros(size(lambda_list));
%}

%for i =  1:N
    %lambda = lambda_list(i);
    lambda = 13.58974359;
    
    %---DCT
    SURE_temp_list = zeros(3,1);
    Yb = dct(Y,[],3); %dctした�?ンソルをYbとする
    Ybd = zeros(size(Yb)); %Ybの同じサイズの空の�?ンソルを定義?��ここに処�?後�?�値を�?�れる?�?
    % 特異値�?解してソフト閾値処�?
    for j = 1:3
        Ybd(:,:,j) = SVDandThresh(Yb(:,:,j), lambda);
        %SUREの計�?
        [~, SURE] = calcSURE(Yb(:,:,j), lambda, noise_sigma);
        SURE_temp_list(j) = SURE/numel(Yb(:,:,j));
    end
    dct_Yd = idct(Ybd,[],3); %閾値処�?した�?ンソルを�??dctしたも�?�をYdとする

    %MSE_dct_list(i)  = mean((X - Yd).^2, 'all');
    %SURE_dct_list(i) = mean(SURE_temp_list,'all');
    %PSNR_dct_list(i) = PSNR(X,Yd,1.0);

    %---DFT
    SURE_temp_list = zeros(3,1);
    Yb = fft(Y,[],3)/sqrt(3); %dftした�?ンソルをYbとする
    Ybd = zeros(size(Yb)); %Ybの同じサイズの空の�?ンソルを定義?��ここに処�?後�?�値を�?�れる?�?
    % 特異値�?解してソフト閾値処�?
    for j = 1:3
        Ybd(:,:,j) = SVDandThresh(Yb(:,:,j), lambda);
        %SUREの計�?
        [~, SURE] = calcSURE(Yb(:,:,j), lambda, noise_sigma);
        SURE_temp_list(j) = SURE/numel(Yb(:,:,j));
    end
    dft_Yd = ifft(Ybd,[],3)*sqrt(3); %閾値処�?した�?ンソルを�??dftしたも�?�をYdとする

    %MSE_dft_list(i)  = mean((X - Yd).^2, 'all');
    %SURE_dft_list(i) = mean(SURE_temp_list,'all');
    %PSNR_dft_list(i) = PSNR(X,Yd,1.0);

    %---PCA
    % --- PCA前準備
    [n1,n2,n3] = size(Y);
    UfY = Unfold(Y, [n1,n2,n3], 3);%観測のunfold
    cUfY = UfY-mean(UfY,2);
    C = cUfY*cUfY'/(n1*n2);
    [Eig,~] = eig(C); %固有�?�クトルを取�?
    % ---  PCA t-SVD �?ノイズ
    SURE_temp_list = zeros(3,1);
    
    PY = Eig'*UfY; %固有�?�クトル行�?�かける�?ンソルunfold
    Yb = Fold(PY,[n1,n2,n3],3); %PCAした�?ンソルをYbとする

    Ybd = zeros(size(Yb)); %Ybの同じサイズの空の�?ンソルを定義?��ここに処�?後�?�値を�?�れる?�?
    % 特異値�?解してソフト閾値処�?
    for j = 1:3
        Ybd(:,:,j) = SVDandThresh(Yb(:,:,j), lambda);
        %SUREの計�?
        [~, SURE] = calcSURE(Yb(:,:,j), lambda, noise_sigma);
        SURE_temp_list(j) = SURE/numel(Yb(:,:,j));
    end

    UfX = Unfold(Ybd, [n1,n2,n3], 3); %閾値処�?した�?ンソルをunfoldしたも�?�をYdとする
    iPX = Eig*UfX; %unfoldしたYbdに固有�?�クトル行�?�をかけてる�?pcaしてる？�?
    pca_Yd = Fold(iPX, [n1,n2,n3], 3); %�?ンソルに戻�?
    %MSE_pca_list(i)  = mean((X - Yd).^2, 'all');
    %SURE_pca_list(i) = mean(SURE_temp_list,'all');
    %PSNR_pca_list(i) = PSNR(X,Yd,1.0);
    
    
%end

%{
plot(lambda_list, SURE_dct_list,'b-.')
hold on;
plot(lambda_list, MSE_dct_list, 'b-')
plot(lambda_list, SURE_dft_list, 'r-.')
plot(lambda_list, MSE_dft_list, 'r-')
plot(lambda_list, SURE_pca_list, 'g-.')
plot(lambda_list, MSE_pca_list, 'g-')
legend('SURE(dct)','MSE(dct)','SURE(dft)','MSE(dft)','SURE(pca)','MSE(pca)')
%}
figure(2)
subplot(1,5,1)
imshow(dft_Yd)
subplot(1,5,2)
imshow(dct_Yd)
subplot(1,5,3)
imshow(pca_Yd)
subplot(1,5,4)
imshow(X)
subplot(1,5,5)
imshow(Y)



% SVTによる行�?�デノイズ

function X = SVDandThresh(Y, lambda)
    [U,S,V] = svd(Y); %特異値�?解を実�?
    S=max(S-lambda, 0); %Sを閾値処�?し�?�Sを�?�定義?��閾値はlambda)
    X = U*S*V'; %閾値処�?した特異値で�?ンソルXを定義
end

%SURE?��行�?�SVT専用?�?
function [C_df, SURE] = calcSURE(Y, lambda, nsig)
    [U,SigmaY,V]=svd(Y);
    sY  = diag(SigmaY);
    [n, m] = size(Y); 
    lambda;
    minnm = min(n,m);
    %idx_w = (sY-lambda >0);%しき�?値処�?で0にされな�?�?囲?�?f_iの微�?値がある�?�はここ�?け�?
    term1 = abs(m-n) * sum(max((1- lambda./sY),0)); %Eq.(5)第1�?
    term2 = sum(sY > lambda); %Eq.(5)第2�? 
    numer = repmat(sY.*max((sY-lambda),0), [1,minnm]); %Eq.(5)第3�?の�?�?
    denom = repmat(sY.^2, [1,minnm]) - repmat((sY.^2)', [minnm,1]); %Eq.(5)第3�?の�?�?
    DD = numer./denom; % i=jのときNaNが�?�るが…
    DD( isnan(DD) | isinf(DD) | abs(DD) > 1E6 )  =   0;%ここで消すので問題なし�?�割と邪道なのでち�?んと�?ったほ�?がい�?
    term3 = 2*sum(DD(:)); %Eq. (5)第3�?
    C_df = term1+term2+term3;
    SURE = - m*n*nsig^2 + sum(min(lambda^2,sY.^2)) + 2*nsig^2*C_df;
end
function psnr = PSNR(Xfull,Xrecover,maxP)

% Written by Canyi Lu (canyilu@gmail.com)

Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
[n1,n2,n3] = size(Xrecover);
MSE = norm(Xfull(:)-Xrecover(:))^2/(n1*n2*n3);
psnr = 10*log10(maxP^2/MSE);
end
function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);%shiftdimで三次�?�?1次�?方向に持ってこられて�?る（多�??�?,reshapeで?��次�?にして�?る�??3�?65536にして�?ると思われる�?
end

function [X] = Fold(X, dim, i)
dim = circshift(dim, [1-i, 1-i]);
X = shiftdim(reshape(X, dim), length(dim)+1-i);
end