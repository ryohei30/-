%function I = tsvd_pca_version(image)
% t-SVDによるガウシアンテンソルデノイズ
% 疑似乱数のシード
%rng(1)

% パラメータ
%  原テンソルのパラメータ
%N = 256;
%M = 256;
%K = 3;

% 観測テンソルのパラメータ
%noise_sigma = 0.2;

%   t-SVDのパラメータ
%lambda = 0.1;

%最適なlambdaを探す場合
%lambda_list = [];
%psnr_list = [];
%for i=0:0.1:5
%    lambda = 0.1+i;
%    lambda_list(end+1) = lambda;


% 適当な低ランクテンソルを作る
%a = randn(N,1,1);
%b = randn(1,M,1);
%c = randn(1,1,K);
%X = a.*b.*c;

%画像を使う場合
%pic_name = './Mandrill.bmp';
%pic_name = image;
%X = double(imread(pic_name));
%X = image;
%X = X/255;

%Y = X + randn(size(X))*noise_sigma;

function [SigmaY,SURE_pca,psnr_pca,MSE_pca] = tsvd_pca_version(Y,K,noise_sigma,X)
lambda =4.475;

%第３軸方向にPCA
[n1,n2,n3] = size(Y);
UfY = Unfold(Y, [n1,n2,n3], 3);%観測のunfold
cUfY = UfY-mean(UfY,2);
C = cUfY*cUfY'/(n1*n2);
%UfXh = Unfold(X, [n1,n2,n3], 3);%原テンソルのunfold
%cUfXh = UfXh-mean(UfXh,2); %中心化
%C = cUfXh*cUfXh'/(n1*n2); %共分散行列を獲得
[Eig,~] = eig(C); %固有ベクトルを取得
PY = Eig'*UfY; %固有ベクトル行列かけるテンソルunfold
Yb = Fold(PY,[n1,n2,n3],3); %PCAしたテンソルをYbとする
Ybd = zeros(size(Yb)); %Ybの同じサイズの空のテンソルを定義（ここに処理後の値を入れる）

SURE_list = [];
    
% 特異値分解してソフト閾値処理
for i = 1:K
    Ybd(:,:,i) = SVDandThresh(Yb(:,:,i), lambda);
    %SUREの計算
    [SigmaY,~, SURE] = calcSURE(Ybd(:,:,i), lambda, noise_sigma);
    SURE_list(end+1) = SURE;
end
 UfX = Unfold(Ybd, [n1,n2,n3], 3); %閾値処理したテンソルをunfoldしたものをYdとする
 iPX = Eig*UfX; %unfoldしたYbdに固有ベクトル行列をかけてる（pcaしてる？）
 
 Yd = Fold(iPX, [n1,n2,n3], 3); %テンソルに戻す

% MSEをはかる
MSE_pca = mean(sum((X - Yd).^2),'all'); 

%PSNRとSUREを算出
%psnr_list(end+1) = PSNR(X,Yd,1.0);
psnr_pca = PSNR(X,Yd,1.0);
SURE_pca = sum(SURE_list);
%end

%返す値
%[M,I]=max(psnr_list);

% 結果表示
%fprintf("MSE = %f\n", MSE) %MSEの結果
%fprintf("psnr = %f\n", psnr) %PSNRの結果

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
    X = U*S*V'; %閾値処理した特異値で行列Xを定義
end

function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);%shiftdimで三次元が1次元方向に持ってこられている（多分）,reshapeで２次元にしている。3×65536にしていると思われる。
end

function psnr = PSNR(Xfull,Xrecover,maxP)

% Written by Canyi Lu (canyilu@gmail.com)

Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
[n1,n2,n3] = size(Xrecover);
MSE = norm(Xfull(:)-Xrecover(:))^2/(n1*n2*n3);
psnr = 10*log10(maxP^2/MSE);
end

function [X] = Fold(X, dim, i)
dim = circshift(dim, [1-i, 1-i]);
X = shiftdim(reshape(X, dim), length(dim)+1-i);
end

function [SigmaY,C_df, SURE] = calcSURE(Y, lambda, nsig)
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