% t-SVDによるガウシアンテンソルデノイズ
% 疑似乱数のシード
rng('default')


K = 3;

% 観測テンソルのパラメータ
noise_sigma = 0.2;



%原テンソル
%pic_name = './Airplane.bmp';
%X = double(imread(pic_name));
%X = X/255;

%原テンソル12枚全部
bmpFiles = dir('*.bmp');
numfiles = length(bmpFiles);
mydata = cell(1,numfiles);

SUREs = cell(numfiles,1);
psnrs = cell(numfiles,1);
MSEs = cell(numfiles,1);
SUREs{1} = ['dft,dct,pca'];
psnrs{1} = ['dft,dct,pca'];
MSEs{1} = ['dft,dct,pca'];

%観測テンソル
%Y = X + randn(size(X))*noise_sigma;

%SUREとPSNRを求める
%[SURE_dft,psnr_dft] = tsvd_dft(Y,K,noise_sigma,X);
%[SURE_dct,psnr_dct] = tsvd_dct_version(Y,K,noise_sigma,X);
%[SURE_pca,psnr_pca] = tsvd_pca_version(Y,K,noise_sigma,X);

%SUREとPSNRを求める12枚全部
for k = 1: numfiles
    X = double(imread(bmpFiles(k).name));
    X = X/255;
    Y = X + randn(size(X))*noise_sigma;
   [SigmaY_dft,SURE_dft,psnr_dft,MSE_dft] = tsvd_dft_version(Y,K,noise_sigma,X);
   [SURE_dct,psnr_dct,MSE_dct] = tsvd_dct_version(Y,K,noise_sigma,X);
   [SigmaY_pca,SURE_pca,psnr_pca,MSE_pca] = tsvd_pca_version(Y,K,noise_sigma,X);

    %SURE_dft = cast(SURE_dft,'single');
    %SURE_dct = cast(SURE_dct,'single');
    %SURE_pca = cast(SURE_pca,'single');
    %psnr_dft = cast(psnr_dft,'single');
    %psnr_dct = cast(psnr_dct,'single');
    %psnr_pca = cast(psnr_pca,'single');
    
    bmpFiles(k).name
    SUREs{k+1} = [SURE_dft,SURE_dct,SURE_pca];
    psnrs{k+1} = [psnr_dft,psnr_dct,psnr_pca];
    MSEs{k+1} = [MSE_dft,MSE_dct,MSE_pca];

    %一番低いSUREと一番高いPSNR
    %{
    if SURE_dct > SURE_dft && SURE_pca > SURE_dft
        SURE_dft
    elseif SURE_dft > SURE_dct && SURE_pca > SURE_dct
        SURE_dct
    elseif SURE_dft > SURE_pca && SURE_dct > SURE_pca
        SURE_pca
    end

    if psnr_dct < psnr_dft && psnr_pca < psnr_dft
        psnr_dft
    elseif psnr_dft < psnr_dct && psnr_pca < psnr_dct
        psnr_dct
    elseif psnr_dft < psnr_pca && psnr_dct < psnr_pca
        psnr_pca
    end
    %}



end