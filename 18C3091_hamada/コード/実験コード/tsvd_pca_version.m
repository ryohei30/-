%function I = tsvd_pca_version(image)
% t-SVD�ɂ��K�E�V�A���e���\���f�m�C�Y
% �^�������̃V�[�h
%rng(1)

% �p�����[�^
%  ���e���\���̃p�����[�^
%N = 256;
%M = 256;
%K = 3;

% �ϑ��e���\���̃p�����[�^
%noise_sigma = 0.2;

%   t-SVD�̃p�����[�^
%lambda = 0.1;

%�œK��lambda��T���ꍇ
%lambda_list = [];
%psnr_list = [];
%for i=0:0.1:5
%    lambda = 0.1+i;
%    lambda_list(end+1) = lambda;


% �K���Ȓ჉���N�e���\�������
%a = randn(N,1,1);
%b = randn(1,M,1);
%c = randn(1,1,K);
%X = a.*b.*c;

%�摜���g���ꍇ
%pic_name = './Mandrill.bmp';
%pic_name = image;
%X = double(imread(pic_name));
%X = image;
%X = X/255;

%Y = X + randn(size(X))*noise_sigma;

function [SigmaY,SURE_pca,psnr_pca,MSE_pca] = tsvd_pca_version(Y,K,noise_sigma,X)
lambda =4.475;

%��R��������PCA
[n1,n2,n3] = size(Y);
UfY = Unfold(Y, [n1,n2,n3], 3);%�ϑ���unfold
cUfY = UfY-mean(UfY,2);
C = cUfY*cUfY'/(n1*n2);
%UfXh = Unfold(X, [n1,n2,n3], 3);%���e���\����unfold
%cUfXh = UfXh-mean(UfXh,2); %���S��
%C = cUfXh*cUfXh'/(n1*n2); %�����U�s����l��
[Eig,~] = eig(C); %�ŗL�x�N�g�����擾
PY = Eig'*UfY; %�ŗL�x�N�g���s�񂩂���e���\��unfold
Yb = Fold(PY,[n1,n2,n3],3); %PCA�����e���\����Yb�Ƃ���
Ybd = zeros(size(Yb)); %Yb�̓����T�C�Y�̋�̃e���\�����`�i�����ɏ�����̒l������j

SURE_list = [];
    
% ���ْl�������ă\�t�g臒l����
for i = 1:K
    Ybd(:,:,i) = SVDandThresh(Yb(:,:,i), lambda);
    %SURE�̌v�Z
    [SigmaY,~, SURE] = calcSURE(Ybd(:,:,i), lambda, noise_sigma);
    SURE_list(end+1) = SURE;
end
 UfX = Unfold(Ybd, [n1,n2,n3], 3); %臒l���������e���\����unfold�������̂�Yd�Ƃ���
 iPX = Eig*UfX; %unfold����Ybd�ɌŗL�x�N�g���s��������Ă�ipca���Ă�H�j
 
 Yd = Fold(iPX, [n1,n2,n3], 3); %�e���\���ɖ߂�

% MSE���͂���
MSE_pca = mean(sum((X - Yd).^2),'all'); 

%PSNR��SURE���Z�o
%psnr_list(end+1) = PSNR(X,Yd,1.0);
psnr_pca = PSNR(X,Yd,1.0);
SURE_pca = sum(SURE_list);
%end

%�Ԃ��l
%[M,I]=max(psnr_list);

% ���ʕ\��
%fprintf("MSE = %f\n", MSE) %MSE�̌���
%fprintf("psnr = %f\n", psnr) %PSNR�̌���

%result2 = max(psnr_list);
%fprintf("result2 = %f\n", result2)
%max(lambda_list)


%�v���b�g
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
    [U,S,V] = svd(Y); %���ْl���������s
    S=max(S-lambda, 0); %S��臒l�������AS���Ē�`�i臒l��lambda)
    X = U*S*V'; %臒l�����������ْl�ōs��X���`
end

function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);%shiftdim�ŎO������1���������Ɏ����Ă����Ă���i�����j,reshape�łQ�����ɂ��Ă���B3�~65536�ɂ��Ă���Ǝv����B
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
    %idx_w = (sY-lambda >0);%�������l������0�ɂ���Ȃ��͈́if_i�̔����l������̂͂��������j
    term1 = abs(m-n) * sum(max((1- lambda./sY),0)); %Eq. (9)��1��
    term2 = sum(sY > lambda); %Eq. (9)��2�� �т�����
    numer = repmat(sY.*max((sY-lambda),0), [1,minnm]); %Eq. (9)��3���̕��q
    denom = repmat(sY.^2, [1,minnm]) - repmat((sY.^2)', [minnm,1]); %Eq. (9)��3���̕���
    DD = numer./denom; % i=j�̂Ƃ�NaN���o�邪�c
    DD( isnan(DD) | isinf(DD) | abs(DD) > 1E6 )  =   0;%�����ŏ����̂Ŗ��Ȃ������Ǝד��Ȃ̂ł����Ƃ�����ق�������
    term3 = 2*sum(DD(:)); %Eq. (5)��3��
    C_df = term1+term2+term3;
    SURE = - m*n*nsig^2 + sum(min(lambda^2,sY.^2)) + 2*nsig^2*C_df;
end

end