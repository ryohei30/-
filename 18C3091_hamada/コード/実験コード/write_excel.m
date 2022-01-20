filename = 'study_data.xlsx';
MIN_SURE_dft = min(SURE_dft_list);
MIN_SURE_dct = min(SURE_dct_list);
MIN_SURE_pca = min(SURE_pca_list);
MIN_MSE_dft = min(MSE_dft_list);
MIN_MSE_dct = min(MSE_dct_list);
MIN_MSE_pca = min(MSE_pca_list);
MAX_PSNR_dft = max(PSNR_dft_list);
MAX_PSNR_dct = max(PSNR_dct_list);
MAX_PSNR_pca = max(PSNR_pca_list);

%writecell(SUREs,filename,'Sheet',4,'Range','B2');
%writecell(psnrs,filename,'Sheet',4,'Range','G2');
%writecell(MSEs,filename,'Sheet',4,'Range','K2'); 
writematrix(MIN_SURE_dft,filename,'Sheet',4,'Range','B66');
writematrix(MIN_SURE_dct,filename,'Sheet',4,'Range','C66');
writematrix(MIN_SURE_pca,filename,'Sheet',4,'Range','D66');
writematrix(MIN_MSE_dft,filename,'Sheet',4,'Range','G66');
writematrix(MIN_MSE_dct,filename,'Sheet',4,'Range','H66');
writematrix(MIN_MSE_pca,filename,'Sheet',4,'Range','I66');

writematrix(MAX_PSNR_dft,filename,'Sheet',4,'Range','L66');
writematrix(MAX_PSNR_dct,filename,'Sheet',4,'Range','M66');
writematrix(MAX_PSNR_pca,filename,'Sheet',4,'Range','N66');

