function F = plot_Msaccade(tarpos,spike,title_define)
% % % plot the response in memory saccade task
% % % the input argument spike is the mean firing rate in a certain time
% % % window, not something like PSTH


% % % do a 2-D convolution at first 
x_range = max(abs(tarpos(1,:)))+1;
y_range = max(abs(tarpos(2,:)))+1;
mu=[0, 0];
Sigma=[0.5 0;0 0.5];
[X,Y]=meshgrid(-x_range:0.1:x_range,-y_range:0.1:y_range);
kernal=mvnpdf([X(:) Y(:)],mu, Sigma);
kernal = kernal/max(max(kernal));
kernal=reshape(kernal,size(X));
pluse = zeros(size(X));
response = zeros([length(tarpos),size(X)]);
weight_kernal = response;
weight_mu = tarpos';
x_pos = round(vpa((tarpos(1,:)+x_range)/0.1,2))+1;
y_pos = round(vpa((tarpos(2,:)+y_range)/0.1,2))+1;
for i  = 1:size(tarpos,2)
    weight_kernal(i,:,:) = reshape(mvnpdf([X(:) Y(:)], weight_mu(i,:),Sigma),size(X));
    pluse(:) = 0;
    pluse(y_pos(i),x_pos(i)) = spike(i);
    response(i,:,:) = conv2(pluse,kernal,'same');
end

reponse_final = squeeze(sum(weight_kernal.*response,1))./squeeze(sum(weight_kernal,1));
% % % % % the response of the position we didn't test is set to -2
% % % reponse_final(squeeze(sum(weight_kernal,1))<=0.0001) = -2; 
%%
F = figure('Paperunits','Normalized','Paperposition',[0.1 0.1 1 0.5],...
    'units','Normalized','position',[0.1 0.1 0.815 0.8]);
set(gcf,'Position',get(gcf,'Position').*[1 1 1.3 1])

subplot(2,3,[1 2 4 5])
surf(X,Y,reponse_final);axis tight;
grid on;shading flat;hidden off

title(['response of memory saccade--',title_define])
subplot(2,3,3)
surf(X,Y,reponse_final),view(2),axis tight,title('XOY')
axis tight;grid on;shading flat;hidden off
subplot(2,3,6)
surf(X,Y,reponse_final),view([0 0]),axis tight,title('XOZ')
axis tight;grid on;shading flat;hidden off

end
