%% 基于WSNM的RNIPT算法
tic
clc;
clear;
close all;

% setup parameters
lambdaL = 5;
C=5;
p=0.8;
L=4;%L

strDir='data\';%input data

%%
k=20; %total image number
for i=1:k
    fprintf('%d/%d: %s\n', 120, i);
    picname=[strDir  num2str(i),'.bmp'];
    I=imread(picname);
    [~, ~, ch]=size(I);
    if ch==3
        I=rgb2gray(I);
    end
    D(:,:,i)=I;
end
tenD=double(D);
[n1,n2,n3]=size(tenD);
n_1=max(n1,n2);%n(1)
n_2=min(n1,n2);%n(2)
patch_frames=L;% temporal slide parameter
patch_num=n3/patch_frames;
for l=1:patch_num
    l
    for i=1:patch_frames
        temp(:,:,i)=tenD(:,:,patch_frames*(l-1)+i);
    end           
        T=C*sqrt(n1*n2);
        lambda =lambdaL / sqrt(min(n_1*patch_frames));
        mu = 1e-2;
        [tenB, tenT,change] = WSNMSTIPT(temp, lambda, mu,T,p);  
 for i=1:patch_frames
      tarImg=tenT(:,:,i);
      backImg=tenB(:,:,i);
        a=uint8(tarImg);        
        figure;
        imshow(a, []);
end 
end
toc