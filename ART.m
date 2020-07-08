%% SB Recon
% Proj: original projection
% Img: Recon image
% B: Projection Matrix
% 
% By Mianyi Chen, 20180907
%% loading data and sampling projection
% load('ProjHead.mat')
% load('systemMatrix.mat')



[DetNum,angles]=size(Proj);
sampling=5;
Proj=Proj(:,1:sampling:end);

a=1;
for i=1:sampling:angles
    i
    B1((DetNum*(a-1)+1:DetNum*a),:)=B(DetNum*(i-1)+1:DetNum*i,:);
    a=a+1;
end

% B1 = SystemMatrix;
% Proj = y;

B1 = B;


%% recon:ART
ImgSize=256;
W=sum(B1.^2,2);
BT = B1';
Img=zeros(ImgSize,ImgSize);


for i=1:500
    i
    P=(Proj(:)-B1*Img(:));
%     P1=sparse(P');
    Img(:)=Img(:)+(BT*(P./(W).*0.001));
    Img(Img(:)<0)=0;
    
    if mod(i,2000) == 1
        save(strcat('result',num2str(i)),'Img'); 
    end
end

for i=1:ImgSize
    for j=1:ImgSize
        if (i-ImgSize/2)^2+(j-ImgSize/2)^2>(ImgSize/2-13)^2
            Img(i,j)=0;
        end
    end
end