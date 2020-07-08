function [ x ] = SIRT( x,y,A,alpha,iter )

ImageSize = size(x,1);
x = x(:);
y = y(:);


x_pre = zeros(ImageSize*ImageSize,1);
for ind=1:iter       
    ind
    y_pre=A*x_pre;
    dB=y-y_pre;
    x=x_pre+alpha.*(A'*dB);
    NegInd=find(x<0);
    x(NegInd)=0;
%     if mod(ind,20)==1 && ind <700
%         x=reshape(x,ImageSize,ImageSize);
%         x=medfilt2(x,[3,3]);
%         x=reshape(x,ImageSize*ImageSize,1);
%     end
    x_pre=x;
end
x = reshape(x,[256,256]);
end

