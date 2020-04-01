function out=block_avg(in,pix)


out=zeros(size(in,1)./pix,size(in,2)./pix,size(in,3));
for k=1:size(in,3)
    rr=1;
    cc=1;
    for i=1:pix:size(in,1)-(pix-1)
        for j=1:pix:size(in,2)-(pix-1)
            out(rr,cc,k)=sum(sum(in(i:i+(pix-1),j:j+(pix-1),k)));
            cc=cc+1;
        end
        rr=rr+1;
        cc=1;
    end
end
out=out./pix^2;