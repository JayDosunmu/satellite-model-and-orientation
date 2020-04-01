function out=block_rep(in,ssp)

out=zeros(size(in,1)*ssp,size(in,2)*ssp,size(in,3));

for k=1:size(in,3)
    ii=1;
    for i=1:size(in,1)
        jj=1;
        for j=1:size(in,2)
            out(ii:ii+(ssp-1),jj:jj+(ssp-1),k)=in(i,j,k);
            jj=jj+ssp;
        end
        ii=ii+ssp;
    end
end

            