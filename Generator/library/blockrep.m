function b = blockrep(a,n)

% Block replicate the pixels of image cube a. Each pixel of each slice is replicated into a block of nxn pixels

b = zeros(size(a,1)*n,size(a,2)*n,size(a,3));
for k=1:n; b(1:n:end,k:n:end,:) = a; end
for k=2:n; b(k:n:end,:,:)=b(1:n:end,:,:); end

