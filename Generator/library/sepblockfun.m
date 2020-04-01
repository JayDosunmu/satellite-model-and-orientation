function X=sepblockfun(X,blockdims,fun)

if issparse(X)&& exist('ndSparse','class') && ~isa(X,'ndSparse')
     X=ndSparse(X);
end
if ischar(fun)
  switch fun
    case 'max'
         fun=@(b,d) max(b,[],d);
    case 'min'
         fun=@(b,d) min(b,[],d);
    case 'sum'
         fun=@sum;
    case 'mean'
         fun=@mean;
    case 'prod'
        fun=@prod;
    otherwise
     error 'Unrecognized fun() selection'
  end
end

nn=max(length(blockdims),ndims(X));
blockdims(end+1:nn)=1;

[sz{1:nn}]=size(X); %M is the original array
sz=[sz{:}];

idx=~isfinite(blockdims);
blockdims(idx)=sz(idx);

newdims=(sz./blockdims);
args=(num2cell([blockdims;newdims]));
X=reshape(X,args{:});
for ii=1:nn
 X=fun(X,2*ii-1);
end
X=reshape(X,newdims);