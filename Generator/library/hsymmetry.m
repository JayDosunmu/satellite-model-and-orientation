function cout = hsymmetry( cin )

% use Hermitian symmetry ( for real signal ) to compute the full spectrum
% from the half spectrum - assumes DC is centered !

mdim = size(cin,2);
icen=mdim/2 + 1;

cout = cin;


%** conjugate the zero frequencies change
cout(icen,icen+1:mdim) = conj(cout(icen,(icen-1):-1:2));
cout(icen+1:mdim,icen) = conj(cout((icen-1):-1:2,icen));


% stripe the high frequencies 
cout(icen+1:mdim,1) = conj(cout(icen-1:-1:2,1));


%** fill in other 2 quadrants 
cout((icen+1):mdim,(icen+1):mdim) = conj(cout((icen-1):-1:2,(icen-1):-1:2));
cout((icen+1):1:mdim,(icen-1):-1:2) = conj(cout((icen-1):-1:2,(icen+1):1:mdim));

cout  = fftshift(cout);

