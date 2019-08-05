function N = kmer_norm_weights(W)

bgD = 0.2;

% calculate positional normalization factors
x1 = W;
x1(W<0) = 0;
bg1 = sum(abs(x1),1);
bg1 = bg1./median(bg1);
bg1(abs(bg1-1)<bgD) = 1;
bg1(bg1<1) = 1;

x2 = W;
x2(W>0) = 0;
bg2 = sum(abs(x2),1);
bg2 = bg2./median(bg2);
bg2(abs(bg2-1)<bgD) = 1;
bg2(bg2<1) = 1;

b1 = repmat(bg1,size(W,1),1);
b2 = repmat(bg2,size(W,1),1);
BG = zeros(size(W));
BG(W>0) = b1(W>0);
BG(W<0) = b2(W<0);

N = W./BG;
N(isnan(N)) = 0;
