function out=chasi_v6_tr38901_lsp_field(mode,arg1,arg2,arg3,arg4)
% Joint stationary 2-D Gaussian spectral field for correlated LSPs.
% Marginal spatial spectra are exponential with their own distances;
% the zero-separation cross covariance is the 7.5-6 LSP matrix.
if strcmp(mode,'create')
 stream=arg1;distance=double(arg2(:).');C=double(arg3);
 features=arg4;P=numel(distance);
 assert(P==7 && isequal(size(C),[P,P]) && ...
     all(isfinite(distance) & distance>0) && ...
     features>=128 && features==fix(features));
 selected=randi(stream,P,features,1);
 u=min(1-eps,max(eps,rand(stream,features,1)));
 selectedDistance=distance(selected); selectedDistance=selectedDistance(:);
 radius=sqrt((1-u).^(-2)-1)./selectedDistance;
 azimuth=2*pi*rand(stream,features,1);
 kx=radius.*cos(azimuth);ky=radius.*sin(azimuth);
 spectral=zeros(features,P);
 for j=1:P
  spectral(:,j)=distance(j)^2./ ...
      (1+(radius*distance(j)).^2).^(3/2);
 end
 proposal=mean(spectral,2);
 weight=sqrt(spectral./proposal);
 overlap=(weight.'*weight)/features;
 coupling=C./overlap;coupling=(coupling+coupling.')/2;
 [basis,eigenvalue]=eig(coupling);
 spectrum=max(1e-10,real(diag(eigenvalue)));
 coupling=basis*diag(spectrum)*basis.';
 scale=sqrt(diag(coupling).*diag(overlap));
 coupling=coupling./(scale*scale.');
 coupling=(coupling+coupling.')/2;
 factor=chol(coupling+1e-12*eye(P),'lower');
 realizedCross=coupling.*overlap;
 out=struct('kx',kx,'ky',ky,'weight',weight, ...
     'cosCoefficient',randn(stream,features,P)*factor.', ...
     'sinCoefficient',randn(stream,features,P)*factor.', ...
     'crossCorrelation',C,'realizedCrossCorrelation',realizedCross, ...
     'crossCorrelationMismatch',max(abs(realizedCross-C),[],'all'), ...
     'distanceMeters',distance, ...
     'features',features);
elseif strcmp(mode,'evaluate')
 field=arg1;xy=double(arg2(:).');
 assert(numel(xy)==2 && all(isfinite(xy)));
 phase=field.kx*xy(1)+field.ky*xy(2);
 out=(sum(field.weight.*(field.cosCoefficient.*cos(phase)+ ...
     field.sinCoefficient.*sin(phase)),1).'/sqrt(field.features));
else
 error('Unknown joint LSP field mode');
end
end
