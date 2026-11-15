function out=chasi_v6_tr38901_spectral_field(mode,arg1,arg2,arg3,arg4)
% Continuous 2-D exponential-correlation field using frozen Fourier modes.
% The finite-mode approximation is explicit; no new random draws at read.
if strcmp(mode,'create')
 stream=arg1;distanceMeters=arg2;variables=arg3;features=arg4;
 assert(distanceMeters>0 && variables>=1 && features>=32);
 u=min(1-eps,max(eps,rand(stream,features,variables)));
 radius=sqrt((1-u).^(-2)-1)/distanceMeters;
 angle=2*pi*rand(stream,features,variables);
 out=struct('kx',radius.*cos(angle),'ky',radius.*sin(angle), ...
     'cosWeight',randn(stream,features,variables), ...
     'sinWeight',randn(stream,features,variables), ...
     'distanceMeters',distanceMeters,'features',features);
elseif strcmp(mode,'evaluate')
 field=arg1;xy=double(arg2(:).');
 assert(numel(xy)==2 && all(isfinite(xy)));
 argument=field.kx*xy(1)+field.ky*xy(2);
 out=sum(field.cosWeight.*cos(argument)+ ...
     field.sinWeight.*sin(argument),1).'/sqrt(field.features);
else
 error('Unknown spatial field mode');
end
end
