function z=chasi_v6_tr38901_spatial_normal(stream,positionsMeters, ...
    correlationDistanceMeters)
% TR 38.901 v17.1.0 7.6.3.1: exponential 2-D field at supplied locations.
% Each column is one independently seeded unit-normal random field.
% Cross-LSP correlation is applied subsequently by chasi_v6_tr38901_lsp.
T=size(positionsMeters,1);P=numel(correlationDistanceMeters);
assert(T>=1 && size(positionsMeters,2)==2 && ...
    all(isfinite(positionsMeters(:))) && ...
    all(isfinite(correlationDistanceMeters(:)) & ...
        correlationDistanceMeters(:)>0));
deltaX=positionsMeters(:,1)-positionsMeters(:,1).';
deltaY=positionsMeters(:,2)-positionsMeters(:,2).';
distance=hypot(deltaX,deltaY);
z=zeros(T,P);uniqueDistance=unique(correlationDistanceMeters(:).');
for d=uniqueDistance
 columns=find(correlationDistanceMeters==d);
 covariance=exp(-distance/d);
 [factor,flag]=chol(covariance+1e-10*eye(T),'lower');
 assert(flag==0,'Spatial covariance is not positive definite');
 z(:,columns)=factor*randn(stream,T,numel(columns));
end
end
