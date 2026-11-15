function panel=chasi_tr38901_array_panel(bearingRad,azimuthRad,zenithRad)
% TR 38.901 7.3 panel geometry and single-element field pattern.
% Eight horizontal by four vertical locations, each with +/-45 degree
% polarization, give 64 ports. Coordinates are in wavelengths.
assert(isscalar(bearingRad) && isfinite(bearingRad));
side=[-sin(bearingRad);cos(bearingRad);0];up=[0;0;1];
position=zeros(3,64);slant=zeros(1,64);j=0;
for row=1:4
 for col=1:8
  for pol=1:2
   j=j+1;
   position(:,j)=.5*((col-4.5)*side+(row-2.5)*up);
   slant(j)=(-1)^pol*pi/4;
  end
 end
end
relativeAz=atan2(sin(azimuthRad-bearingRad), ...
    cos(azimuthRad-bearingRad))*180/pi;
zenithDeg=zenithRad*180/pi;
horizontal=min(12*(relativeAz/65)^2,30);
vertical=min(12*((zenithDeg-90)/65)^2,30);
attenuationDb=min(horizontal+vertical,30);
fieldAmplitude=10^(-attenuationDb/20);
panel=struct('positionWavelength',position, ...
    'polarizationSlantRad',slant, ...
    'fieldTheta',fieldAmplitude*cos(slant), ...
    'fieldPhi',fieldAmplitude*sin(slant), ...
    'elementPowerGain',fieldAmplitude^2, ...
    'attenuationDb',attenuationDb, ...
    'reference','TR 38.901 7.3 single-element 65-degree, 30-dB cap and 8x4x2 panel; array response not yet used by the v6 channel sampler');
end
