function [lastEvent,count]=chasi_v6_los_transition_events( ...
    field,initialXY,velocity,trp,elapsed)
% Enumerate every hard LOS boundary along the frozen straight trajectory.
% A Lipschitz bound certifies intervals without zero crossings; unresolved
% intervals are bisected, including intervals with equal endpoint states.
initialXY=double(initialXY(:).');velocity=double(velocity(:).');
trp=double(trp(:).');assert(elapsed>=0 && isfinite(elapsed));
v=velocity(1:2);speed=norm(v);
lastEvent=0;count=0;
if elapsed==0 || speed==0,return;end
losField=field.los;
waveSpeed=abs(losField.kx*v(1)+losField.ky*v(2));
amplitude=hypot(losField.cosWeight,losField.sinWeight);
% Gaussian CDF derivative <= 1/sqrt(2*pi); outdoor LOS probability
% derivative for d>=10 m is bounded by 0.3 per meter.
bound=sum(amplitude.*waveSpeed,'all')/ ...
    sqrt(losField.features*2*pi)+.3*speed;
fa=margin(0);fb=margin(elapsed);
visit(0,elapsed,fa,fb,0);
    function value=margin(t)
        xy=initialXY+t*v;d=max(10,norm(xy-trp));
        if strcmp(field.scenario,'UMi'),scale=36;else,scale=63;end
        probability=min(1,18/d+(1-18/d)*exp(-d/scale));
        normal=chasi_v6_tr38901_spectral_field( ...
            'evaluate',losField,xy);
        value=.5*(1+erf(normal/sqrt(2)))-probability;
    end
    function visit(a,b,ma,mb,depth)
        if abs(ma)>bound*(b-a),return;end
        if depth>=48 || speed*(b-a)<1e-9
            assert((ma<0)~=(mb<0), ...
                'LOS boundary cannot be resolved to 1 nm; inspect trajectory');
            lastEvent=b;count=count+1;return
        end
        middle=(a+b)/2;mm=margin(middle);
        visit(a,middle,ma,mm,depth+1);
        visit(middle,b,mm,mb,depth+1);
    end
end
