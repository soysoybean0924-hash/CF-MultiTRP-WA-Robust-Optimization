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
zBound=sum(amplitude,'all')/sqrt(losField.features);
zFirstBound=sum(amplitude.*waveSpeed,'all')/sqrt(losField.features);
zSecondBound=sum(amplitude.*waveSpeed.^2,'all')/sqrt(losField.features);
bound=zFirstBound/sqrt(2*pi)+.3*speed;
% This also bounds the curvature of the outdoor LOS probability for d>18 m.
curvatureBound=(zSecondBound+zBound*zFirstBound^2)/sqrt(2*pi)+speed^2;
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
    function value=distanceAt(t)
        value=norm(initialXY+t*v-trp);
    end
    function value=marginSlope(t)
        xy=initialXY+t*v;d=norm(xy-trp);
        theta=losField.kx*xy(1)+losField.ky*xy(2);
        z=sum(losField.cosWeight.*cos(theta)+ ...
            losField.sinWeight.*sin(theta),'all')/sqrt(losField.features);
        dz=sum((-losField.cosWeight.*sin(theta)+ ...
            losField.sinWeight.*cos(theta)).* ...
            (losField.kx*v(1)+losField.ky*v(2)),'all')/ ...
            sqrt(losField.features);
        scale=36;if ~strcmp(field.scenario,'UMi'),scale=63;end
        e=exp(-d/scale);
        dp=-18/d^2*(1-e)-(1-18/d)*e/scale;
        dd=dot(xy-trp,v)/d;
        value=exp(-z^2/2)*dz/sqrt(2*pi)-dp*dd;
    end
    function visit(a,b,ma,mb,depth)
        width=b-a;
        if abs(ma)>bound*width,return;end
        % A constant derivative sign certifies either no crossing or
        % exactly one crossing, even when the global bound is loose.
        smooth=min(distanceAt(a),distanceAt(b))-speed*width>18;
        if smooth && abs(marginSlope(a))>curvatureBound*width
            if (ma<0)==(mb<0),return;end
            low=a;high=b;lowState=ma<0;
            for bisect=1:48
                middle=(low+high)/2;
                if (margin(middle)<0)==lowState
                    low=middle;
                else
                    high=middle;
                end
            end
            lastEvent=high;count=count+1;return
        end
        assert(depth<48 && width>eps(max(1,abs(b))), ...
            'LOS event requires higher precision; inspect field and trajectory');
        middle=(a+b)/2;mm=margin(middle);
        visit(a,middle,ma,mm,depth+1);
        visit(middle,b,mm,mb,depth+1);
    end
end
