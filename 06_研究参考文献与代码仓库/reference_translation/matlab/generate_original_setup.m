function scene=generate_original_setup(seed,realizations,apFixed)
% Independent MATLAB generation of one original nonideal-CSI setup.
% Implements original main_DNN.py:26-298 and APLocation_Generation.py.
% Same parameters/distribution, not identical random samples across languages.
if ~isempty(seed), rng(seed,'twister'); end
L=16; K=20; M=4; tau_p=10; p=100; Pmax=1000;
squareLength=1000; gridSide=4; spacing=250;
if nargin>=3 && ~isempty(apFixed)
    ap=apFixed;
else
    [gx,gy]=meshgrid(linspace(spacing/2,squareLength-spacing/2,gridSide));
    ap=complex(gx(:),gy(:));
    ap=ap+(rand(L,1)-0.5)*spacing/2+1i*(rand(L,1)-0.5)*spacing/2;
end
ue=complex(squareLength*rand(K,1),squareLength*rand(K,1));
noiseDbm=-174+10*log10(20e6)+7;
beta=zeros(K,L); CorrR=complex(zeros(M,M,K,L));
for k=1:K
    for l=1:L
        dx=real(ue(k)-ap(l)); dy=imag(ue(k)-ap(l));
        if abs(dx)>squareLength/2
            dx=(squareLength-abs(dx))*sign(-dx);
        end
        if abs(dy)>squareLength/2
            dy=(squareLength-abs(dy))*sign(-dy);
        end
        distance=sqrt(10^2+dx^2+dy^2);
        gainDb=-30.5-36.7*log10(distance)-noiseDbm;
        beta(k,l)=10^(gainDb/10);
        theta=angle(dx+1i*dy); asd=10*pi/180;
        column=0:M-1;
        firstRow=exp(1i*pi*sin(theta)*column).* ...
            exp(-(asd^2/2)*(pi*cos(theta)*column).^2);
        CorrR(:,:,k,l)=beta(k,l)*toeplitz(firstRow);
    end
end
pilot=-ones(K,1);
pilot(1:tau_p)=randperm(tau_p).';
for k=tau_p+1:K
    [~,l]=max(beta(k,:));
    interference=zeros(tau_p,1);
    for t=1:tau_p
        interference(t)=sum(beta(pilot(1:k-1)==t,l));
    end
    [~,pilot(k)]=min(interference);
end
z=(randn(M,realizations,K,L)+1i*randn(M,realizations,K,L))/sqrt(2);
H=complex(zeros(M,realizations,K,L));
for l=1:L
    for k=1:K
        H(:,:,k,l)=sqrtm(CorrR(:,:,k,l))*z(:,:,k,l);
    end
end
Np=(randn(M,realizations,L,tau_p)+ ...
    1i*randn(M,realizations,L,tau_p))/sqrt(2);
Hhat=pilot_mmse_from_channels(H,CorrR,Np,pilot-1,p,tau_p);
scene.H=H; scene.Hhat=Hhat; scene.CorrR=CorrR;
scene.ap=ap; scene.ue=ue; scene.beta=beta; scene.pilot=pilot;
scene.Np=Np; scene.p=p; scene.Pmax=Pmax;
scene.prelog=(200-tau_p)/200;
end
